// C++ dependencies
#include <cmath>
#include <iostream>

// ROOT dependencies
#include "TGraph.h"

// X17 dependencies
#include "CircleFit3D.h"
#include "Field.h"
#include "Points.h"
#include "Vector.h"
#include "Utilities.h"
#include "X17Utilities.h"

namespace X17
{
    CircleFit3D* CircleFit3D::lastfit = nullptr;

    //// Public methods.

    CircleFit3D::CircleFit3D(const Vector& orig, const Vector& orient)
    {
        this->m_origin      = orig;
        this->m_orientation = orient;
        this->m_orientation.Normalize();

        this->m_fitter = lastfit->m_fitter; // Setting the fitter again every time would take much more time.
        lastfit = this;
    }

    void CircleFit3D::SetFitter(int parameters, bool print)
    {        
        m_fitter = TVirtualFitter::Fitter(nullptr,parameters);
        m_fitter->SetFCN(this->_Eval);

        if(!print)
        {
            double arg = -1;
            m_fitter->ExecuteCommand("SET PRINTOUT",&arg,1);
            m_fitter->ExecuteCommand("SET NOW", &arg ,1);
        }
    }

    void CircleFit3D::FitCircle3D(double max_iter, double toleration)
    {
        m_fitter->SetParameter(0,"length",m_length,0.01,-10,5);
        m_fitter->SetParameter(1,"alpha",m_alpha,0.001,-M_PI/2,(3/2)*M_PI);
        m_fitter->SetParameter(2,"radius",m_radius,0.01,10,50);
        m_fitter->SetParameter(3,"phi_max",m_phi_max,0.001,0.15,M_PI/1.5);

        double arglist[2] = {max_iter,toleration};    // Maximal number of iterations, step size (toleration).
        m_fitter->ExecuteCommand("MIGRAD",arglist,2); // Last parameter is the number of prints (verbosity).

        m_length  = m_fitter->GetParameter(0);
        m_l_err   = m_fitter->GetParError(0);
        m_alpha   = m_fitter->GetParameter(1);
        m_a_err   = m_fitter->GetParError(1);
        m_radius  = m_fitter->GetParameter(2);
        m_r_err   = m_fitter->GetParError(2);
        m_phi_max = m_fitter->GetParameter(3);
        m_phi_err = m_fitter->GetParError(3);

        _UpdateCurve();
    }
    
    void CircleFit3D::PrintFitParams()
    {
        std::cout << "\nCIRCLE FIT PARAMETERS:\n";
        std::cout << "Length:  " << m_length  << " +- " << m_l_err   << "\n";
        std::cout << "Alpha:   " << m_alpha   << " +- " << m_a_err   << "\n";
        std::cout << "Radius:  " << m_radius  << " +- " << m_r_err   << "\n";
        std::cout << "Phi_max: " << m_phi_max << " +- " << m_phi_err << "\n\n";
    }

    TGraph2D* CircleFit3D::GetGraph(double step, double dist)
    {
        TGraph2D* fit_graph = new TGraph2D();
        
        // 1st line
        double param = 0;
        while (param < m_length)
        {
            Vector cur_pos = _GetLinePoint(param,true);
            if(IsInSector(cur_pos,dist)) 
                fit_graph->AddPoint(cur_pos.x,cur_pos.y,cur_pos.z);
            param += step;
        }

        // circle
        param = 0;
        while (param/m_radius < m_phi_max)
        {
            Vector cur_pos = _GetCirclePoint(param / m_radius);
            if(IsInSector(cur_pos,dist)) 
                fit_graph->AddPoint(cur_pos.x,cur_pos.y,cur_pos.z);
            param += step;
        }

        // 2nd line
        param = 0;
        while (param < 20)
        {
            Vector cur_pos = _GetLinePoint(param,false);
            if(IsInSector(cur_pos,dist)) 
                fit_graph->AddPoint(cur_pos.x,cur_pos.y,cur_pos.z);
            param += step;
        }

        return fit_graph;
    }

    double CircleFit3D::GetEnergy(const Field<Vector>& magfield, bool middle)
    {
        using namespace constants;

        Vector bfield;
        if(middle) bfield = _GetMiddleField(magfield);
        else       bfield = _GetAvgField(magfield);

        double b_proj = m_normal*bfield;

        double betasq = 1 / (1 + pow((E0 / (c * (cm2m * m_radius) * b_proj)), 2));
        double Ekin   = E0 * (1 / sqrt(1 - betasq) - 1);

        return Ekin;
    }

    TGraph* CircleFit3D::GetEnergyGraph(const Field<Vector>& magfield, double step)
    {
        using namespace constants;

        TGraph* graph = new TGraph();
        double param = 0;

        while (param/m_radius < m_phi_max) 
        {
            Vector cur_point = _GetCirclePoint(param/m_radius);
            if(IsInSector(cur_point))
            {
                Vector bfield = magfield.GetField(cur_point);
                double b_proj = m_normal * bfield;

                double betasq = 1 / (1 + pow((E0 / (c * (cm2m * m_radius) * b_proj)), 2));
                double Ekin   = E0 * (1 / sqrt(1 - betasq) - 1);
                graph->AddPoint(param,Ekin);
            }
            param += step;
        }
        
        graph->SetTitle("Energy along the track;Parameter [cm];Energy [eV]");
        return graph;
    }





    //// Private methods.
    
    Vector CircleFit3D::_GetLinePoint(double param, bool first_line)
    {
        if (first_line) return m_origin  + param * m_orientation;
        else            return m_origin2 + param * m_orientation2;
    }

    Vector CircleFit3D::_GetCirclePoint(double varphi)
    {
        double cos_varphi = cos(varphi);
        double sin_varphi = sin(varphi);

        return m_originc + m_radius * Vector{  (1-cos_varphi) * (m_cos_alpha*m_cos_theta*m_cos_phi - m_sin_alpha*m_sin_phi) + sin_varphi*m_sin_theta*m_cos_phi,
                                               (1-cos_varphi) * (m_cos_alpha*m_cos_theta*m_sin_phi + m_sin_alpha*m_cos_phi) + sin_varphi*m_sin_theta*m_sin_phi,
                                              -(1-cos_varphi) *  m_cos_alpha*m_sin_theta + sin_varphi*m_cos_theta};
    }

    void CircleFit3D::_UpdateCurve()
    {
        m_cos_theta = m_orientation.z;
        m_sin_theta = sqrt(1 - m_cos_theta * m_cos_theta);

        m_cos_phi = m_orientation.x / m_sin_theta;
        m_sin_phi = m_orientation.y / m_sin_theta;

        m_cos_alpha = cos(m_alpha);
        m_sin_alpha = sin(m_alpha);

        m_originc = this->_GetLinePoint(m_length,true);

        m_normal = Vector {-m_sin_alpha*m_cos_theta*m_cos_phi - m_cos_alpha*m_sin_phi,
                           -m_sin_alpha*m_cos_theta*m_sin_phi + m_cos_alpha*m_cos_phi,
                            m_sin_alpha*m_sin_theta };
        
        m_center = m_originc + m_radius * Vector { m_cos_alpha*m_cos_theta*m_cos_phi - m_sin_alpha*m_sin_phi,
                                                   m_cos_alpha*m_cos_theta*m_sin_phi + m_sin_alpha*m_cos_phi,
                                                  -m_cos_alpha*m_sin_theta };
        
        m_origin2 = this->_GetCirclePoint(m_phi_max);

        double cos_varphi = cos(m_phi_max);
        double sin_varphi = sin(m_phi_max);

        m_orientation2 = Vector{ sin_varphi * (m_cos_alpha*m_cos_theta*m_cos_phi - m_sin_alpha*m_sin_phi) + cos_varphi*m_sin_theta*m_cos_phi,
                                 sin_varphi * (m_cos_alpha*m_cos_theta*m_sin_phi + m_sin_alpha*m_cos_phi) + cos_varphi*m_sin_theta*m_sin_phi,
                                -sin_varphi *  m_cos_alpha*m_sin_theta + cos_varphi*m_cos_theta };
    }

    double CircleFit3D::_LineSqDist(const RecoPoint& point, bool first_line)
    {
        Vector *orig,*orient;
        if (first_line) { orig = &m_origin;  orient = &m_orientation;  }
        else            { orig = &m_origin2; orient = &m_orientation2; }

        double t_close = *orient*(point.AsVector()-*orig)/orient->SqMagnitude();
        if (first_line  && (t_close > m_length)) t_close = m_length;
        if (!first_line && (t_close < 0))      t_close = 0;

        Vector line_vector = *orig+t_close*(*orient)-point.AsVector();

        return line_vector.SqMagnitude();
    }

    double CircleFit3D::_CircleSqDist(const RecoPoint& point)
    {
        Vector projection = (point.AsVector()-m_center)-((point.AsVector()-m_center)*m_normal)*m_normal;
        projection.Normalize();

        Vector circle_point   = m_radius * projection + m_center;
        Vector circle_vector  = circle_point - point.AsVector();
        Vector arc_beg_vector = m_originc - point.AsVector();
        Vector arc_end_vector = m_origin2 - point.AsVector();

        return find_min(circle_vector.SqMagnitude(),arc_beg_vector.SqMagnitude(),arc_end_vector.SqMagnitude());
    }
    
    double CircleFit3D::_SqDistance(const RecoPoint& point)
    {
        double sqdist1 = _LineSqDist(point,true);
        double sqdist2 = _CircleSqDist(point);
        double sqdist3 = _LineSqDist(point,false);

        return find_min(sqdist1,sqdist2,sqdist3);
    }

    double CircleFit3D::_SumSq()
    {
        double sum = 0;
        for (RecoPoint point : m_fit_data) sum += point.count*_SqDistance(point);
        return sum;
    }

    void CircleFit3D::_EvalSumSq(int& npar, double* gin, double& sumsq, double* par, int iflag)
    {
        m_length  = par[0];
        m_alpha   = par[1];
        m_radius  = par[2];
        m_phi_max = par[3];

        _UpdateCurve();

        sumsq = _SumSq();
    }

    void CircleFit3D::_Eval(int& npar, double* gin, double& sumsq, double* par, int iflag)
    {
        return lastfit->_EvalSumSq(npar,gin,sumsq,par,iflag);
    }

    Vector CircleFit3D::_GetAvgField(const Field<Vector>& magfield, double step)
    {
        int i = 0;
        double param = 0;
        Vector bfield = {0,0,0};

        while (param/m_radius < m_phi_max) 
        {
            Vector cur_point = _GetCirclePoint(param/m_radius);
            if(IsInSector(cur_point))
            {
                bfield += magfield.GetField(cur_point);
                i++;
            }
            param += step;
        }

        return bfield/i;
    }
    
    Vector CircleFit3D::_GetMiddleField(const Field<Vector>& magfield, double tolerance)
    {
        using namespace constants;

        double xmiddle = (xmax + xmin) / 2;
        double low  = 0;
        double high = m_phi_max;
        double mid  = (low + high) / 2;
        Vector vmid = _GetCirclePoint(mid);

        while ((abs(xmiddle - vmid.x) > tolerance) && (low * (1 + 1e-15) < high))
        {
            mid = (low + high) / 2;
            vmid  = _GetCirclePoint(mid);

            if (vmid.x < xmiddle) low  = mid;
            else                  high = mid;
        }
        
        return magfield.GetField(vmid);
    }
} // namespace X17