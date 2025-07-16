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
    //// Public methods.

    CircleFit3D::CircleFit3D(Vector orig, Vector orient)
    {
        SetOrigOrient(orig, orient);

        this->m_radius = 0;
        // this->m_radius.setCallback(coutCallback<double>("CircleFit3D::m_radius"));

        // Setting the fitter again every time would take much more time.
        // if (lastfit != nullptr) this->m_fitter = lastfit->m_fitter;
        lastfit = this;
    }

    void CircleFit3D::SetAlpha(bool electron)
    {
        if (electron) 
        {
            m_alpha = 0;
            // m_alpha_min = -M_PI/3;
            // m_alpha_max = M_PI/3;
        }
        else
        {
            m_alpha = M_PI;
            // m_alpha_min = 2*M_PI/3;
            // m_alpha_max = 4*M_PI/3;
        }
    }

    void CircleFit3D::SetFitter(int parameters, bool print)
    {
        if (!m_fitter) //(!m_fitter || m_fitter->GetNumberTotalParameters() != parameters)
        {
            m_fitter = TVirtualFitter::Fitter(nullptr,parameters);
        }
        
        if (!print)
        {
            double arg = -1;
            m_fitter->ExecuteCommand("SET PRINTOUT",&arg,1);
            m_fitter->ExecuteCommand("SET NOW", &arg ,1);
        }
        
        m_fitter->SetFCN(this->_Eval);
        m_fitter->SetParameter(0,"length",m_length,0.01,-5,7);
        m_fitter->SetParameter(1,"alpha",m_alpha,0.001,m_alpha_min-0.3,m_alpha_max+0.3);
        m_fitter->SetParameter(2,"radius",m_radius,0.01,5,5000);
        m_fitter->SetParameter(3,"phi_max",m_phi_max,0.001,0.03,M_PI/1.125);

        if (parameters == 6)
        {
            m_fitter->SetParameter(4,"m_theta",m_theta,0.001,M_PI/6,5*M_PI/6);
            m_fitter->SetParameter(5,"m_varphi",m_varphi,0.001,-M_PI/6,M_PI/6);
        }

        std::cout << "SetFitter: Fitter set with " << m_fitter->GetNumberFreeParameters() << " free parameters.\n";
    }

    void CircleFit3D::SetOrigOrient(Vector orig, Vector orient)
    {
        m_origin      = orig;
        m_orientation = orient;
        m_orientation.Normalize();

        m_theta = std::acos(m_orientation.z);
        m_varphi = std::acos(m_orientation.x / std::sin(m_theta));
    }

    void CircleFit3D::SetParameters(double length, double radius, double phi_max, bool electron)
    {
        SetAlpha(electron);

        m_length  = length;
        m_radius  = radius;
        m_phi_max = phi_max;

        _UpdateCurve();
    }

    void CircleFit3D::FitCircle3D(double max_iter, double toleration)
    {
        // Make sure that the fitter calls the correct function.
        lastfit = this;

        double count = 0;
        for (auto& p : m_fit_data) count += p.count;
        m_total_count = count;
        m_total_count.setCallback(coutCallback<int>("CircleFit3D::m_total_count"));

        double arglist[2] = {max_iter,toleration};    // Maximal number of iterations, step size (toleration).
        // m_fitter->ExecuteCommand("SIMPLEX", arglist, 2);
        int err = m_fitter->ExecuteCommand("MIGRAD",arglist,2); // Last parameter is the number of parameters in arglist.
        std::cout << "FitCircle3D: MIGRAD returned " << err << "\n";

        m_length  = m_fitter->GetParameter(0);
        m_l_err   = m_fitter->GetParError(0);
        m_alpha   = m_fitter->GetParameter(1);
        m_a_err   = m_fitter->GetParError(1);
        m_radius  = m_fitter->GetParameter(2);
        m_r_err   = m_fitter->GetParError(2);
        m_phi_max = m_fitter->GetParameter(3);
        m_phi_err = m_fitter->GetParError(3);

        if (m_fitter->GetNumberFreeParameters() == 6)
        {
            m_theta   = m_fitter->GetParameter(4);
            m_theta_err = m_fitter->GetParError(4);
            m_varphi  = m_fitter->GetParameter(5);
            m_varphi_err = m_fitter->GetParError(5);
        }

        _UpdateCurve();
    }
    
    void CircleFit3D::PrintFitParams() const
    {
        std::cout << "\nCIRCLE FIT PARAMETERS:\n";
        std::cout << "   Length:  " << m_length           << " +- " << m_l_err            << " cm\n";
        std::cout << "   Alpha:   " << m_alpha*180/M_PI   << " +- " << m_a_err*180/M_PI   << " deg\n";
        std::cout << "   Radius:  " << m_radius           << " +- " << m_r_err            << " cm\n";
        std::cout << "   Phi_max: " << m_phi_max*180/M_PI << " +- " << m_phi_err*180/M_PI << " deg\n";
        
        if (m_fitter->GetNumberFreeParameters() == 6)
        {
            std::cout << "   Theta:   " << m_theta*180/M_PI  << " +- " << m_theta_err*180/M_PI  << " deg\n";
            std::cout << "   Varphi:  " << m_varphi*180/M_PI << " +- " << m_varphi_err*180/M_PI << " deg\n";
        }

        std::cout << "Final sum of squares: " << _SumSq() << "\n\n";
    }

    TGraph2D* CircleFit3D::GetGraph(double step, double dist) const
    {
        TGraph2D* fit_graph = new TGraph2D();
        
        // 1st line
        double param = 0;
        while (param < m_length)
        {
            Vector cur_pos = _GetLinePoint(param,true);
            if(IsInTPC(cur_pos,dist)) 
                fit_graph->AddPoint(cur_pos.x,cur_pos.y,cur_pos.z);
            param += step;
        }

        // circle
        param = 0;
        while (param/m_radius < m_phi_max)
        {
            Vector cur_pos = _GetCirclePoint(param / m_radius);
            if(IsInTPC(cur_pos,dist)) 
                fit_graph->AddPoint(cur_pos.x,cur_pos.y,cur_pos.z);
            param += step;
        }

        // 2nd line
        param = 0;
        while (param < 20)
        {
            Vector cur_pos = _GetLinePoint(param,false);
            if(IsInTPC(cur_pos,dist)) 
                fit_graph->AddPoint(cur_pos.x,cur_pos.y,cur_pos.z);
            param += step;
        }

        return fit_graph;
    }

    double CircleFit3D::GetEnergy(const Field<Vector>& magfield, bool middle) const
    {
        using namespace constants;

        Vector bfield;
        if(middle) bfield = _GetMiddleField(magfield);
        else       bfield = _GetAvgField(magfield);

        double b_proj = m_normal*bfield;

        double betasq = 1 / (1 + pow((E0 / (c * (cm2m * m_radius) * b_proj)), 2));
        double Ekin   = E0 * (1 / std::sqrt(1 - betasq) - 1);

        if (!middle) Ekin /= 0.9263; // corection given by the test on RK tracks

        return Ekin;
    }

    TGraph* CircleFit3D::GetEnergyGraph(const Field<Vector>& magfield, double step) const
    {
        using namespace constants;

        TGraph* graph = new TGraph();
        double param = 0;

        while (param/m_radius < m_phi_max) 
        {
            Vector cur_point = _GetCirclePoint(param/m_radius);
            if(IsInTPC(cur_point))
            {
                Vector bfield = magfield.GetField(cur_point);
                double b_proj = m_normal * bfield;

                double betasq = 1 / (1 + pow((E0 / (c * (cm2m * m_radius) * b_proj)), 2));
                double Ekin   = E0 * (1 / std::sqrt(1 - betasq) - 1);
                graph->AddPoint(param,Ekin);
            }
            param += step;
        }
        
        graph->SetTitle("Energy along the track;Parameter [cm];Energy [eV]");
        return graph;
    }





    //// Private methods.
    
    Vector CircleFit3D::_GetLinePoint(double param, bool first_line) const
    {
        if (first_line) return m_origin  + param * m_orientation;
        else            return m_origin2 + param * m_orientation2;
    }

    Vector CircleFit3D::_GetCirclePoint(double varphi) const
    {
        double cos_varphi = std::cos(varphi);
        double sin_varphi = std::sin(varphi);

        return m_originc + m_radius * Vector{  (1-cos_varphi) * (m_cos_alpha*m_cos_theta*m_cos_phi - m_sin_alpha*m_sin_phi) + sin_varphi*m_sin_theta*m_cos_phi,
                                               (1-cos_varphi) * (m_cos_alpha*m_cos_theta*m_sin_phi + m_sin_alpha*m_cos_phi) + sin_varphi*m_sin_theta*m_sin_phi,
                                              -(1-cos_varphi) *  m_cos_alpha*m_sin_theta + sin_varphi*m_cos_theta};
    }

    void CircleFit3D::_UpdateCurve()
    {
        if (m_fitter && m_fitter->GetNumberFreeParameters() == 6)
            m_orientation = { std::sin(m_theta) * std::cos(m_varphi), std::sin(m_theta) * std::sin(m_varphi), std::cos(m_theta) };

        m_cos_theta = m_orientation.z;
        m_sin_theta = std::sqrt(1 - m_cos_theta * m_cos_theta);

        if (m_sin_theta == 0) throw std::runtime_error("CircleFit3D: Orientation vector cannot be parallel to z.");

        m_cos_phi = m_orientation.x / m_sin_theta;
        m_sin_phi = m_orientation.y / m_sin_theta;

        m_cos_alpha = std::cos(m_alpha);
        m_sin_alpha = std::sin(m_alpha);

        m_originc = this->_GetLinePoint(m_length,true);

        m_normal = Vector {-m_sin_alpha*m_cos_theta*m_cos_phi - m_cos_alpha*m_sin_phi,
                           -m_sin_alpha*m_cos_theta*m_sin_phi + m_cos_alpha*m_cos_phi,
                            m_sin_alpha*m_sin_theta };
        
        m_center = m_originc + m_radius * Vector { m_cos_alpha*m_cos_theta*m_cos_phi - m_sin_alpha*m_sin_phi,
                                                   m_cos_alpha*m_cos_theta*m_sin_phi + m_sin_alpha*m_cos_phi,
                                                  -m_cos_alpha*m_sin_theta };
        
        m_origin2 = this->_GetCirclePoint(m_phi_max);

        double cos_varphi = std::cos(m_phi_max);
        double sin_varphi = std::sin(m_phi_max);

        m_orientation2 = Vector{ sin_varphi * (m_cos_alpha*m_cos_theta*m_cos_phi - m_sin_alpha*m_sin_phi) + cos_varphi*m_sin_theta*m_cos_phi,
                                 sin_varphi * (m_cos_alpha*m_cos_theta*m_sin_phi + m_sin_alpha*m_cos_phi) + cos_varphi*m_sin_theta*m_sin_phi,
                                -sin_varphi *  m_cos_alpha*m_sin_theta + cos_varphi*m_cos_theta };
    }

    double CircleFit3D::_LineSqDist(RecoPoint point, bool first_line) const
    {
        Vector orig,orient;
        if (first_line) { orig = m_origin;  orient = m_orientation;  }
        else            { orig = m_origin2; orient = m_orientation2; }

        double t_close = orient * (point.AsVector() - orig) / orient.SqMagnitude();
        if (first_line  && (t_close > m_length)) t_close = m_length;
        if (!first_line && (t_close < 0))        t_close = 0;

        Vector line_vector = orig + t_close * orient - point.AsVector();

        return line_vector.SqMagnitude();
    }

    double CircleFit3D::_CircleSqDist(RecoPoint point) const
    {
        Vector projection = (point.AsVector()-m_center)-((point.AsVector()-m_center)*m_normal)*m_normal;
        projection.Normalize();

        Vector circle_point   = m_radius * projection + m_center;
        Vector circle_vector  = circle_point - point.AsVector();
        Vector arc_beg_vector = m_originc - point.AsVector();
        Vector arc_end_vector = m_origin2 - point.AsVector();

        // Checking if the point is on the arc.
        // double angle1 = (m_originc - m_center).Angle(projection);
        // double angle2 = (m_origin2 - m_center).Angle(projection);
        // if (angle1 + angle2 > m_phi_max)
        // {
        //     circle_vector = angle1 < angle2 ? arc_beg_vector : arc_end_vector;
        // }

        return find_min(circle_vector.SqMagnitude(),arc_beg_vector.SqMagnitude(),arc_end_vector.SqMagnitude());
    }
    
    double CircleFit3D::_SqDistance(RecoPoint point) const
    {
        double sqdist1 = _LineSqDist(point,true);
        double sqdist2 = _CircleSqDist(point);
        double sqdist3 = _LineSqDist(point,false);

        return find_min(sqdist1,sqdist2,sqdist3);
    }

    double CircleFit3D::_SumSq() const
    {
        double sum = 0;
        for (RecoPoint point : m_fit_data)
            sum += point.count*_SqDistance(point);
        // sum /= m_total_count;

        // Penalize m_origin2 outside of the OFTPC.
        // if (m_origin2.x > constants::xmax)
        // {
        //     sum += std::pow((m_origin2.x - constants::xmax)/2.5e+9,2) / m_total_count;
        // }

        return sum;
    }

    void CircleFit3D::_EvalSumSq(int& npar, double* gin, double& sumsq, double* par, int iflag)
    {
        m_length  = par[0];
        m_alpha   = par[1];
        m_radius  = par[2];
        m_phi_max = par[3];

        if (npar == 6)
        {
            m_theta   = par[4];
            m_varphi  = par[5];
        }

        _UpdateCurve();

        sumsq = _SumSq();
    }

    void CircleFit3D::_Eval(int& npar, double* gin, double& sumsq, double* par, int iflag)
    {
        return lastfit->_EvalSumSq(npar,gin,sumsq,par,iflag);
    }

    Vector CircleFit3D::_GetAvgField(const Field<Vector>& magfield, double step) const
    {
        int i = 0;
        double param = 0;
        Vector bfield = {0,0,0};

        while (param/m_radius < m_phi_max) 
        {
            Vector cur_point = _GetCirclePoint(param/m_radius);
            if(IsInTPC(cur_point))
            {
                bfield += magfield.GetField(cur_point);
                i++;
            }
            param += step;
        }

        return bfield/i;
    }
    
    Vector CircleFit3D::_GetMiddleField(const Field<Vector>& magfield, double tolerance) const
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