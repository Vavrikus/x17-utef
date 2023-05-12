// X17 dependencies
#include "CircleFit3D.h"

namespace X17
{
    CircleFit3D* CircleFit3D::lastfit = nullptr;

    Vector CircleFit3D::GetLinePoint(double param, bool first_line)
    {
        if (first_line) return origin  + param * orientation;
        else            return origin2 + param * orientation2;
    }

    Vector CircleFit3D::GetCirclePoint(double varphi)
    {
        double cos_varphi = cos(varphi);
        double sin_varphi = sin(varphi);

        return originc + radius*Vector{ (1-cos_varphi)*(cos_alpha*cos_theta*cos_phi - sin_alpha*sin_phi) + sin_varphi*sin_theta*cos_phi,
                                        (1-cos_varphi)*(cos_alpha*cos_theta*sin_phi + sin_alpha*cos_phi) + sin_varphi*sin_theta*sin_phi,
                                        -(1-cos_varphi)*cos_alpha*sin_theta + sin_varphi*cos_theta};
    }

    void CircleFit3D::UpdateCurve()
    {
        cos_theta = orientation.z;
        sin_theta = sqrt(1-cos_theta*cos_theta);

        cos_phi = orientation.x/sin_theta;
        sin_phi = orientation.y/sin_theta;

        cos_alpha = cos(alpha);
        sin_alpha = sin(alpha);

        originc = this->GetLinePoint(length,true);

        normal = Vector {-sin_alpha*cos_theta*cos_phi - cos_alpha*sin_phi,
                            -sin_alpha*cos_theta*sin_phi + cos_alpha*cos_phi,
                            sin_alpha*sin_theta };
        
        center = originc + radius * Vector { cos_alpha*cos_theta*cos_phi - sin_alpha*sin_phi,
                                                cos_alpha*cos_theta*sin_phi + sin_alpha*cos_phi,
                                            -cos_alpha*sin_theta };
        
        origin2 = this->GetCirclePoint(phi_max);

        double cos_varphi = cos(phi_max);
        double sin_varphi = sin(phi_max);

        orientation2 = Vector{ sin_varphi*(cos_alpha*cos_theta*cos_phi - sin_alpha*sin_phi) + cos_varphi*sin_theta*cos_phi,
                                sin_varphi*(cos_alpha*cos_theta*sin_phi + sin_alpha*cos_phi) + cos_varphi*sin_theta*sin_phi,
                                -sin_varphi*cos_alpha*sin_theta + cos_varphi*cos_theta };
    }

    double CircleFit3D::LineSqDist(const RecoPoint& point, bool first_line)
    {
        Vector *orig,*orient;
        if (first_line) { orig = &origin;  orient = &orientation;  }
        else            { orig = &origin2; orient = &orientation2; }

        double t_close = *orient*(point.AsVector()-*orig)/orient->SqMagnitude();
        if (first_line  && (t_close > length)) t_close = length;
        if (!first_line && (t_close < 0))      t_close = 0;

        Vector line_vector = *orig+t_close*(*orient)-point.AsVector();

        return line_vector.SqMagnitude();
    }

    double CircleFit3D::CircleSqDist(const RecoPoint& point)
    {
        Vector projection = (point.AsVector()-center)-((point.AsVector()-center)*normal)*normal;
        projection.Normalize();

        Vector circle_point   = radius * projection + center;
        Vector circle_vector  = circle_point - point.AsVector();
        Vector arc_beg_vector = originc - point.AsVector();
        Vector arc_end_vector = origin2 - point.AsVector();

        return find_min(circle_vector.SqMagnitude(),arc_beg_vector.SqMagnitude(),arc_end_vector.SqMagnitude());
    }
    
    double CircleFit3D::SqDistance(const RecoPoint& point)
    {
        double sqdist1 = LineSqDist(point,true);
        double sqdist2 = CircleSqDist(point);
        double sqdist3 = LineSqDist(point,false);

        return find_min(sqdist1,sqdist2,sqdist3);
    }

    double CircleFit3D::SumSq()
    {
        double sum = 0;
        for (RecoPoint point : fit_data) sum += point.count*SqDistance(point);
        return sum;
    }

    void CircleFit3D::EvalSumSq(int& npar, double* gin, double& sumsq, double* par, int iflag)
    {
        length  = par[0];
        alpha   = par[1];
        radius  = par[2];
        phi_max = par[3];

        UpdateCurve();

        sumsq = SumSq();
    }

    void CircleFit3D::Eval(int& npar, double* gin, double& sumsq, double* par, int iflag)
    {
        return lastfit->EvalSumSq(npar,gin,sumsq,par,iflag);
    }

    Vector CircleFit3D::GetAvgField(const Field<Vector>& magfield, double step)
    {
        int i = 0;
        double param = 0;
        Vector bfield = {0,0,0};

        while (param/radius < phi_max) 
        {
            Vector cur_point = GetCirclePoint(param/radius);
            if(IsInSector(cur_point))
            {
                bfield += magfield.GetField(cur_point);
                i++;
            }
            param += step;
        }

        return bfield/i;
    }
    
    Vector CircleFit3D::GetMiddleField(const Field<Vector>& magfield, double tolerance)
    {
        using namespace constants;

        double xmiddle = (xmax + xmin) / 2;
        double low  = 0;
        double high = phi_max;
        double mid  = (low + high) / 2;
        Vector vmid = GetCirclePoint(mid);

        while ((abs(xmiddle-vmid.x) > tolerance) && (low*(1+1e-15) < high))
        {
            mid = (low + high) / 2;
            vmid  = GetCirclePoint(mid);

            if (vmid.x < xmiddle) low  = mid;
            else                   high = mid;
        }
        
        return magfield.GetField(vmid);
    }

    CircleFit3D::CircleFit3D(const Vector& orig, const Vector& orient)
    {
        this->origin      = orig;
        this->orientation = orient;
        this->orientation.Normalize();

        this->gFitter = lastfit->gFitter; // Setting the fitter again every time would take much more time.
        lastfit = this;
    }

    void CircleFit3D::SetFitter(int parameters, bool print)
    {        
        gFitter = TVirtualFitter::Fitter(nullptr,parameters); // the second number is number of parameters
        gFitter->SetFCN(this->Eval);

        if(!print)
        {
            double arg = -1;
            gFitter->ExecuteCommand("SET PRINTOUT",&arg,1);
            gFitter->ExecuteCommand("SET NOW", &arg ,1);
        }
    }

    void CircleFit3D::FitCircle3D(double max_iter, double toleration)
    {
        gFitter->SetParameter(0,"length",length,0.01,-10,5);
        gFitter->SetParameter(1,"alpha",alpha,0.001,-M_PI/2,(3/2)*M_PI);
        gFitter->SetParameter(2,"radius",radius,0.01,10,50);
        gFitter->SetParameter(3,"phi_max",phi_max,0.001,0.15,M_PI/1.5);

        double arglist[2] = {max_iter,toleration};   // max iterations, step size (toleration)
        gFitter->ExecuteCommand("MIGRAD",arglist,2); // last one num of prints (verbosity)

        length  = gFitter->GetParameter(0);
        l_err   = gFitter->GetParError(0);
        alpha   = gFitter->GetParameter(1);
        a_err   = gFitter->GetParError(1);
        radius  = gFitter->GetParameter(2);
        r_err   = gFitter->GetParError(2);
        phi_max = gFitter->GetParameter(3);
        phi_err = gFitter->GetParError(3);

        UpdateCurve();
    }
    
    void CircleFit3D::PrintFitParams()
    {
        std::cout << "\nCIRCLE FIT PARAMETERS:\n";
        std::cout << "Length:  " << length  << " +- " << l_err   << "\n";
        std::cout << "Alpha:   " << alpha   << " +- " << a_err   << "\n";
        std::cout << "Radius:  " << radius  << " +- " << r_err   << "\n";
        std::cout << "Phi_max: " << phi_max << " +- " << phi_err << "\n\n";
    }

    TGraph2D* CircleFit3D::GetGraph(double step, double dist)
    {
        TGraph2D* fit_graph = new TGraph2D();
        
        // 1st line
        double param = 0;
        while (param < length)
        {
            Vector cur_pos = GetLinePoint(param,true);
            if(IsInSector(cur_pos,dist)) 
                fit_graph->AddPoint(cur_pos.x,cur_pos.y,cur_pos.z);
            param += step;
        }

        // circle
        param = 0;
        while (param/radius < phi_max)
        {
            Vector cur_pos = GetCirclePoint(param/radius);
            if(IsInSector(cur_pos,dist)) 
                fit_graph->AddPoint(cur_pos.x,cur_pos.y,cur_pos.z);
            param += step;
        }

        // 2nd line
        param = 0;
        while (param < 20)
        {
            Vector cur_pos = GetLinePoint(param,false);
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
        if(middle) bfield = GetMiddleField(magfield);
        else       bfield = GetAvgField(magfield);

        double b_proj = normal*bfield;

        double betasq = 1 / (1 + pow((E0 / (c * (cm2m * radius) * b_proj)), 2));
        double Ekin   = E0 * (1 / sqrt(1 - betasq) - 1);

        return Ekin;
    }

    TGraph* CircleFit3D::GetEnergyGraph(const Field<Vector>& magfield, double step)
    {
        using namespace constants;

        TGraph* graph = new TGraph();
        double param = 0;

        while (param/radius < phi_max) 
        {
            Vector cur_point = GetCirclePoint(param/radius);
            if(IsInSector(cur_point))
            {
                Vector bfield = magfield.GetField(cur_point);
                double b_proj = normal * bfield;

                double betasq = 1 / (1 + pow((E0 / (c * (cm2m * radius) * b_proj)), 2));
                double Ekin   = E0 * (1 / sqrt(1 - betasq) - 1);
                graph->AddPoint(param,Ekin);
            }
            param += step;
        }
        
        graph->SetTitle("Energy along the track;Parameter [cm];Energy [eV]");
        return graph;
    }
} // namespace X17