#pragma once

// X17 dependencies
#include "RK4.h"

namespace X17
{
    template<int N>        
    RK4<N>::RK4(double start, double step, VecFn equation, VectorN initial, EndFn end_condition)
        : param(start),step(step),dif_eq(equation),end_fn(end_condition),current(initial)
    {
        results.push_back(current);
    }

    template<int N>
    void RK4<N>::Run()
    {
        while(!end_fn(step,current))
        {
            VectorN k1,k2,k3,k4;
            dif_eq(param,        current,      k1);
            k1 *= step;
            dif_eq(param+step/2, current+k1/2, k2);
            k2 *= step;
            dif_eq(param+step/2, current+k2/2, k3);
            k3 *= step;
            dif_eq(param+step,   current+k3,   k4);
            k4 *= step;

            VectorN diff = (1.0/6.0)*(k1+2*k2+2*k3+k4);
            current += diff;
            param   += step;
            results.push_back(current);
        }
    }

    Matrix<8,1> GetInitParams(const double& kin_en, const Vector& origin, const Vector& orientation)
    {
        using namespace constants;

        Vector n_orient = orientation;
        n_orient.Normalize();

        Vector start = origin/m2cm;

        double gamma    = 1 + kin_en/E0;
        double velocity = sqrt((1-pow(1/gamma,2))*pow(c,2));

        return Matrix<8,1>({
            0,
            start.x,
            start.y,
            start.z,
            gamma*c,
            gamma*velocity*n_orient.x,
            gamma*velocity*n_orient.y,
            gamma*velocity*n_orient.z
        });
    }

    Matrix<4,4> GetEMtensor(const Field<Vector>& magfield, const Vector& position)
    {
        using namespace constants;

        Vector b = magfield.GetField(position);

        return Matrix<4,4>({0,                -m2cm*efield.x/c, -m2cm*efield.y/c, -m2cm*efield.z/c,
                            m2cm*efield.x/c,  0,                 b.z,             -b.y,
                            m2cm*efield.y/c, -b.z,              0,                 b.x,
                            m2cm*efield.z/c,  b.y,             -b.x,              0               });
    }

    void EMMotion(const Field<Vector>& magfield, const bool& electron, const double& tau, const Matrix<8,1>& params, Matrix<8,1>& output)
    {
        using namespace constants;

        const double charge = electron ? -e : e;
        Matrix<4,1> fourvelocity({params.at(4,0),params.at(5,0),params.at(6,0),params.at(7,0)});
        Vector position = {params.at(1,0),params.at(2,0),params.at(3,0)};

        Matrix<4,4> EMtensor = GetEMtensor(magfield,position);
        Matrix<4,1> out1 = (charge/m0)*EMtensor*fourvelocity;
        output =  Matrix<8,1>({
            params.at(4,0),
            params.at(5,0),
            params.at(6,0),
            params.at(7,0),
            out1.at(0,0),
            out1.at(1,0),
            out1.at(2,0),
            out1.at(3,0)
        });
    }

    bool IsOutOfSector(const double& tau, const Matrix<8,1>& params)
    {
        using namespace constants;
        return (!IsInSector(m2cm*params.at(1,0),m2cm*params.at(2,0),m2cm*params.at(3,0),-0.5)) && (m2cm*params.at(1,0) > xmin);
    }

    RK4<8>* GetTrackRK(const Field<Vector>& magfield, const bool& electron, const double& step, const double& kin_en, const Vector& origin, const Vector& orientation)
    {
        Matrix<8,1> init = GetInitParams(kin_en,origin,orientation);
        using VecFn = std::function<void(const double&, const Matrix<8, 1>&, Matrix<8, 1>&)>;
        VecFn f = std::bind(EMMotion,magfield,electron,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);
        return new RK4<8>(0,step,f,init,IsOutOfSector);
    }

    TGraph2D* GetGraphRK(const RK4<8>* rk_track)
    {
        TGraph2D* output = new TGraph2D();
        auto results = rk_track->GetResults();

        for (auto r : results)
        {
            using namespace constants;
            Vector point = {m2cm*r.at(1,0),m2cm*r.at(2,0),m2cm*r.at(3,0)};
            if (IsInSector(point)) output->AddPoint(point.x,point.y,point.z);
        }
        return output;
    }

    RKFit* RKFit::lastfit = nullptr;

    const RKPoint& RKFit::GetPoint(const int& index) const
    {
        using namespace constants;
        auto rk_pts = curr_rk->GetResults();
        return {1e+9*rk_pts[index].at(0,0),m2cm*rk_pts[index].at(1,0),m2cm*rk_pts[index].at(2,0),m2cm*rk_pts[index].at(3,0)};
    }

    double RKFit::SqDist(const int& index, const Vector& point)
    {
        RKPoint rk_point = GetPoint(index);
        return rk_point.AsVector().SqDist(point);
    }

    double RKFit::GetSqDist(const Vector& point)
    {
        int min_index = 0;
        int max_index = curr_rk->GetSize()-2;

        double min_der = DistDerivative(min_index,point);
        double max_der = DistDerivative(max_index,point);

        // convex function minimum search
        while (min_index + 1 != max_index)
        {
            int mid_index = (min_index+max_index)/2;
            double mid_der = DistDerivative(mid_index,point);

            if (mid_der*max_der > 0) {max_der = mid_der; max_index = mid_index;}
            else                     {min_der = mid_der; min_index = mid_index;}
        }

        // minimal distance from two lines
        Vector orig    = GetPoint(max_index).AsVector();
        Vector orient1 = GetPoint(max_index + 1).AsVector() - orig;
        Vector orient2 = GetPoint(max_index - 1).AsVector() - orig;

        double dist1 = LineSqDist(orig,orient1,orient1.Magnitude(),point);
        double dist2 = LineSqDist(orig,orient2,orient2.Magnitude(),point);

        return std::min(dist1,dist2);
    }

    void RKFit::SumSqDist(int& npar, double* gin, double& sumsq, double* par, int iflag)
    {
        kin_en = par[0];
        
        if (curr_rk != nullptr) delete curr_rk;
        curr_rk = GetTrackRK(*magfield,electron,step,kin_en,origin,orientation);
        curr_rk->Run();
        
        sumsq = 0;
        for (auto p : fit_data) sumsq += p.count*GetSqDist(p.AsVector());
    }

    RKFit::RKFit(Field<Vector>* magfield, const bool& electron, const double& step, const Vector& origin,
            const Vector& orientation, const std::vector<RecoPoint>& fit_data)
        : magfield(magfield), electron(electron), step(step), origin(origin), orientation(orientation), fit_data(fit_data)
    { lastfit = this; }

    void RKFit::SetFitter(int parameters = 1, bool print = true)
    {        
        fitter = TVirtualFitter::Fitter(nullptr,parameters); // the second number is number of parameters
        fitter->SetFCN(this->Eval);

        if(!print)
        {
            double arg = -1;
            fitter->ExecuteCommand("SET PRINTOUT",&arg,1);
            fitter->ExecuteCommand("SET NOW", &arg ,1);
        }
    }

    void RKFit::FitRK(double max_iter = 500, double toleration = 0.001)
    {
        fitter->SetParameter(0,"kin_en",kin_en,1000,1000000,16000000);

        double arglist[2] = {max_iter,toleration};  // max iterations, step size (toleration)
        fitter->ExecuteCommand("MIGRAD",arglist,2); // last one num of prints (verbosity)

        kin_en = fitter->GetParameter(0);
        e_err  = fitter->GetParError(0);
    }

    void RKFit::PrintFitParams()
    {
        std::cout << "\nRK FIT RESULT:\n";
        std::cout << "Kinetic energy:  " << kin_en << " +- " << e_err << "\n\n";
    }
} // namespace X17