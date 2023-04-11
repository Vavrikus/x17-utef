#pragma once

#include <functional>
#include <iostream>
#include <vector>

#include "TGraph2D.h"

#include "CircleFit3D.h"
#include "VectorField.h"
#include "X17Utilities.h"

template<int N>
class RK4
{
    typedef Matrix<N,1> VectorN;
    typedef std::function<void(const double&,const VectorN&,VectorN&)> VecFn;
    typedef std::function<bool(const double&,const VectorN&)> EndFn;

private:
    VecFn dif_eq;
    VectorN current;
    double param;
    double step;
    EndFn end_fn;

    std::vector<VectorN> results;
    
public:
    RK4(double start, double step, VecFn equation, VectorN initial, EndFn end_condition)
        : param(start),step(step),dif_eq(equation),end_fn(end_condition),current(initial)
    {
        results.push_back(current);
    }

    void Run()
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
            dif_eq(param+step  , current+k3,   k4);
            k4 *= step;

            VectorN diff = (1.0/6.0)*(k1+2*k2+2*k3+k4);
            current += diff;
            param += step;
            results.push_back(current);
        }
    }

    const std::vector<VectorN>& GetResults() const {return results;}
    int GetSize() const {return results.size();}
    double GetStep() const {return step;}
};

/// @brief Function for setting initial parameters of RK4 simulation (currently taking origin literally)
/// @param kin_en      Kinetic energy of the particle [eV]
/// @param origin      Coordinates of the origin [cm]
/// @param orientation Initial direction of motion
/// @return Vector of position and four-velocity
Matrix<8,1> GetInitParams(const double& kin_en, const Vector& origin, const Vector& orientation)
{
    Vector n_orient = orientation;
    n_orient.Normalize();

    Vector start = origin/X17::m2cm;//(origin + (X17::xmin/n_orient.vx)*n_orient)/X17::m2cm;

    double gamma    = 1 + kin_en/X17::E0;
    double velocity = sqrt((1-pow(1/gamma,2))*pow(X17::c,2));

    return Matrix<8,1>({
        0,
        start.vx,
        start.vy,
        start.vz,
        gamma*X17::c,
        gamma*velocity*n_orient.vx,
        gamma*velocity*n_orient.vy,
        gamma*velocity*n_orient.vz
    });
}

Matrix<4,4> GetEMtensor(VectorField* magfield, const Vector& position)
{
    using namespace X17;

    Vector b = magfield->GetField(position);

    return Matrix<4,4>({0,                -m2cm*efield.vx/c, -m2cm*efield.vy/c, -m2cm*efield.vz/c,
                        m2cm*efield.vx/c,  0,                 b.vz,             -b.vy,
                        m2cm*efield.vy/c, -b.vz,              0,                 b.vx,
                        m2cm*efield.vz/c,  b.vy,             -b.vx,              0               });
}

void EMMotion(VectorField* magfield, const bool& electron, const double& tau, const Matrix<8,1>& params, Matrix<8,1>& output)
{
    const double charge = electron ? -X17::e : X17::e;
    Matrix<4,1> fourvelocity({params.at(4,0),params.at(5,0),params.at(6,0),params.at(7,0)});
    Vector position = {params.at(1,0),params.at(2,0),params.at(3,0)};

    Matrix<4,4> EMtensor = GetEMtensor(magfield,position);
    Matrix<4,1> out1 = (charge/X17::m0)*EMtensor*fourvelocity;
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
    using namespace X17;
    return (!X17::IsInSector(m2cm*params.at(1,0),m2cm*params.at(2,0),m2cm*params.at(3,0),-0.5)) && (m2cm*params.at(1,0) > X17::xmin);
}

RK4<8>* GetTrackRK(VectorField* magfield, const bool& electron, const double& step, const double& kin_en, const Vector& origin, const Vector& orientation)
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
        using namespace X17;
        Vector point = {m2cm*r.at(1,0),m2cm*r.at(2,0),m2cm*r.at(3,0)};
        if (X17::IsInSector(point)) output->AddPoint(point.vx,point.vy,point.vz);
    }
    return output;
}

class RKFit
{
    static RKFit* lastfit;

    VectorField* magfield;
    bool electron;
    double step;
    Vector origin, orientation;

    std::vector<DataPoint> fit_data;
    RK4<8>* curr_rk = nullptr;

    TVirtualFitter* fitter;
    double kin_en,e_err;

    Vector GetPoint(const int& index)
    {
        using namespace X17;
        auto rk_pts = curr_rk->GetResults();
        return {m2cm*rk_pts[index].at(1,0),m2cm*rk_pts[index].at(2,0),m2cm*rk_pts[index].at(3,0)};
    }

    double SqDist(const int& index, const Vector& point)
    {
        Vector rk_point = GetPoint(index);
        return rk_point.SqDist(point);
    }

    // isn't scaled by step (not currently necessary)
    double DistDerivative(const int& index, const Vector& point)
    {
        return SqDist(index+1,point)-SqDist(index,point);
    }

    double GetSqDist(const Vector& point)
    {
        int min_index = 0;
        int max_index = curr_rk->GetSize()-2;

        double min_der = DistDerivative(min_index,point);
        double max_der = DistDerivative(max_index,point);

        // convex function minimum search
        while (min_index+1 != max_index)
        {
            int mid_index = (min_index+max_index)/2;
            double mid_der = DistDerivative(mid_index,point);

            if (mid_der*max_der > 0) {max_der = mid_der; max_index = mid_index;}
            else                     {min_der = mid_der; min_index = mid_index;}
        }

        // minimal distance from two lines
        Vector orig    = GetPoint(max_index);
        Vector orient1 = GetPoint(max_index+1) - orig;
        Vector orient2 = GetPoint(max_index-1) - orig;

        double dist1 = LineSqDist(orig,orient1,orient1.Magnitude(),point);
        double dist2 = LineSqDist(orig,orient2,orient2.Magnitude(),point);

        return min(dist1,dist2);
    }

    void SumSqDist(int& npar, double* gin, double& sumsq, double* par, int iflag)
    {
        kin_en = par[0];
        
        if (curr_rk != nullptr) delete curr_rk;
        curr_rk = GetTrackRK(magfield,electron,step,kin_en,origin,orientation);
        curr_rk->Run();
        
        sumsq = 0;
        for (auto p : fit_data) sumsq += p.count*GetSqDist(p.ToVector());
        // std::cout << "Kin_en: " << kin_en << " sumsq: " << sumsq << "\n";
    }

    static void Eval(int& npar, double* gin, double& sumsq, double* par, int iflag)
    {
        lastfit->SumSqDist(npar,gin,sumsq,par,iflag);
    }

public:
    RKFit(VectorField* magfield, const bool& electron, const double& step, const Vector& origin, const Vector& orientation, const std::vector<DataPoint>& fit_data)
        : magfield(magfield), electron(electron), step(step), origin(origin), orientation(orientation), fit_data(fit_data)
    {lastfit = this;}

    void SetEnergy(const double& kin_en) {this->kin_en = kin_en;}

    void SetFitter(int parameters = 1, bool print = true)
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

    void FitRK(double max_iter = 500, double toleration = 0.001)
    {
        fitter->SetParameter(0,"kin_en",kin_en,1000,1000000,16000000);

        double arglist[2] = {max_iter,toleration};  // max iterations, step size (toleration)
        fitter->ExecuteCommand("MIGRAD",arglist,2); // last one num of prints (verbosity)

        kin_en = fitter->GetParameter(0);
        e_err  = fitter->GetParError(0);
    }

    void PrintFitParams()
    {
        cout << "\nRK FIT RESULT:\n";
        cout << "Kinetic energy:  " << kin_en << " +- " << e_err << "\n\n";
    }

    TGraph2D* GetFitGraph() {return GetGraphRK(curr_rk);}
};

RKFit* RKFit::lastfit = nullptr;