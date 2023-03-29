#include <functional>
#include <vector>

#include "TGraph2D.h"

#include "../VectorField.h"
#include "../X17Utilities.h"

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

    std::vector<VectorN> GetResults() const {return results;}

    int GetSize() const {return results.size();}
};

Matrix<8,1> GetInitParams(const double& kin_en, const Vector& origin, const Vector& orientation)
{
    Vector n_orient = orientation;
    n_orient.Normalize();

    Vector start = origin/100.0;//(origin + (X17::xmin/n_orient.vx)*n_orient)/100.0;

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

    constexpr double m2cm = 100;
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
    return (!X17::IsInSector(100*params.at(1,0),100*params.at(2,0),100*params.at(3,0),-0.5)) && (100*params.at(1,0) > X17::xmin);
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

    constexpr double m2cm = 100;
    for (auto r : results) output->AddPoint(m2cm*r.at(1,0),m2cm*r.at(2,0),m2cm*r.at(3,0));
    return output;
}