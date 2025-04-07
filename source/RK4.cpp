// C++ dependencies
#include <cmath>
#include <iostream>

// ROOT dependencies
#include "TGraph2D.h"

// X17 dependencies
#include "Field.h"
#include "Matrix.h"
#include "Points.h"
#include "RK4.h"
#include "Vector.h"
#include "X17Utilities.h"

namespace X17
{
    //// Functions related to the RK4 templated class.

    Matrix<8,1> GetInitParams(double kin_en, Vector origin, Vector orientation)
    {
        using namespace constants;

        Vector n_orient = orientation;
        n_orient.Normalize();

        Vector start = origin / m2cm;

        double gamma    = 1 + kin_en / E0;
        double velocity = c * sqrt(1 - pow(1.0 / gamma,2));

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

    Matrix<4,4> GetEMtensor(const Field<Vector>& magfield, Vector position)
    {
        using namespace constants;

        Vector b = magfield.GetField(m2cm * position);

        // The electric field in the TPC has units of [V/cm] so we transform to [V/m].
        return Matrix<4,4>({
            0,                -m2cm*efield.x/c, -m2cm*efield.y/c, -m2cm*efield.z/c,
            m2cm*efield.x/c,  0,                 b.z,             -b.y,
            m2cm*efield.y/c, -b.z,               0,                b.x,
            m2cm*efield.z/c,  b.y,              -b.x,              0
        });
    }

    void EMMotion(const Field<Vector>& magfield, bool electron, double tau, const Matrix<8,1>& params, Matrix<8,1>& output)
    {
        using namespace constants;

        const double charge = electron ? -e : e;
        Matrix<4,1> fourvelocity({params.at(4,0),params.at(5,0),params.at(6,0),params.at(7,0)});
        Vector position = {params.at(1,0),params.at(2,0),params.at(3,0)};

        Matrix<4,4> EMtensor = GetEMtensor(magfield,position);
        Matrix<4,1> out1 = (charge/m0)*EMtensor*fourvelocity;
        output = Matrix<8,1>({
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

    bool IsOutOfSector(double tau, const Matrix<8,1>& params)
    {
        using namespace constants;
        return (!IsInSector(m2cm * params.at(1,0), m2cm * params.at(2,0), m2cm * params.at(3,0),-0.5)); // && (m2cm * params.at(1,0) > xmin);
    }

    RK4<8>* GetTrackRK(const Field<Vector>& magfield, bool electron, double step, double kin_en, Vector origin, Vector orientation)
    {
        Matrix<8,1> init = GetInitParams(kin_en,origin,orientation);
        auto f = [&magfield, electron](double t, const Matrix<8, 1>& y, Matrix<8, 1>& dydt) {
            EMMotion(magfield, electron, t, y, dydt);
        };
        return new RK4<8>(0,step,f,init,IsOutOfSector);
    }

    TGraph2D* GetGraphRK(const RK4<8>* rk_track)
    {
        TGraph2D* output = new TGraph2D();
        auto results = rk_track->GetResults();

        for (auto r : results)
        {
            using namespace constants;
            Vector point = {m2cm * r.at(1,0), m2cm * r.at(2,0), m2cm * r.at(3,0)};
            if (IsInSector(point)) output->AddPoint(point.x,point.y,point.z);
        }
        return output;
    }

    //// Class RKFit static variables.

    RKFit* RKFit::lastfit = nullptr;

    
    
    
    
    //// Public methods of the RKFit class.

    RKFit::RKFit(Field<Vector>* magfield, bool electron, double step, Vector origin, Vector orientation, const std::vector<RecoPoint>& fit_data)
        : m_magfield(magfield), m_electron(electron), m_step(step), m_origin(origin), m_orientation(orientation), m_fit_data(fit_data)
    { lastfit = this; }

    void RKFit::SetFitter(int parameters, bool print)
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

    void RKFit::FitRK(double max_iter, double toleration)
    {
        m_fitter->SetParameter(0,"kin_en",m_kin_en,1000,1000000,16000000);

        double arglist[2] = {max_iter,toleration};    // Maximum iterations, step size (toleration).
        m_fitter->ExecuteCommand("MIGRAD",arglist,2); // The last parameter is the num of prints (verbosity).

        m_kin_en = m_fitter->GetParameter(0);
        m_e_err  = m_fitter->GetParError(0);
    }

    void RKFit::PrintFitParams() const
    {
        std::cout << "\nRK FIT RESULT:\n";
        std::cout << "Kinetic energy:  " << m_kin_en << " +- " << m_e_err << "\n\n";
    }





    //// Private methods of the RKFit class.

    RKPoint RKFit::_GetPoint(int index) const
    {
        using namespace constants;
        auto rk_pts = m_curr_rk->GetResults();
        return RKPoint(m2cm * rk_pts[index].at(1,0), m2cm * rk_pts[index].at(2,0), m2cm * rk_pts[index].at(3,0), 1e+9 / c * rk_pts[index].at(0,0));
    }

    double RKFit::_SqDist(int index, Vector point) const
    {
        RKPoint rk_point = _GetPoint(index);
        return rk_point.AsVector().SqDist(point);
    }

    double RKFit::_GetSqDist(Vector point) const
    {
        int min_index = 0;
        int max_index = m_curr_rk->GetSize() - 2;

        double min_der = _DistDerivative(min_index,point);
        double max_der = _DistDerivative(max_index,point);

        // Convex function minimum search.
        while (min_index + 1 != max_index)
        {
            int mid_index = (min_index + max_index) / 2;
            double mid_der = _DistDerivative(mid_index,point);

            if (mid_der*max_der > 0) {max_der = mid_der; max_index = mid_index;}
            else                     {min_der = mid_der; min_index = mid_index;}
        }

        // Minimal distance from two lines.
        Vector orig    = _GetPoint(max_index).AsVector();
        Vector orient1 = _GetPoint(max_index + 1).AsVector() - orig;
        Vector orient2 = _GetPoint(max_index - 1).AsVector() - orig;

        double dist1 = LineSqDist(orig,orient1,orient1.Magnitude(),point);
        double dist2 = LineSqDist(orig,orient2,orient2.Magnitude(),point);

        return std::min(dist1,dist2);
    }

    void RKFit::_SumSqDist(int& npar, double* gin, double& sumsq, double* par, int iflag)
    {
        m_kin_en = par[0];
        
        if (m_curr_rk != nullptr) delete m_curr_rk;
        m_curr_rk = GetTrackRK(*m_magfield,m_electron,m_step,m_kin_en,m_origin,m_orientation);
        m_curr_rk->Integrate();
        
        sumsq = 0;
        for (auto p : m_fit_data) sumsq += p.count * _GetSqDist(p.AsVector());
    }
} // namespace X17