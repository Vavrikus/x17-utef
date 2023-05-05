#pragma once

#include <cmath>

#include "TGraph2D.h"
#include "TVirtualFitter.h"

#include "Matrix.h"
#include "Points.h"
#include "Vector.h"
#include "X17Utilities.h"

namespace X17
{
    /// @brief A templated class implementing the Runge-Kutta method of order 4 for solving differential equations.
    /// @tparam N The dimension of the system of differential equations.
    template<int N>
    class RK4
    {
        typedef Matrix<N,1> VectorN;
        typedef std::function<void(const double&,const VectorN&,VectorN&)> VecFn;
        typedef std::function<bool(const double&,const VectorN&)> EndFn;

    private:
        VecFn dif_eq;    // The function object that returns the derivative of the dependent variable with respect to the independent variable.
        VectorN current; // The current value of the dependent variable.
        double param;    // The current value of the independent variable.
        double step;     // The step size for the integration.
        EndFn end_fn;    // The function object that determines the end condition for the integration.

        std::vector<VectorN> results; // A vector that stores the values of the dependent variable at different points in the integration.
        
    public:
        /// @brief Constructs an RK4 object with the specified initial conditions and integration parameters.
        /// @param start The initial value of the independent variable.
        /// @param step The step size for the integration.
        /// @param equation The function object that returns the derivative of the dependent variable with respect to the independent variable.
        /// @param initial The initial value of the dependent variable.
        /// @param end_condition The function object that determines the end condition for the integration.
        RK4(double start, double step, VecFn equation, VectorN initial, EndFn end_condition)
            : param(start),step(step),dif_eq(equation),end_fn(end_condition),current(initial)
        {
            results.push_back(current);
        }

        /// @brief Performs the integration using the fourth-order Runge-Kutta method.
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
                dif_eq(param+step,   current+k3,   k4);
                k4 *= step;

                VectorN diff = (1.0/6.0)*(k1+2*k2+2*k3+k4);
                current += diff;
                param   += step;
                results.push_back(current);
            }
        }

        /// @brief Returns the vector of values of the dependent variable at different points in the integration.
        /// @return The vector of values of the dependent variable at different points in the integration.
        const std::vector<VectorN>& GetResults() const { return results; }

        /// @brief Returns the number of points in the integration.
        /// @return The number of points in the integration.
        int GetSize() const { return results.size(); }

        /// @brief Returns the step size used for the integration.
        /// @return The step size used for the integration.
        double GetStep() const { return step; }
    };

    /// @brief Computes the initial parameters for a simulation using the RK4 algorithm.
    /// @param kin_en      Kinetic energy of the particle [eV].
    /// @param origin      Coordinates of the origin [cm].
    /// @param orientation Initial direction of motion.
    /// @return A matrix containing the position and four-velocity of the particle.
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

    /// @brief Calculates the electromagnetic field tensor at a given position.
    /// @param magfield Vector field representing the magnetic field.
    /// @param position The position at which to calculate the electromagnetic field.
    /// @return A 4x4 matrix representing the electromagnetic field tensor.
    Matrix<4,4> GetEMtensor(const Field<Vector>& magfield, const Vector& position)
    {
        using namespace constants;

        Vector b = magfield.GetField(position);

        return Matrix<4,4>({0,                -m2cm*efield.x/c, -m2cm*efield.y/c, -m2cm*efield.z/c,
                            m2cm*efield.x/c,  0,                 b.z,             -b.y,
                            m2cm*efield.y/c, -b.z,              0,                 b.x,
                            m2cm*efield.z/c,  b.y,             -b.x,              0               });
    }

    /// @brief Calculates the electromagnetic motion of a charged particle using the Lorentz force equation.
    /// @param magfield Reference to the magnetic field.
    /// @param electron Boolean value indicating whether the particle is an electron or a positron.
    /// @param tau The independent variable (proper time) [s].
    /// @param params Vector containing the initial position and four-velocity of the particle.
    /// @param output Vector containing the final position and four-velocity of the particle.
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

    /// @brief Checks if the particle is out of the sector after it has moved beyond the minimum x position.
    /// @param tau Proper time elapsed (independent variable) [s].
    /// @param params Vector of position and four-velocity of the particle.
    /// @return True if the particle is out of the sector after it has moved beyond the minimum x position, false otherwise.
    bool IsOutOfSector(const double& tau, const Matrix<8,1>& params)
    {
        using namespace constants;
        return (!IsInSector(m2cm*params.at(1,0),m2cm*params.at(2,0),m2cm*params.at(3,0),-0.5)) && (m2cm*params.at(1,0) > xmin);
    }

    /// @brief Returns a pointer to an instance of RK4 class for tracking a charged particle in a magnetic field
    /// @param magfield Vector field representing the magnetic field.
    /// @param electron Boolean indicating if the particle is an electron or a positron.
    /// @param step Step size for the RK4 integration [s].
    /// @param kin_en Kinetic energy of the particle [eV].
    /// @param origin Coordinates of the origin [cm].
    /// @param orientation Initial direction of motion.
    /// @return Pointer to an instance of RK4 class.
    RK4<8>* GetTrackRK(const Field<Vector>& magfield, const bool& electron, const double& step, const double& kin_en, const Vector& origin, const Vector& orientation)
    {
        Matrix<8,1> init = GetInitParams(kin_en,origin,orientation);
        using VecFn = std::function<void(const double&, const Matrix<8, 1>&, Matrix<8, 1>&)>;
        VecFn f = std::bind(EMMotion,magfield,electron,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);
        return new RK4<8>(0,step,f,init,IsOutOfSector);
    }

    /// @brief Returns a TGraph2D object representing the track of a particle obtained from a Runge-Kutta 4th order simulation.
    /// @param rk_track A pointer to a RK4 object containing the simulation results.
    /// @return A pointer to a TGraph2D object representing the particle track.
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

    /// @brief A class for fitting a 4th order Runge-Kutta generated points to a particle's trajectory.
    class RKFit
    {
        static RKFit* lastfit; // A static pointer to the last instance of RKFit class created.

        Field<Vector>* magfield; // A pointer to the magnetic field that the particle is traveling through.
        bool electron;           // A boolean indicating whether the particle is an electron (true) or a positron (false).
        double step;             // The step size used in the Runge-Kutta 4th order method.
        Vector origin;           // The initial position of the particle.
        Vector orientation;      // The initial orientation of the particle.

        std::vector<RecoPoint> fit_data; // A vector containing the reconstructed points used for fitting.
        RK4<8>* curr_rk = nullptr;       // A pointer to the RK4 object used to compute the trajectory of the particle.

        TVirtualFitter* fitter; // A pointer to the ROOT fitter used for fitting.

        double kin_en; // The kinetic energy of the particle [eV].
        double e_err;  // The error of the kinetic energy [eV].

        /// @brief Gets a point on the computed trajectory at a specified index.
        /// @param index The index of the desired point on the trajectory.
        /// @return The RKPoint object representing the position (with time) of the particle at the specified index on the trajectory.
        const RKPoint& GetPoint(const int& index) const
        {
            using namespace constants;
            auto rk_pts = curr_rk->GetResults();
            return {1e+9*rk_pts[index].at(0,0),m2cm*rk_pts[index].at(1,0),m2cm*rk_pts[index].at(2,0),m2cm*rk_pts[index].at(3,0)};
        }

        /// @brief Computes the square of the distance between a point and a point on the computed trajectory at a specified index.
        /// @param index The index of the point on the trajectory to compare the distance to.
        /// @param point The Vector object representing the point to calculate the distance from.
        /// @return The square of the distance between the specified point and the point on the trajectory at the specified index.
        double SqDist(const int& index, const Vector& point)
        {
            RKPoint rk_point = GetPoint(index);
            return rk_point.AsVector().SqDist(point);
        }

        /// @brief Computes the scaled derivative of the square of the distance between a point and the point on the computed trajectory with respect to the index.
        /// @param index The index of the point on the trajectory to compute the derivative at.
        /// @param point The Vector object representing the point to calculate the derivative for.
        /// @return The derivative of the square of the distance between the specified point and the computed trajectory at the specified index with respect to the index.
        double DistDerivative(const int& index, const Vector& point)
        {
            return SqDist(index + 1, point) - SqDist(index,point);
        }

        /// @brief Calculates the squared distance between the given point and the curve.
        /// @param point The point to calculate the distance to.
        /// @return The squared distance between the point and the curve.
        double GetSqDist(const Vector& point)
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

        /// @brief Calculates the sum of squared distances for a set of fit data.
        /// @param npar The number of parameters to fit.
        /// @param gin The gradient of the fit parameters.
        /// @param sumsq The calculated sum of squared distances.
        /// @param par The parameters to fit.
        /// @param iflag The fitting mode flag.
        void SumSqDist(int& npar, double* gin, double& sumsq, double* par, int iflag)
        {
            kin_en = par[0];
            
            if (curr_rk != nullptr) delete curr_rk;
            curr_rk = GetTrackRK(*magfield,electron,step,kin_en,origin,orientation);
            curr_rk->Run();
            
            sumsq = 0;
            for (auto p : fit_data) sumsq += p.count*GetSqDist(p.AsVector());
        }
        
        /// @brief A static wrapper function used to get the correct function pointer type needed for ROOT. Simply calls lastfit->SumSqDist().
        /// @param npar The number of parameters to fit.
        /// @param gin The gradient of the fit parameters.
        /// @param sumsq The calculated sum of squared distances.
        /// @param par The parameters to fit.
        /// @param iflag The fitting mode flag.
        static void Eval(int& npar, double* gin, double& sumsq, double* par, int iflag)
        {
            lastfit->SumSqDist(npar,gin,sumsq,par,iflag);
        }

    public:
        /// @brief Constructor for the RKFit class.
        /// @param magfield A pointer to the magnetic data.
        /// @param electron A boolean value indicating if the particle is electron or positron.
        /// @param step A double representing the step size.
        /// @param origin A Vector object representing the origin.
        /// @param orientation A Vector object representing the orientation.
        /// @param fit_data A vector of RecoPoint objects representing the fit data.
        RKFit(Field<Vector>* magfield, const bool& electron, const double& step, const Vector& origin,
              const Vector& orientation, const std::vector<RecoPoint>& fit_data)
            : magfield(magfield), electron(electron), step(step), origin(origin), orientation(orientation), fit_data(fit_data)
        { lastfit = this; }

        /// @brief Sets the kinetic energy of the particle being tracked.
        /// @param kin_en The kinetic energy to set.
        void SetEnergy(const double& kin_en) { this->kin_en = kin_en; }

        /// @brief Sets up the TVirtualFitter for the circle fitting with the number of parameters and printout options.
        /// @param parameters Number of fitting parameters.
        /// @param print If true, the fit will printout the fitting status. If false, the fit will be silent.
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

        /// @brief Fits the Runge-Kutta object to the given data using the MINUIT fitter with specified maximum number of iterations and toleration
        /// @param max_iter The maximum number of iterations for the fitting procedure (default = 500)
        /// @param toleration The toleration level for the fitting procedure (default = 0.001)
        void FitRK(double max_iter = 500, double toleration = 0.001)
        {
            fitter->SetParameter(0,"kin_en",kin_en,1000,1000000,16000000);

            double arglist[2] = {max_iter,toleration};  // max iterations, step size (toleration)
            fitter->ExecuteCommand("MIGRAD",arglist,2); // last one num of prints (verbosity)

            kin_en = fitter->GetParameter(0);
            e_err  = fitter->GetParError(0);
        }

        /// @brief Prints the fitting parameters of the RK4 fit with their errors.
        void PrintFitParams()
        {
            std::cout << "\nRK FIT RESULT:\n";
            std::cout << "Kinetic energy:  " << kin_en << " +- " << e_err << "\n\n";
        }

        TGraph2D* GetFitGraph() {return GetGraphRK(curr_rk);}
    };

    RKFit* RKFit::lastfit = nullptr;
} // namespace X17
