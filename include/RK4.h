#pragma once

// C++ dependencies
#include <functional>
#include <vector>

// ROOT dependencies
#include "TGraph2D.h"
#include "TVirtualFitter.h"

// X17 dependencies
#include "Matrix.h"
#include "Points.h"

namespace X17
{
    /// @brief A templated class implementing the Runge-Kutta method of order 4 for solving differential equations.
    /// @tparam N The dimension of the system of differential equations.
    template<int N>
    class RK4
    {
        typedef Matrix<N,1> VectorN;
        typedef std::function<void(double,const VectorN&,VectorN&)> VecFn;
        typedef std::function<bool(double,const VectorN&)> EndFn;

    private:
        VecFn m_dif_eq;    // The function object that returns the derivative of the dependent variable with respect to the independent variable.
        VectorN m_current; // The current value of the dependent variable.
        double m_param;    // The current value of the independent variable.
        double m_step;     // The step size for the integration.
        EndFn m_end_fn;    // The function object that determines the end condition for the integration.

        std::vector<VectorN> m_results; // A vector that stores the values of the dependent variable at different points in the integration.
        
    public:
        /// @brief Constructs an RK4 object with the specified initial conditions and integration parameters.
        /// @param start The initial value of the independent variable.
        /// @param step The step size for the integration.
        /// @param equation The function object that returns the derivative of the dependent variable with respect to the independent variable.
        /// @param initial The initial value of the dependent variable.
        /// @param end_condition The function object that determines the end condition for the integration.
        RK4(double start, double step, VecFn equation, VectorN initial, EndFn end_condition);

        /// @brief Performs the integration using the fourth-order Runge-Kutta method.
        void Integrate();

        /// @brief Returns the vector of values of the dependent variable at different points in the integration.
        /// @return The vector of values of the dependent variable at different points in the integration.
        const std::vector<VectorN>& GetResults() const { return m_results; }

        /// @brief Returns the number of points in the integration.
        /// @return The number of points in the integration.
        int GetSize() const { return m_results.size(); }

        /// @brief Returns the step size used for the integration.
        /// @return The step size used for the integration.
        double GetStep() const { return m_step; }
    };

    /// @brief Computes the initial parameters for a simulation using the RK4 algorithm.
    /// @param kin_en      Kinetic energy of the particle [eV].
    /// @param origin      Coordinates of the origin [cm].
    /// @param orientation Initial direction of motion.
    /// @return A matrix containing the position and four-velocity of the particle.
    Matrix<8,1> GetInitParams(double kin_en, Vector origin, Vector orientation);

    /// @brief Calculates the electromagnetic field tensor at a given position.
    /// @param magfield Vector field representing the magnetic field.
    /// @param position The position at which to calculate the electromagnetic field [m].
    /// @return A 4x4 matrix representing the electromagnetic field tensor.
    Matrix<4,4> GetEMtensor(const Field<Vector>& magfield, Vector position);

    /// @brief Calculates the electromagnetic motion of a charged particle using the Lorentz force equation.
    /// @param magfield Reference to the magnetic field.
    /// @param electron Boolean value indicating whether the particle is an electron or a positron.
    /// @param tau The independent variable (proper time) [s].
    /// @param params Vector containing the initial position and four-velocity of the particle.
    /// @param output Vector containing the final position and four-velocity of the particle.
    void EMMotion(const Field<Vector>& magfield, bool electron, double tau, const Matrix<8,1>& params, Matrix<8,1>& output);

    /// @brief Checks if the particle is out of the first sector TPC after it has moved beyond the minimum x position.
    /// @param tau Proper time elapsed (independent variable) [s].
    /// @param params Vector of position and four-velocity of the particle.
    /// @return True if the particle is out of the TPC after it has moved beyond the minimum x position, false otherwise.
    bool IsOutOfTPC(double tau, const Matrix<8,1>& params);

    /// @brief Checks if the particle is out of the first sector TPC and outside a rectangular volume in front of it (x < xmin).
    /// @param tau Proper time elapsed (independent variable) [s].
    /// @param params Vector of position and four-velocity of the particle.
    /// @return True if the particle is out of the TPC and outside the rectangular volume, false otherwise.
    bool IsOutOfTPCandRect(double tau, const Matrix<8,1>& params);

    /// @brief Returns a pointer to an instance of RK4 class for tracking a charged particle in a magnetic field.
    /// @param magfield Vector field representing the magnetic field.
    /// @param electron Boolean indicating if the particle is an electron or a positron.
    /// @param step Step size for the RK4 integration [s].
    /// @param kin_en Kinetic energy of the particle [eV].
    /// @param origin Coordinates of the origin [cm].
    /// @param orientation Initial direction of motion.
    /// @param big_volume Should a volume be added in front of the TPC?
    /// @return Pointer to an instance of RK4 class.
    RK4<8>* GetTrackRK(const Field<Vector>& magfield, bool electron, double step, double kin_en, Vector origin, Vector orientation, bool big_volume = false);

    /// @brief Gets a point on the computed trajectory at a specified index.
    /// @param track Pointer to a simulated RK4 trajectory.
    /// @param index The index of the desired point on the trajectory.
    /// @return The RKPoint object representing the position (with time) of the particle at the specified index on the trajectory.
    RKPoint GetTrackRKPoint(const RK4<8>* track, int index);

    /// @brief Computes the square of the distance between a point and a point on the computed trajectory at a specified index.
    /// @param track Pointer to a simulated RK4 trajectory.
    /// @param index The index of the point on the trajectory to compare the distance to.
    /// @param point The Vector object representing the point to calculate the distance from.
    /// @return The square of the distance between the specified point and the point on the trajectory at the specified index.
    double GetTrackRKSqDist(const RK4<8>* track, int index, Vector point);

    /// @brief Computes the scaled derivative of the square of the distance between a point and the point on the computed trajectory with respect to the index.
    /// @param track Pointer to a simulated RK4 trajectory.
    /// @param index The index of the point on the trajectory to compute the derivative at.
    /// @param point The Vector object representing the point to calculate the derivative for.
    /// @return The derivative of the square of the distance between the specified point and the computed trajectory at the specified index with respect to the index.
    inline double TrackRKDistDer(const RK4<8>* track, int index, Vector point)
    {
        return GetTrackRKSqDist(track, index + 1, point) - GetTrackRKSqDist(track, index, point);
    }

    /// @brief Calculates the squared distance between the given point and the curve.
    /// @param track Pointer to a simulated RK4 trajectory.
    /// @param point The point to calculate the distance to.
    /// @return The squared distance between the point and the curve.
    double GetTrackRKSqDist(const RK4<8>* track, Vector point);

    /// @brief Calculates the squared distance between the given point and the curve and returns the closest point on the curve.
    /// @param track Pointer to a simulated RK4 trajectory.
    /// @param point The point to calculate the distance to.
    /// @param closest_point Will be set to the closest point on the curve to the given point.
    /// @return The squared distance between the point and the curve.
    double GetTrackRKSqDistAndCP(const RK4<8>* track, Vector point, Vector& closest_point);

    /// @brief Returns a TGraph2D object representing the track of a particle obtained from a Runge-Kutta 4th order simulation.
    /// @param rk_track A pointer to a RK4 object containing the simulation results.
    /// @return A pointer to a TGraph2D object representing the particle track.
    TGraph2D* GetGraphRK(const RK4<8>* rk_track);

    /// @brief A class for fitting a 4th order Runge-Kutta generated points to a particle's trajectory.
    class RKFit
    {
    private:
        static RKFit* lastfit; // A static pointer to the last instance of RKFit class created.

        Field<Vector>* m_magfield; // A pointer to the magnetic field that the particle is traveling through.
        bool m_electron;           // A boolean indicating whether the particle is an electron (true) or a positron (false).
        double m_step;             // The step size used in the Runge-Kutta 4th order method.
        Vector m_origin;           // The initial position of the particle.
        Vector m_orientation;      // The initial orientation of the particle.
        bool m_big_volume;         // A boolean indicating if a volume should be added in front of the TPC.

        std::vector<RecoPoint> m_fit_data; // A vector containing the reconstructed points used for fitting.
        RK4<8>* m_curr_rk = nullptr;       // A pointer to the RK4 object used to compute the trajectory of the particle.

        TVirtualFitter* m_fitter; // A pointer to the ROOT fitter used for fitting.

        double m_kin_en; // The kinetic energy of the particle [eV].
        double m_e_err;  // The error of the kinetic energy [eV].

    public:
        /// @brief Constructor for the RKFit class.
        /// @param magfield A pointer to the magnetic data.
        /// @param electron A boolean value indicating if the particle is electron or positron.
        /// @param step A double representing the step size.
        /// @param origin A Vector object representing the origin.
        /// @param orientation A Vector object representing the orientation.
        /// @param fit_data A vector of RecoPoint objects representing the fit data.
        /// @param big_volume A boolean value indicating if a volume should be added in front of the TPC.
        RKFit(Field<Vector>* magfield, bool electron, double step, Vector origin, Vector orientation, const std::vector<RecoPoint>& fit_data, bool big_volume = false);

        /// @brief Getter for the kinetic energy.
        /// @return Kinetic energy [eV].
        double GetEnergy() { return m_kin_en; }

        /// @brief Getter for the kinetic energy error.
        /// @return Error of the kinetic energy [eV].
        double GetEnergyError() { return m_e_err; }

        /// @brief Sets the kinetic energy of the particle being tracked.
        /// @param kin_en The kinetic energy to set.
        void SetEnergy(double kin_en) { this->m_kin_en = kin_en; }

        /// @brief Sets up the TVirtualFitter for the circle fitting with the number of parameters and printout options.
        /// @param parameters Number of fitting parameters.
        /// @param print If true, the fit will printout the fitting status. If false, the fit will be silent.
        void SetFitter(int parameters = 1, bool print = true);

        /// @brief Fits the Runge-Kutta object to the given data using the MINUIT fitter with specified maximum number of iterations and toleration.
        /// @param max_iter The maximum number of iterations for the fitting procedure (default = 500).
        /// @param toleration The toleration level for the fitting procedure (default = 0.001).
        void FitRK(double max_iter = 500, double toleration = 0.001);

        /// @brief Prints the fitting parameters of the RK4 fit with their errors.
        void PrintFitParams() const;

        TGraph2D* GetFitGraph() { return GetGraphRK(m_curr_rk); }

    private:
        /// @brief Calculates the sum of squared distances for a set of fit data.
        /// @param npar The number of parameters to fit.
        /// @param gin The gradient of the fit parameters.
        /// @param sumsq The calculated sum of squared distances.
        /// @param par The parameters to fit.
        /// @param iflag The fitting mode flag.
        void _SumSqDist(int& npar, double* gin, double& sumsq, double* par, int iflag);
        
        /// @brief A static wrapper function used to get the correct function pointer type needed for ROOT. Simply calls lastfit->_SumSqDist().
        /// @param npar The number of parameters to fit.
        /// @param gin The gradient of the fit parameters.
        /// @param sumsq The calculated sum of squared distances.
        /// @param par The parameters to fit.
        /// @param iflag The fitting mode flag.
        static void _Eval(int& npar, double* gin, double& sumsq, double* par, int iflag)
        {
            lastfit->_SumSqDist(npar,gin,sumsq,par,iflag);
        }
    };

    inline RKFit* RKFit::lastfit = nullptr;
} // namespace X17

// Templated function definitions.
#include "RK4.inl"