#pragma once

// C++ dependencies
#include <vector>

// ROOT dependencies
#include "TGraph.h"
#include "TGraph2D.h"
#include "TVirtualFitter.h"

// X17 dependencies
#include "Field.h"
#include "Points.h"
#include "Vector.h"

namespace X17
{
    /// @brief Class for fitting of the reconstructed track with circular arc with smoothly attached lines.
    class CircleFit3D
    {
        typedef void(*EvalFn)(Int_t&, Double_t*, Double_t&, Double_t*, Int_t);

    private:
        static CircleFit3D* lastfit; // Pointer to the last instance.

        std::vector<RecoPoint> m_fit_data; // A std::vector of reconstructed points to be fitted.
        Vector m_origin;                   // The origin of the first line [cm].
        Vector m_orientation;              // The orientation vector of the first line. Used for theta and phi calculation.
        double m_length;                   // The length of the first line (also max parameter value) [cm].
        double m_alpha;                    // The rotation angle of the arc around the first line (if 0 curves towards negative z) [rad].
        double m_radius;                   // The radius of the arc [rad].
        double m_phi_max;                  // The maximal angle on the arc [rad].

        double m_cos_theta; // The cosine of theta.
        double m_sin_theta; // The sine of theta.
        double m_cos_phi;   // The cosine of phi.
        double m_sin_phi;   // The sine of phi.
        double m_cos_alpha; // The cosine of m_alpha.
        double m_sin_alpha; // The sine of m_alpha.

        Vector m_originc; // The origin of the circular arc (to be calculated).
        Vector m_center;  // The center of the circular arc (to be calculated).
        Vector m_normal;  // The normal to the plane of the circular arc (to be calculated).

        Vector m_origin2;      // The origin of the second line (to be calculated).
        Vector m_orientation2; // The orientation vector of the second line (to be calculated).

        TVirtualFitter* m_fitter;  // Contains TVirtualFitter instance used for the fit.
        double m_l_err;            // The length fit error.
        double m_a_err;            // The m_alpha fit error.
        double m_r_err;            // The m_radius fit error.
        double m_phi_err;          // The phi fit error.

    public:
        /// @brief The default constructor for the 3D circle fitting algorithm that smoothly attaches lines.
        CircleFit3D() { lastfit = this; }

        /// @brief Constructor for the 3D circle fitting algorithm that smoothly attaches lines.
        /// @param orig The origin of the first line [cm].
        /// @param orient The orientation vector of the first line.
        CircleFit3D(Vector orig, Vector orient);

        /// @brief Deleted copy constructor. 
        CircleFit3D(const CircleFit3D&) = delete;

        /// @brief Deleted assignment operator. 
        CircleFit3D& operator=(const CircleFit3D&) = delete;

        /// @brief Get the current fit data points.
        /// @return The vector of RecoPoints used for fitting.
        std::vector<RecoPoint> GetData() const { return m_fit_data; }

        /// @brief Add a new RecoPoint to the fit data.
        /// @param x The x-coordinate of the RecoPoint.
        /// @param y The y-coordinate of the RecoPoint.
        /// @param z The z-coordinate of the RecoPoint.
        /// @param count The count associated with the RecoPoint.
        void AddPoint(double x, double y, double z, int count) { m_fit_data.emplace_back(x,y,z,count); }

        /// @brief Add a new RecoPoint to the fit data.
        /// @param p The RecoPoint to add.
        void AddPoint(RecoPoint p) { m_fit_data.push_back(p); }

        /// @brief Add a new RecoPoint to the fit data.
        /// @param p The RKPoint with necessary coordinate information.
        void AddPoint(RKPoint p) { m_fit_data.emplace_back(p.x,p.y,p.z,1); }

        /// @brief Set the m_alpha angle for the circle fit.
        /// @param electron If true, set m_alpha to 0; if false, set m_alpha to pi.
        void SetAlpha(bool electron)
        {
            if (electron) m_alpha = 0;
            else          m_alpha = M_PI;
        }

        /// @brief Sets up the TVirtualFitter for the circle fitting with the number of parameters and printout options.
        /// @param parameters Number of fitting parameters.
        /// @param print If true, the fit will printout the fitting status. If false, the fit will be silent.
        void SetFitter(int parameters = 4, bool print = true);

        /// @brief Fits a 3D circle using a maximum of max_iter iterations and the given toleration.
        /// @param max_iter The maximum number of iterations for the fit. Defaults to 500.
        /// @param toleration The toleration for the fit. Defaults to 0.001.
        void FitCircle3D(double max_iter = 500, double toleration = 0.001);

        /// @brief Prints the fitting parameters of the circle with their errors.
        void PrintFitParams() const;

        /// @brief Generates a TGraph2D representing the fitted circle with attached lines.
        /// @param step The step size between points in the TGraph2D.
        /// @param dist The distance from the TPC walls within which the points are included.
        /// @return A TGraph2D representing the fitted circle with attached lines.
        TGraph2D* GetGraph(double step = 0.1, double dist = 0) const;

        /// @brief Calculates the energy of a charged particle moving along the 3D circle fit.
        /// @param magfield The magnetic field data used for the calculation.
        /// @param middle Boolean flag indicating whether to use the middle field of the circle (true) or the average field (false).
        /// @return The kinetic energy of the particle.
        double GetEnergy(const Field<Vector>& magfield, bool middle = true) const;

        /// @brief Returns a graph representing the energy along the track in the given magnetic field.
        /// @param magfield The magnetic field data used for the calculation.
        /// @param step The step size for the parameter along the track.
        /// @return A TGraph object representing the energy as a function of the parameter.
        TGraph* GetEnergyGraph(const Field<Vector>& magfield, double step = 0.1) const;

    private:
        /// @brief Returns a point on the line.
        /// @param param The parameter that defines the position on the line.
        /// @param first_line Whether to return a point on the first line or the second line.
        /// @return A point on the specified line.
        Vector _GetLinePoint(double param, bool first_line) const;

        /// @brief Calculates a point on the circle described by the current instance.
        /// @param varphi The angle (in radians) at which to calculate the point. Zero coresponds to the beginning of the arc.
        /// @return A vector representing the point on the circle.
        Vector _GetCirclePoint(double varphi) const;

        /// @brief Updates the curve parameters based on the current m_orientation, m_length, m_radius, m_alpha, and m_phi_max values.
        /// This function calculates the new values of m_cos_theta, m_sin_theta, m_cos_phi, m_sin_phi, m_cos_alpha, m_sin_alpha, m_originc,
        /// m_normal, m_center, m_origin2, and m_orientation2. These values are used to update the properties of the curve.
        /// @note This function assumes that the m_orientation, m_length, m_radius, m_alpha, and m_phi_max values have already been set.
        void _UpdateCurve();

        /// @brief Calculates the squared distance between a given point and a line.
        /// @param point The point to calculate the distance from.
        /// @param first_line Flag indicating whether to use the first line (true) or the second line (false).
        /// @return The squared distance between the point and the line.
        double _LineSqDist(RecoPoint point, bool first_line) const;

        /// @brief Calculate the squared distance between a given point and the circular trajectory.
        /// @param point The point to calculate the distance from.
        /// @return The squared distance between the point and the circular trajectory.
        double _CircleSqDist(RecoPoint point) const;

        /// @brief Computes the squared distance between the given point and the curve.
        /// @param point The point to compute the squared distance from.
        /// @return The squared distance between the point and the curve.
        double _SqDistance(RecoPoint point) const;

        /// @brief Calculates the sum of squared distances from each point in the m_fit_data vector to the curve. Each distance is scaled by the count/charge.
        /// @return The sum of squared distances from each point in m_fit_data to the curve.
        double _SumSq() const;

        /// @brief Calculate the sum of squared distances for a given set of parameters.
        /// @param npar The number of parameters.
        /// @param gin An array of size npar used to store the first derivatives of the function with respect to the parameters.
        /// @param sumsq The sum of squared distances, which is calculated by this function.
        /// @param par An array of size npar that contains the values of the parameters.
        /// @param iflag An integer that can be used to control the behavior of the function.
        void _EvalSumSq(int& npar, double* gin, double& sumsq, double* par, int iflag);

        /// @brief A static wrapper function used to get the correct function pointer type needed for ROOT. Simply calls lastfit->_EvalSumSq().
        /// @param npar The number of parameters.
        /// @param gin An array of size npar used to store the first derivatives of the function with respect to the parameters.
        /// @param sumsq The sum of squared distances, which is calculated by this function.
        /// @param par An array of size npar that contains the values of the parameters.
        /// @param iflag An integer that can be used to control the behavior of the function.
        static void _Eval(int& npar, double* gin, double& sumsq, double* par, int iflag);

        /// @brief Computes the average magnetic field vector over the circular path defined by the
        ///        current object's parameters, using the provided magnetic field vector field.
        /// @param magfield Pointer to an object containing the magnetic field information.
        /// @param step Step size to use when sampling the circular path. Default value is 0.1.
        /// @return The average magnetic field vector over the circular path.
        Vector _GetAvgField(const Field<Vector>& magfield, double step = 0.1) const;

        /// @brief Get the magnetic field at the point on the circle that passes through the middle of the X17 detector in the x direction.
        /// @param magfield The magnetic field vector field to sample from.
        /// @param tolerance The maximum deviation from the middle of the detector in the x direction for which to stop the binary search. Default is 0.0001.
        /// @return The magnetic field vector at the point on the circle that passes through the middle of the X17 detector in the x direction.
        Vector _GetMiddleField(const Field<Vector>& magfield, double tolerance = 0.0001) const;
    };
} // namespace X17