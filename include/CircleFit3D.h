#pragma once

// ROOT dependencies
#include "TGraph.h"
#include "TGraph2D.h"
#include "TVirtualFitter.h"

// X17 dependencies
#include "Field.h"
#include "Points.h"
#include "Utilities.h"
#include "Vector.h"
#include "X17Utilities.h"

namespace X17
{
    /// @brief class for fitting reconstructed track with circular arc with smoothly attached lines
    class CircleFit3D
    {
        typedef void(*EvalFn)(Int_t&, Double_t*, Double_t&, Double_t*, Int_t);

    private:
        static CircleFit3D* lastfit; // Pointer to the last instance.

        std::vector<RecoPoint> fit_data; // A std::vector of reconstructed points to be fitted.
        Vector origin;                   // The origin of the first line [cm].
        Vector orientation;              // The orientation vector of the first line. Used for theta and phi calculation.
        double length;                   // The length of the first line (also max parameter value) [cm].
        double alpha;                    // The rotation angle of the arc around the first line (if 0 curves towards negative z) [rad].
        double radius;                   // The radius of the arc [rad].
        double phi_max;                  // The maximal angle on the arc [rad].

        double cos_theta; // The cosine of theta.
        double sin_theta; // The sine of theta.
        double cos_phi;   // The cosine of phi.
        double sin_phi;   // The sine of phi.
        double cos_alpha; // The cosine of alpha.
        double sin_alpha; // The sine of alpha.

        Vector origin2;      // The origin of the second line (to be calculated).
        Vector orientation2; // The orientation vector of the second line (to be calculated).

        Vector originc; // The origin of the circular arc (to be calculated).
        Vector center;  // The center of the circular arc (to be calculated).
        Vector normal;  // The normal to the plane of the circular arc (to be calculated).

        TVirtualFitter* gFitter; // Contains TVirtualFitter instance used for the fit.
        double l_err;            // The length fit error.
        double a_err;            // The alpha fit error.
        double r_err;            // The radius fit error.
        double phi_err;          // The phi fit error.

        /// @brief Returns a point on the line.
        /// @param param The parameter that defines the position on the line.
        /// @param first_line Whether to return a point on the first line or the second line.
        /// @return A point on the specified line.
        Vector GetLinePoint(double param, bool first_line);

        /// @brief Calculates a point on the circle described by the current instance.
        /// @param varphi The angle (in radians) at which to calculate the point. Zero coresponds to the beginning of the arc.
        /// @return A vector representing the point on the circle.
        Vector GetCirclePoint(double varphi);

        /// @brief Updates the curve parameters based on the current orientation, length, radius, alpha, and phi_max values.
        /// This function calculates the new values of cos_theta, sin_theta, cos_phi, sin_phi, cos_alpha, sin_alpha, originc,
        /// normal, center, origin2, and orientation2. These values are used to update the properties of the curve.
        /// @note This function assumes that the orientation, length, radius, alpha, and phi_max values have already been set.
        void UpdateCurve();

        /// @brief Calculates the squared distance between a given point and a line.
        /// @param point The point to calculate the distance from.
        /// @param first_line Flag indicating whether to use the first line (true) or the second line (false).
        /// @return The squared distance between the point and the line.
        double LineSqDist(const RecoPoint& point, bool first_line);

        /// @brief Calculate the squared distance between a given point and the circular trajectory.
        /// @param point The point to calculate the distance from.
        /// @return The squared distance between the point and the circular trajectory.
        double CircleSqDist(const RecoPoint& point);

        /// @brief Computes the squared distance between the given point and the curve.
        /// @param point The point to compute the squared distance from.
        /// @return The squared distance between the point and the curve.
        double SqDistance(const RecoPoint& point);

        /// @brief Calculates the sum of squared distances from each point in the fit_data vector to the curve. Each distance is scaled by the count/charge.
        /// @return The sum of squared distances from each point in fit_data to the curve.
        double SumSq();

        /// @brief Calculate the sum of squared distances for a given set of parameters.
        /// @param npar The number of parameters.
        /// @param gin An array of size npar used to store the first derivatives of the function with respect to the parameters.
        /// @param sumsq The sum of squared distances, which is calculated by this function.
        /// @param par An array of size npar that contains the values of the parameters.
        /// @param iflag An integer that can be used to control the behavior of the function.
        void EvalSumSq(int& npar, double* gin, double& sumsq, double* par, int iflag);

        /// @brief A static wrapper function used to get the correct function pointer type needed for ROOT. Simply calls lastfit->EvalSumSq().
        /// @param npar The number of parameters.
        /// @param gin An array of size npar used to store the first derivatives of the function with respect to the parameters.
        /// @param sumsq The sum of squared distances, which is calculated by this function.
        /// @param par An array of size npar that contains the values of the parameters.
        /// @param iflag An integer that can be used to control the behavior of the function.
        static void Eval(int& npar, double* gin, double& sumsq, double* par, int iflag);

        /// @brief Computes the average magnetic field vector over the circular path defined by the
        ///        current object's parameters, using the provided magnetic field vector field.
        /// @param magfield Pointer to an object containing the magnetic field information.
        /// @param step Step size to use when sampling the circular path. Default value is 0.1.
        /// @return The average magnetic field vector over the circular path.
        Vector GetAvgField(const Field<Vector>& magfield, double step = 0.1);

        /// @brief Get the magnetic field at the point on the circle that passes through the middle of the X17 detector in the x direction.
        /// @param magfield The magnetic field vector field to sample from.
        /// @param tolerance The maximum deviation from the middle of the detector in the x direction for which to stop the binary search. Default is 0.0001.
        /// @return The magnetic field vector at the point on the circle that passes through the middle of the X17 detector in the x direction.
        Vector GetMiddleField(const Field<Vector>& magfield, double tolerance = 0.0001);

    public:
        /// @brief The default constructor for the 3D circle fitting algorithm that smoothly attaches lines.
        CircleFit3D() { lastfit = this; }

        /// @brief Constructor for the 3D circle fitting algorithm that smoothly attaches lines.
        /// @param orig The origin of the first line [cm].
        /// @param orient The orientation vector of the first line.
        CircleFit3D(const Vector& orig, const Vector& orient);

        /// @brief Deleted copy constructor. 
        CircleFit3D(const CircleFit3D&) = delete;

        /// @brief Deleted assignment operator. 
        CircleFit3D& operator=(const CircleFit3D&) = delete;

        /// @brief Get the current fit data points.
        /// @return The vector of RecoPoints used for fitting.
        std::vector<RecoPoint> GetData() { return fit_data; }

        /// @brief Add a new RecoPoint to the fit data.
        /// @param x The x-coordinate of the RecoPoint.
        /// @param y The y-coordinate of the RecoPoint.
        /// @param z The z-coordinate of the RecoPoint.
        /// @param count The count associated with the RecoPoint.
        void AddPoint(double x, double y, double z, int count) { fit_data.emplace_back(x,y,z,count); }

        /// @brief Add a new RecoPoint to the fit data.
        /// @param p The RecoPoint to add.
        void AddPoint(RecoPoint p) { fit_data.push_back(p); }

        /// @brief Add a new RecoPoint to the fit data.
        /// @param p The RKPoint with necessary coordinate information.
        void AddPoint(RKPoint p) { fit_data.emplace_back(p.x,p.y,p.z,1); }

        /// @brief Set the alpha angle for the circle fit.
        /// @param electron If true, set alpha to 0; if false, set alpha to pi.
        void SetAlpha(bool electron)
        {
            if (electron) alpha = 0;
            else          alpha = M_PI;
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
        void PrintFitParams();

        /// @brief Generates a TGraph2D representing the fitted circle with attached lines.
        /// @param step The step size between points in the TGraph2D.
        /// @param dist The distance from the TPC walls within which the points are included.
        /// @return A TGraph2D representing the fitted circle with attached lines.
        TGraph2D* GetGraph(double step = 0.1, double dist = 0);

        /// @brief Calculates the energy of a charged particle moving along the 3D circle fit.
        /// @param magfield The magnetic field data used for the calculation.
        /// @param middle Boolean flag indicating whether to use the middle field of the circle (true) or the average field (false).
        /// @return The kinetic energy of the particle.
        double GetEnergy(const Field<Vector>& magfield, bool middle = true);

        /// @brief Returns a graph representing the energy along the track in the given magnetic field.
        /// @param magfield The magnetic field data used for the calculation.
        /// @param step The step size for the parameter along the track.
        /// @return A TGraph object representing the energy as a function of the parameter.
        TGraph* GetEnergyGraph(const Field<Vector>& magfield, double step = 0.1);
    };
} // namespace X17