// ROOT dependencies
#include "TF1.h"
#include "TGraph.h"
#include "TSpline.h"

// X17 dependencies
#include "Field.h"

namespace X17
{
    /// @brief A function for simple circlular arc fit (using function for half of circle, curves up).
    /// @param graph The TGraph object to be fitted.
    /// @param min The lower bound of the fit [cm].
    /// @param max The upper bound of the fit [cm].
    /// @return The fitted function.
    TF1* FitCircle(TGraph* graph, double min, double max);

    /// @brief A function for circular arc with smoothly attached lines at endpoints (nodes).
    /// @param x The variable x [cm].
    /// @param par The set of parameters (radius of the circle, 1st and 2nd node x and y coordinates).
    /// @return The value at x.
    double circle_func(double* x, double* par);

    /// @brief A function for fitting with circular arc with smoothly attached lines at endpoints.
    /// @param graph The TGraph object to be fitted.
    /// @param min The lower bound of the fit [cm].
    /// @param max The upper bound of the fit [cm].
    /// @return The fitted function.
    TF1* FitCircle2(TGraph* graph, double min, double max);

    /// @brief A function for energy reconstruction from spline fit.
    /// @param sp_fit Fitted splines.
    /// @param magfield The magnetic data.
    /// @param energy The output graph for reconstructed energy as function of coordinate.
    /// @param radius The output graph for reconstructed radius as function of coordinate.
    /// @param magnetic The output graph for magnetic field along fitted trajectory.
    /// @param min The lower bound [cm].
    /// @param max The upper bound [cm].
    /// @param step The step between iterations.
    void RecoEnergy(TSpline3* sp_fit, const Field<Vector>& magfield, TGraph* energy, TGraph* radius, TGraph* magnetic, double min, double max, double step);

    /// @brief A function for energy reconstruction from circle with lines fit.
    /// @param fit The fitted circle with lines function.
    /// @param magfield Magnetic data.
    /// @param magnetic The output graph for magnetic field along fitted trajectory.
    /// @param min The lower bound [cm].
    /// @param max The upper bound [cm].
    /// @param step The step between iterations.
    void RecoEnergy(TF1* fit, const Field<Vector>& magfield, TGraph* magnetic, double min, double max, double step);
}