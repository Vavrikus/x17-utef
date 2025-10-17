#pragma once

// ROOT dependencies
#include "TF1.h"
#include "TGraph.h"
#include "TSpline.h"

// X17 dependencies
#include "Field.h"

namespace X17
{
    /// @brief Function for simple circlular arc fit (using function for half of circle, curves up).
    /// @param graph TGraph object to be fitted.
    /// @param min Lower bound of the fit [cm].
    /// @param max Upper bound of the fit [cm].
    /// @return Fitted function.
    TF1* FitCircle(TGraph* graph, double min, double max);

    /// @brief Function for circular arc (upper half) with smoothly attached lines at endpoints (nodes).
    /// @param x Variable x [cm].
    /// @param par Set of parameters (radius of the circle, 1st and 2nd node x coordinates, 1st and 2nd node y coordinates).
    /// @return Value at x.
    double circle_func(double* x, double* par);

    /// @brief Function for fitting with circular arc with smoothly attached lines at endpoints.
    /// @param graph TGraph object to be fitted.
    /// @param min Lower bound of the fit [cm].
    /// @param max Upper bound of the fit [cm].
    /// @return Fitted function.
    TF1* FitCircle2(TGraph* graph, double min, double max);

    /// @brief Function for energy reconstruction from circle with lines fit.
    /// @param fit Fitted circle with lines function.
    /// @param magfield Magnetic data.
    /// @param magnetic Output graph for magnetic field along fitted trajectory.
    /// @param min Lower bound [cm].
    /// @param max Upper bound [cm].
    /// @param step Step between iterations.
    void RecoEnergy(TF1* fit, const Field<Vector>& magfield, TGraph* magnetic, double min, double max, double step);
}