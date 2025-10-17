#pragma once

// C++ dependencies
#include <iostream>

// ROOT dependencies
#include "TDecompLU.h"
#include "TMatrixD.h"
#include "TVectorD.h"

struct FPoint
{
    bool electron;
    double e_sim;
    double e_rec;
    double theta;
    double varphi;
    double dev_reco;
};


int LinearFit3D(const std::vector<FPoint>& data, const int p = 4, bool use_reco = false)
{
    const int n = data.size();     // number of points
    // const int p = 4;            // number of parameters: 1 (bias) + 3 (x1,x2,x3)
    if (p < 2 || p > 4) return -1;
    if (n < p) return -2;

    TMatrixD X(n, p);              // design matrix
    TVectorD Y(n);                 // y values

    for (int i = 0; i < n; ++i)
    {
        X(i, 0) = 1.0;                                      // bias (a0)
        X(i, 1) = use_reco ? data[i].e_rec : data[i].e_sim; // x1
        if (p > 2) X(i, 2) = data[i].theta;                 // x2
        if (p > 3) X(i, 3) = data[i].varphi;                // x3
        Y[i] = data[i].dev_reco;                            // y
    }

    // Solving (X^T X) * a = X^T * Y
    TMatrixD Xt = TMatrixD(TMatrixD::kTransposed, X);
    TMatrixD XtX = Xt * X;
    TVectorD XtY = Xt * Y;

    // Solving the linear system
    TDecompLU lu(XtX);
    TVectorD coeffs = XtY;
    Bool_t ok = lu.Solve(coeffs);

    if (!ok)
    {
        std::cerr << "Fit failed!" << std::endl;
        return 1;
    }

    std::cout << "Fit parameters:" << std::endl;
    for (int i = 0; i < p; ++i)
        std::cout << "a" << i << " = " << coeffs[i] << std::endl;

    return 0;
}
