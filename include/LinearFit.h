#pragma once

namespace X17
{
    /// @brief Class for performing a linear fit on a set of data points.
    /// @tparam N The dimension of the data points.
    template <int N>
    class LinearFit
    {
    private:
    public:
    };
    
} // namespace X17

#include <TDecompLU.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <iostream>

struct FPoint
{
    double e_sim;
    double theta;
    double varphi;
    double dev_reco;
};


int LinearFit3D(std::vector<FPoint> data)
{
    const int n = data.size();     // počet bodů
    const int p = 4;               // počet parametrů: 1 (bias) + 3 (x1,x2,x3)

    TMatrixD X(n, p);              // designová matice
    TVectorD Y(n);                 // hodnoty y

    for (int i = 0; i < n; ++i) {
        X(i, 0) = 1.0;             // bias (a0)
        X(i, 1) = data[i].e_sim;   // x1
        X(i, 2) = data[i].theta;   // x2
        X(i, 3) = data[i].varphi;  // x3
        Y[i] = data[i].dev_reco;   // y
    }

    // Normální rovnice: (A^T A) * a = A^T * Y
    TMatrixD Xt = TMatrixD(TMatrixD::kTransposed, X);
    TMatrixD XtX = Xt * X;
    TVectorD XtY = Xt * Y;

    // Vyřešení soustavy rovnic
    TDecompLU lu(XtX);
    TVectorD coeffs = XtY;
    Bool_t ok = lu.Solve(coeffs);

    if (!ok) {
        std::cerr << "Fit selhal!" << std::endl;
        return 1;
    }

    std::cout << "Fit parametry:" << std::endl;
    for (int i = 0; i < p; ++i) {
        std::cout << "a" << i << " = " << coeffs[i] << std::endl;
    }

    return 0;
}
