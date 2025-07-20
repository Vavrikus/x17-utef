// ROOT dependencies
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

// X17 dependencies
#include "EndPoint.h"
#include "Utilities.h"

ClassImp(X17::EndPoint)

namespace X17
{
    
    Matrix<4, 4> CovarianceMatrix(const std::vector<EndPoint>& points, EndPoint average)
    {
        typedef Matrix<4,1> Vector4;
        std::vector<Vector4> vectors;
        vectors.reserve(points.size());
        for (EndPoint p : points)
        {
            p -= average;
            vectors.push_back(p.ToVector4());
        }

        Matrix<4,4> cov_mat(0);
        for (const Vector4& v : vectors)
            for (int i = 0; i < 4; i++)
                for (int j = i; j < 4; j++)
                    cov_mat.at(i,j) += v.at(i,0) * v.at(j,0);

        cov_mat /= vectors.size() - 1;

        for(int i = 0; i < 4; i++)
            for (int j = i; j < 4; j++)
                cov_mat.at(j,i) = cov_mat.at(i,j);

        return cov_mat;
    }

    std::vector<double> SqMahalanobis(const std::vector<EndPoint>& values, bool no_z)
    {
        EndPoint average = GetAverage(values);
        Matrix<4, 4> cov = CovarianceMatrix(values, average);

        std::vector<double> sq_mahal;
        sq_mahal.reserve(values.size());

        if(no_z)
        {
            TMatrixD cov_m(3,3);
            cov_m(0,0) = cov.at(0,0); cov_m(0,1) = cov.at(0,1); cov_m(0,2) = cov.at(0,3);
            cov_m(1,0) = cov.at(1,0); cov_m(1,1) = cov.at(1,1); cov_m(1,2) = cov.at(1,3);
            cov_m(2,0) = cov.at(3,0); cov_m(2,1) = cov.at(3,1); cov_m(2,2) = cov.at(3,3);

            TMatrixD cov_inv(TMatrixD::kInverted,cov_m);

            for (const EndPoint& p : values)
            {
                TVectorD tv(3);
                tv(0) = p.x() - average.x();
                tv(1) = p.y() - average.y();
                tv(2) = p.t - average.t;
                sq_mahal.push_back(tv * (cov_inv * tv));
            }
        }
        else
        {
            TMatrixD cov_m(4,4,cov.elements);
            TMatrixD cov_inv(TMatrixD::kInverted,cov_m);

            for (const EndPoint& p : values)
            {
                Matrix<4,1> v = p.ToVector4() - average.ToVector4();
                TVectorD tv(4,v.elements);
                sq_mahal.push_back(tv * (cov_inv * tv));
            }
        }

        return sq_mahal;
    }

    void FillQQplot(TGraph* graph, const std::vector<EndPoint>& values, bool no_z)
    {
        int degrees_of_freedom = no_z ? 3 : 4;

        std::vector<double> sq_mahal = SqMahalanobis(values,no_z);
        std::sort(sq_mahal.begin(),sq_mahal.end());
        
        for (int i = 0; i < sq_mahal.size(); i++)
        {
            double chi2_quantile = TMath::ChisquareQuantile((i+1.0)/(sq_mahal.size()+1),degrees_of_freedom);
            graph->AddPoint(chi2_quantile,sq_mahal[i]);
        }
    }

    void MardiaTest(const std::vector<EndPoint>& values, double &outA, double &outB, bool no_z)
    {
        int dof = no_z ? 3 : 4; // degrees of freedom
        int d_fact = dof * (dof + 2);
        int N = values.size();

        // Hanusz et al. correction for the statistic
        double cor1 = (N - 1.0) / (N + 1.0);
        double cor2 = N * (N - 3.0) * (N - dof - 1.0) * (N - dof + 1.0) / ((N + 1.0) * (N + 1.0) * (N + 3.0) * (N + 5.0));

        EndPoint average = GetAverage(values);
        Matrix<4, 4> cov = CovarianceMatrix(values, average);

        double sumA = 0, sumB = 0;

        if(no_z)
        {
            TMatrixD cov_m(3,3);
            cov_m(0,0) = cov.at(0,0); cov_m(0,1) = cov.at(0,1); cov_m(0,2) = cov.at(0,3);
            cov_m(1,0) = cov.at(1,0); cov_m(1,1) = cov.at(1,1); cov_m(1,2) = cov.at(1,3);
            cov_m(2,0) = cov.at(3,0); cov_m(2,1) = cov.at(3,1); cov_m(2,2) = cov.at(3,3);

            // Use the biased maximum likelihood covariance matrix
            cov_m *= (N-1.0)/N;

            TMatrixD cov_inv(TMatrixD::kInverted,cov_m);

            for (EndPoint p : values)
            {
                TVectorD tv(3);
                tv(0) = p.x() - average.x();
                tv(1) = p.y() - average.y();
                tv(2) = p.t - average.t;
                sumB += std::pow(tv * (cov_inv * tv),2);

                for (EndPoint p2 : values)
                {
                    TVectorD tv2(3);
                    tv2(0) = p2.x() - average.x();
                    tv2(1) = p2.y() - average.y();
                    tv2(2) = p2.t - average.t;
                    sumA += std::pow(tv * (cov_inv * tv2),3);
                }
            }
        }
        else assert(false);

        outA = sumA / (6.0*N);
        outB = std::sqrt(N/(8.0*cor2*d_fact))*(sumB/N - cor1*d_fact);

        //std::cout << "Mardia's test: " << "A: " << outA << ", B: " << outB << "\n";
    }
}