// C++ dependencies
#include <iostream>
#include <vector>

// ROOT dependencies
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

// X17 dependencies
#include "Field.h"
#include "Points.h"
#include "X17Utilities.h"

namespace X17
{
    //// Public MicroPoint methods.

    void MicroPoint::MakeTTreeBranches(TTree* tree)
    {
        tree->Branch("x0",&start.point.x);
        tree->Branch("y0",&start.point.y);
        tree->Branch("z0",&start.point.z);
        tree->Branch("t0",&start.t);
        tree->Branch("e0",&e0);
        tree->Branch("x1",&end.point.x);
        tree->Branch("y1",&end.point.y);
        tree->Branch("z1",&end.point.z);
        tree->Branch("t1",&end.t);
        tree->Branch("e1",&e1);
    }

    void MicroPoint::SetTChainBranches(TChain* chain, bool old_data)
    {
        if (old_data)
        {
            std::cout << "MAKE SURE TO CHANGE THE COORDINATE ASSIGNMENT IF USING NEW SIMULATIONS!!!!!!!\n";
            chain->SetBranchAddress("z0",&start.point.x);
            chain->SetBranchAddress("x0",&start.point.y);
            chain->SetBranchAddress("y0",&start.point.z);
            chain->SetBranchAddress("t0",&start.t);
            chain->SetBranchAddress("e0",&e0);
            chain->SetBranchAddress("z1",&end.point.x);
            chain->SetBranchAddress("x1",&end.point.y);
            chain->SetBranchAddress("y1",&end.point.z);
            chain->SetBranchAddress("t1",&end.t);
            chain->SetBranchAddress("e1",&e1);
        }

        else
        {
            chain->SetBranchAddress("x0",&start.point.x);
            chain->SetBranchAddress("y0",&start.point.y);
            chain->SetBranchAddress("z0",&start.point.z);
            chain->SetBranchAddress("t0",&start.t);
            chain->SetBranchAddress("e0",&e0);
            chain->SetBranchAddress("x1",&end.point.x);
            chain->SetBranchAddress("y1",&end.point.y);
            chain->SetBranchAddress("z1",&end.point.z);
            chain->SetBranchAddress("t1",&end.t);
            chain->SetBranchAddress("e1",&e1);
        }
    }

    void MicroPoint::SetTTreeBranches(TTree* tree)
    {
        tree->SetBranchAddress("start.point.x",&start.point.x);
        tree->SetBranchAddress("start.point.y",&start.point.y);
        tree->SetBranchAddress("start.point.z",&start.point.z);
        tree->SetBranchAddress("start.t",&start.t);
        tree->SetBranchAddress("end.point.x",&end.point.x);
        tree->SetBranchAddress("end.point.y",&end.point.y);
        tree->SetBranchAddress("end.point.z",&end.point.z);
        tree->SetBranchAddress("end.t",&end.t);
        tree->SetBranchAddress("e0",&e0);
        tree->SetBranchAddress("e1",&e1);
    }





    //// Public MapPoint methods.

    void MapPoint::Diagonalize(bool no_z)
    {
        if (eigen_vals_set) return;

        TMatrixDSym cov(4);
        cov.SetMatrixArray(cov_mat.elements);
        if (no_z)
        {
            for (int i : {0,1,3})
            {
                cov(2,i) = 0;
                cov(i,2) = 0;
            }
        }
        
        if(cov.IsSymmetric() == false) throw std::runtime_error("Covariance matrix is not symmetric.");
        TMatrixDSymEigen eig(cov);
        eigen_vals = eig.GetEigenValues();
        eigen_vecs = eig.GetEigenVectors();

        eigen_vals_set = true;
    }

    EndPoint MapPoint::GetRandomPoint(TRandom3* rand)
    {
        EndPoint point = this->point;
        this->Diagonalize();

        for (int i = 0; i < 4; i++)
        {
            X17::EndPoint err_vec(eigen_vecs(0,i),eigen_vecs(1,i),eigen_vecs(2,i),eigen_vecs(3,i));
            
            // Biased (no 1/c_4) and maximum likelihood (1/N) stdev should be used for the random generation;
            double err = std::sqrt((this->n-1.0)/this->n * eigen_vals(i));
            // double err = std::sqrt(eigen_vals(i));

            point += rand->Gaus(0,err) * err_vec;
        }
        return point;
    }

    double MapPoint::operator[](int i) const
    {
        switch (i)
        {
        case 0:
            return this->point.point.x;
            break;
        case 1:
            return this->point.point.y;
            break;
        case 2:
            return this->point.point.z;
            break;
        case 3:
            return this->point.t;
            break;
        
        default:
            throw std::out_of_range("MapPoint[] index out of range.\n");
            break;
        }
    }





    //// Functions related to classes from Points.h.

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

    std::vector<TMarker3DBox*> GetDataMarkers(const std::vector<RecoPoint> &data, double zbin_size)
    {
        std::vector<TMarker3DBox*> markers;
        constexpr double max_size = 0.75;

        // Find maximal count.
        int max_count = 0;
        for (RecoPoint p : data) if (p.count > max_count) max_count = p.count;

        // Create markers.
        for (RecoPoint p : data)
        {
            using namespace constants;

            double rel_size = max_size * p.count / max_count;
            double xlen = rel_size * pad_width  / 2.0;
            double ylen = rel_size * pad_height / 2.0;
            double zlen = rel_size * zbin_size  / 2.0;

            markers.push_back(new TMarker3DBox(p.x(),p.y(),p.z(),xlen,ylen,zlen,0,0));
        }
        
        return markers;
    }
} // namespace X17