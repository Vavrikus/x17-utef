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

    void MapPoint::Diagonal(TVectorD &eigen_values, TMatrixD &eigen_vectors, bool no_z) const
    {
        TMatrixDSym cov(no_z ? 3 : 4);
        if(no_z)
        {
            cov = TMatrixDSym(3);
            const int indices[3] = {0,1,3};
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    cov(i,j) = cov_mat.at(indices[i],indices[j]);
        }
        else cov.SetMatrixArray(cov_mat.elements);
        
        if(cov.IsSymmetric() == false) throw std::runtime_error("Covariance matrix is not symmetric.");
        TMatrixDSymEigen eig(cov);
        eigen_values = eig.GetEigenValues();
        eigen_vectors = eig.GetEigenVectors();

        // TMatrixD test = eigen_vectors*TMatrixD(TMatrixD::kTransposed,eigen_vectors);
        // std::cout << "Covariance matrix:\n";
        // cov.Print();
        // std::cout << "Eigenvalues:\n";
        // eigen_values.Print();
        // std::cout << "Eigenvectors:\n";
        // eigen_vectors.Print();
        // std::cout << "Identity matrix?:\n";
        // test.Print();
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

    Matrix<4, 4> CovarianceMatrix(std::vector<EndPoint> points, EndPoint average)
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
        for (Vector4 v : vectors)
            for (int i = 0; i < 4; i++)
                for (int j = i; j < 4; j++)
                    cov_mat.at(i,j) += v.at(i,0) * v.at(j,0);

        cov_mat /= vectors.size() - 1;

        for(int i = 0; i < 4; i++)
            for (int j = i; j < 4; j++)
                cov_mat.at(j,i) = cov_mat.at(i,j);

        return cov_mat;
    }

    std::vector<double> SqMahalanobis(std::vector<EndPoint> values, EndPoint average, Matrix<4, 4> cov)
    {
        std::vector<double> sq_mahal;
        sq_mahal.reserve(values.size());

        TMatrixD cov_m(4,4,cov.elements);
        TMatrixD cov_inv(TMatrixD::kInverted,cov_m);

        for (EndPoint p : values)
        {
            Matrix<4,1> v = p.ToVector4() - average.ToVector4();
            TVectorD tv(4,v.elements);
            sq_mahal.push_back(tv * (cov_inv * tv));
        }

        return sq_mahal;
    }

    void FillQQplot(TGraph* graph, std::vector<EndPoint> values, EndPoint average)
    {
        std::vector<double> sq_mahal = SqMahalanobis(values,average,CovarianceMatrix(values,average));
        std::sort(sq_mahal.begin(),sq_mahal.end());
        
        for (int i = 0; i < sq_mahal.size(); i++)
        {
            double chi2_quantile = TMath::ChisquareQuantile((i+1.0)/(sq_mahal.size()+1),4);
            graph->AddPoint(chi2_quantile,sq_mahal[i]);
        }
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