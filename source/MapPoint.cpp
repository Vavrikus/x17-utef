// ROOT dependencies
#include "TMath.h"
#include "TMatrixDSymEigen.h"

// X17 dependencies
#include "MapPoint.h"

ClassImp(X17::MapPoint)

namespace X17
{
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
        
        if(cov.TMatrixDBase::IsSymmetric() == false) throw std::runtime_error("Covariance matrix is not symmetric.");
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

    TEllipse* MapPoint::GetErrorEllipse(int skip_index, double sigma, bool x_bigger)
    {
        double confidence = TMath::Erf(sigma/std::sqrt(2));
        double factor = TMath::ChisquareQuantile(confidence,2);

        if (skip_index < 0 || skip_index > 3)
            throw std::out_of_range("MapPoint::GetErrorEllipse skip_index out of range [0,3].");
        if (skip_index == 2)
            throw std::runtime_error("MapPoint::GetErrorEllipse z-coordinate already excluded.");

        int min_index = skip_index == 0 ? 1 : 0;
        int max_index = skip_index == 3 ? 1 : 3;
        if (x_bigger) std::swap(min_index,max_index);

        double a = this->cov_mat.at(min_index,min_index);
        double b = this->cov_mat.at(min_index,max_index);
        double d = this->cov_mat.at(max_index,max_index);

        double lambda_1 = (a + d + std::sqrt((a-d)*(a-d) + 4*b*b)) / 2;
        double lambda_2 = (a + d - std::sqrt((a-d)*(a-d) + 4*b*b)) / 2;

        X17::Vector v1(b, lambda_1 - a, 0);
        v1.Normalize();

        double theta = sign(v1.y)*v1.Angle({1,0,0})*180/M_PI;

        return new TEllipse((*this)[min_index],(*this)[max_index],std::sqrt(lambda_1*factor),std::sqrt(lambda_2*factor),0,360,theta);
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
            throw std::out_of_range("MapPoint[] index out of range [0,3].\n");
            break;
        }
    }
}