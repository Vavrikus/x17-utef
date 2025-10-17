#pragma once

// ROOT dependencies
#include "Rtypes.h"
#include "TEllipse.h"
#include "TMatrixD.h"
#include "TRandom3.h"
#include "TVectorD.h"

// X17 dependencies
#include "EndPoint.h"
#include "Utilities.h"

namespace X17
{
    /// @brief A struct for storing the results of map simulation (multiple ionization electrons with same initial coordinates).
    struct MapPoint
    {
        int n;               // The number of ionization electrons in the simulation.
        EndPoint point;      // The average final point of the ionization electrons ([cm] and [ns]).
        Matrix<4,4> cov_mat; // The covariance matrix of individual points ([cm] and [ns]).

        bool eigen_vals_set = false; // True if eigenvalues and eigenvectors are set.
        TMatrixD eigen_vecs; // The eigenvectors of the covariance matrix.
        TVectorD eigen_vals; // The eigenvalues of the covariance matrix.

        double mardia_A = -1; // Mardia's skewness statistic p-value (-1 if not calculated).
        double mardia_B = -1; // Mardia's kurtosis statistic p-value (-1 if not calculated).

        /// @brief Default constructor. Coordinates x,y,z set to zero, time to -1.
        MapPoint() : n(0), point(), cov_mat(0), eigen_vecs(4,4), eigen_vals(4) { }

        /// @brief Constructor for initializing a MapPoint object with given values.
        /// @param n The number of ionization electrons in the simulation.
        /// @param p The EndPoint with (x,y,z,t) coordinates averages ([cm] and [ns]).
        /// @param m The 4x4 covariance matrix of (x,y,z,t) coordinates ([cm] and [ns]).
        MapPoint(int n, EndPoint p, Matrix<4,4> m) : n(n), point(p), cov_mat(m), eigen_vecs(4,4), eigen_vals(4) { }

        /// @brief Diagonalize the covariance matrix and store the eigenvalues and eigenvectors in eigen_values and eigen_vectors.
        /// @param no_z If true, the z-coordinate is not included in the diagonalization.
        void Diagonalize(bool no_z = true);

        /// @brief Returns a random point using the covariance matrix for multivariate normal distribution
        /// @return The random point.
        EndPoint GetRandomPoint(TRandom3* rand);

        /// @brief Caculates the error ellipse of the point for X, Y, T coordinates.
        /// @param skip_index The index of the coordinate to skip (must be in range [0,3] -- see operator[]).
        /// @param sigma Scales the error ellipse using the normal distribution quantiles.
        /// @param x_bigger Should the coordinate with bigger index be used for the x-axis?
        /// @return Pointer to the created TEllipse object.
        TEllipse* GetErrorEllipse(int skip_index, double sigma = 1.0, bool x_bigger = false);

        /// @brief Getter for the x variable.
        /// @return The x-coordinate [cm].
        double x() const { return point.x(); }

        /// @brief Getter for the y variable.
        /// @return The y-coordinate [cm].
        double y() const { return point.y(); }

        /// @brief Getter for the z variable.
        /// @return The z-coordinate [cm].
        double z() const { return point.z(); }

        /// @brief Getter for the t variable.
        /// @return Time [ns].
        double t() const { return point.t; }

        /// @brief Getter for the x variable deviation.
        /// @return The x-coordinate deviation [cm].
        double xdev() const { return StdevBiasFactor(n)*std::sqrt(cov_mat.at(0,0)); }

        /// @brief Getter for the y variable deviation.
        /// @return The y-coordinate deviation [cm].
        double ydev() const { return StdevBiasFactor(n)*std::sqrt(cov_mat.at(1,1)); }

        /// @brief Getter for the z variable deviation.
        /// @return The z-coordinate deviation [cm].
        double zdev() const { return StdevBiasFactor(n)*std::sqrt(cov_mat.at(2,2)); }

        /// @brief Getter for the t variable deviation.
        /// @return Time deviation [ns].
        double tdev() const { return StdevBiasFactor(n)*std::sqrt(cov_mat.at(3,3)); }

        /// @brief Assignment operator.
        /// @param p The MapPoint to be assigned to this MapPoint.
        void operator=(const MapPoint& p)
        {
            n = p.n;
            point = p.point;
            cov_mat   = p.cov_mat;
        }

        /// @brief Add another map point to this map point.
        /// @param[in] p The map point to add.
        void operator+=(const MapPoint& p)
        {
            if (n != p.n)
                std::cerr << "WARNING: Adding map points with different number of electrons!\n";
            this->point += p.point;
            this->cov_mat += p.cov_mat;
        }

        /// @brief Accesses the individual components of the endpoint.
        /// @param i Index of the component to access. Must be in the range [0, 3].
        /// @return The value of the component.
        /// @throws std::out_of_range if i is out of range.
        double operator[](int i) const;

        /// @brief Adds two map points component-wise and returns the result.
        /// @param p2 The second map point.
        /// @return The sum of the two map points.
        MapPoint operator+(const MapPoint& p2) const
        {
            if (n != p2.n)
                std::cerr << "WARNING: Adding map points with different number of electrons!\n";
            return MapPoint{n, point + p2.point, cov_mat + p2.cov_mat};
        }

        /// @brief Subtracts two map points component-wise.
        /// @param p2 The second map point.
        /// @return MapPoint The difference between the two map points.
        MapPoint operator-(const MapPoint& p2) const
        {
            if (n != p2.n)
                std::cerr << "WARNING: Subtracting map points with different number of electrons!\n";
            return MapPoint{n, point - p2.point, cov_mat - p2.cov_mat};
        }

        ClassDefNV(MapPoint, 1)
    };

    /// @brief Scalar multiplication of a map point.
    /// @param d The scalar to multiply with.
    /// @param p The map point to multiply.
    /// @return The scaled map point.
    inline MapPoint operator*(double d, MapPoint p)
    {
        return MapPoint{p.n, d*p.point, d*p.cov_mat};
    }
}