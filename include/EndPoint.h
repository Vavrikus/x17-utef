#pragma once

// C++ dependencies
#include <cassert>

// ROOT dependencies
#include "Rtypes.h"
#include "TGraph.h"

// X17 dependencies
#include "Matrix.h"
#include "Vector.h"

namespace X17
{    
    /// @brief A struct for storing the coordinates and the time of the final point of an ionization electron.
    struct EndPoint
    {
        Vector point; // Final coordinates [cm].
        double t;     // Final time [ns].

        /// @brief Default constructor that initializes the coordinates to 0 and time to -1.
        EndPoint() : point(), t(-1) { }

        /// @brief Constructor that takes individual double arguments for the coordinates and time.
        /// @param x The x-coordinate [cm].
        /// @param y The y-coordinate [cm].
        /// @param z The z-coordinate [cm].
        /// @param t The final time [ns].
        EndPoint(double x, double y, double z, double t) : point(x, y, z), t(t) { }
        
        /// @brief Constructor that takes a Vector object for the coordinates and a double for the time.
        /// @param point The final coordinates [cm].
        /// @param t The final time [ns].
        EndPoint(Vector point, double t) : point(point), t(t) { }

        /// @brief Conversion constructor for static casting zero to EndPoint.
        EndPoint(int zero) : point(Vector(0, 0, 0)), t(0) { assert(zero == 0); }

        /// @brief Default copy constructor.
        /// @param p The endpoint to be copied.
        EndPoint(const EndPoint& p) = default;

        /// @brief Getter for the x variable.
        /// @return The x-coordinate [cm].
        double x() const { return point.x; }

        /// @brief Getter for the y variable.
        /// @return The y-coordinate [cm].
        double y() const { return point.y; }

        /// @brief Getter for the z variable.
        /// @return The z-coordinate [cm].
        double z() const { return point.z; }

        /// @brief Squares the Endpoint object.
        /// @return EndPoint object with all coordinates squared.
        EndPoint Square() const { return EndPoint(point.x*point.x,point.y*point.y,point.z*point.z,t*t); }

        /// @brief Calculates the square root of the EndPoint object.
        /// @return EndPoint object with all coordinates equal to the square root of the original coordinates.
        EndPoint SquareRoot() const { return EndPoint(std::sqrt(point.x),std::sqrt(point.y),std::sqrt(point.z),std::sqrt(t)); }

        Matrix<4,1> ToVector4() const { return Matrix<4,1>({point.x,point.y,point.z,t}); }

        /// @brief Assignment operator.
        /// @param p The endpoint to be assigned.
        void operator=(const EndPoint& other)
        {
            point = other.point;
            t = other.t;
        }

        /// @brief Add another endpoint to this endpoint.
        /// @param[in] p The endpoint to add.
        void operator+=(EndPoint p)
        {
            this->point += p.point;
            this->t += p.t;
        }

        /// @brief Subtract another endpoint to this endpoint.
        /// @param[in] p The endpoint to subtract.
        void operator-=(EndPoint p)
        {
            this->point -= p.point;
            this->t -= p.t;
        }

        /// @brief Division by a scalar component-wise.
        /// @param d The scalar to divide with.
        /// @return The scaled endpoint.
        EndPoint operator/(double d) { return EndPoint(point / d, t / d); }

        /// @brief Division by a scalar component-wise.
        /// @param d The scalar to divide with.
        void operator/=(double d)
        {
            point /= d;
            t /= d;
        }

        /// @brief Adds two endpoints component-wise and returns the result.
        /// @param p2 The second endpoint.
        /// @return The sum of the two endpoints.
        EndPoint operator+(EndPoint p2) const
        {
            return EndPoint{point + p2.point, t + p2.t};
        }

        /// @brief Subtracts two endpoints component-wise.
        /// @param p2 The second endpoint.
        /// @return EndPoint The difference between the two endpoints.
        EndPoint operator-(EndPoint p2) const
        {
            return EndPoint{point - p2.point, t - p2.t};
        }

        ClassDefNV(EndPoint, 1);
    };

    /// @brief Scalar multiplication of a endpoint.
    /// @param d The scalar to multiply with.
    /// @param p The endpoint to multiply.
    /// @return The scaled endpoint.
    inline EndPoint operator*(double d, EndPoint p)
    {
        return EndPoint{d*p.point, d*p.t};
    }

    /// @brief Calculates the unbiased sample covariance matrix of a set of points.
    /// @param values The set of points.
    /// @param average The average of the points.
    /// @return The covariance matrix of the points.
    Matrix<4,4> CovarianceMatrix(const std::vector<EndPoint>& values, EndPoint average);

    /// @brief Calculates the squared Mahalanobis distances of a set of points.
    /// @param values The set of points.
    /// @param no_z If true, the z-coordinate is ignored.
    /// @return The squared Mahalanobis distances of the points.
    std::vector<double> SqMahalanobis(const std::vector<EndPoint>& values, bool no_z = true);

    /// @brief Fills a TGraph with the squared Mahalanobis distances of a set of points vs chi2 quantiles.
    /// @param graph TGraph to be filled.
    /// @param values The set of points.
    /// @param no_z If true, the z-coordinate is ignored.
    void FillQQplot(TGraph* graph, const std::vector<EndPoint>& values, bool no_z = true);

    /// @brief Performs a Mardia's test on a set of points.
    /// @param values The set of points.
    /// @param outA The output skewness.
    /// @param outB The output kurtosis.
    /// @param no_z If true, the z-coordinate is ignored.
    void MardiaTest(const std::vector<EndPoint>& values, double& outA, double& outB, bool no_z = true);
}