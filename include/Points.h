#pragma once

// C++ dependencies
#include <cmath>

// ROOT dependencies
#include "TChain.h"
#include "TEllipse.h"
#include "TGraph.h"
#include "TMarker3DBox.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TTree.h"

// X17 dependencies
#include "Matrix.h"
#include "Utilities.h"
#include "Vector.h"

namespace X17
{
    /// @brief A struct for storing the coordinates and the time of the initial point of an ionization electron.
    struct StartPoint
    {
        Vector point; // Initial coordinates [cm].
        double t;     // Initial time [ns] (should be close to 0).

        /// @brief Default constructor that initializes the coordinates to 0 and time to -1.
        StartPoint() : point(), t(-1) { }

        /// @brief Constructor that takes individual double arguments for the coordinates and time.
        /// @param x The x-coordinate [cm].
        /// @param y The y-coordinate [cm].
        /// @param z The z-coordinate [cm].
        /// @param t The starting time [ns].
        StartPoint(double x, double y, double z, double t) : point(x, y, z), t(t) { }

        /// @brief Getter for the x variable.
        /// @return The x-coordinate [cm].
        double x() const { return point.x; }

        /// @brief Getter for the y variable.
        /// @return The y-coordinate [cm].
        double y() const { return point.y; }

        /// @brief Getter for the z variable.
        /// @return The z-coordinate [cm].
        double z() const { return point.z; }

        /// @brief Assignment operator.
        /// @param other The object to be assigned to this object.
        void operator=(const StartPoint& other)
        {
            point = other.point;
            t = other.t;
        }
    };

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
        EndPoint operator+(EndPoint p2)
        {
            return EndPoint{point + p2.point, t + p2.t};
        }

        /// @brief Subtracts two endpoints component-wise.
        /// @param p2 The second endpoint.
        /// @return EndPoint The difference between the two endpoints.
        EndPoint operator-(EndPoint p2)
        {
            return EndPoint{point - p2.point, t - p2.t};
        }
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

    /// @brief A struct for storing the coordinates and the time of the final point of an ionization electron using pads.
    struct EndPointDiscrete
    {
        int n_pad;    // The channel (pad) number.
        int time_bin; // Number of the time bin (time bin size 100 ns).
    };

    /// @brief A struct for storing the results of microscopic simulation of ionization electrons.
    struct MicroPoint
    {
        StartPoint start; // Initial coordinates of the electron ([cm] and [ns]).
        EndPoint end;     // Final coordinates of the electron ([cm] and [ns]).
        double e0;        // Initial energy [eV].
        double e1;        // Final energy [eV].

        /// @brief The default constructor. Initializes time to -1, everything else to 0.
        MicroPoint() : start(), end(), e0(0), e1(0) { }

        /// @brief Getter for the initial x variable.
        /// @return Initial x-coordinate [cm].
        double x0() const { return start.x(); }

        /// @brief Getter for the initial y variable.
        /// @return Initial y-coordinate [cm].
        double y0() const { return start.y(); }

        /// @brief Getter for the initial z variable.
        /// @return Initial z-coordinate [cm].
        double z0() const { return start.z(); }

        /// @brief Getter for the initial t variable.
        /// @return Initial time [ns].
        double t0() const { return start.t; }

        /// @brief Getter for the final x variable.
        /// @return Final x-coordinate [cm].
        double x1() const { return end.x(); }

        /// @brief Getter for the final y variable.
        /// @return Final y-coordinate [cm].
        double y1() const { return end.y(); }

        /// @brief Getter for the final z variable.
        /// @return Final z-coordinate [cm].
        double z1() const { return end.z(); }

        /// @brief Getter for the final t variable.
        /// @return Final time [ns].
        double t1() const { return end.t; }

        /// @brief Returns the initial position of the ionization electron.
        /// @return The initial position of the electron.
        Vector GetInitPos() { return this->start.point; }

        /// @brief Creates the branches of the given tree used for output.
        /// @param tree TTree used for output.
        void MakeTTreeBranches(TTree* tree);

        /// @brief Sets the branches of a TChain to this MicroPoint object.
        /// @param chain Pointer to TChain object.
        /// @param old_data Boolean indicating whether the data is from old simulations with different coordinate system (zxy). Defaults to true.
        void SetTChainBranches(TChain* chain, bool old_data = true);

        /// @brief Sets the branches of a TTree to this MicroPoint object. Uses different branch structure than MakeTTreeBranches.
        /// @param tree Pointer to a TTree object.
        void SetTTreeBranches(TTree* tree);
    };
    
    /// @brief A struct for storing a Runge-Kutta generated track point.
    struct RKPoint
    {
        StartPoint point; // Track point coordinates ([cm] and [ns]).

        /// @brief Default constructor.
        RKPoint() : point() { }

        /// @brief Constructor that takes individual double arguments for the coordinates and time.
        /// @param x The x-coordinate [cm].
        /// @param y The y-coordinate [cm].
        /// @param z The z-coordinate [cm].
        /// @param t The time [ns].
        RKPoint(double x, double y, double z, double t) : point(x,y,z,t) { }

        /// @brief Getter for the x variable.
        /// @return The x-coordinate [cm].
        double x() const { return point.x(); }

        /// @brief Getter for the y variable.
        /// @return The y-coordinate [cm].
        double y() const { return point.y(); }

        /// @brief Getter for the z variable.
        /// @return The z-coordinate [cm].
        double z() const { return point.z(); }

        /// @brief Assignment operator.
        /// @param other The RKPoint to assign.
        void operator=(const RKPoint& other) { this->point = other.point; }
        
        /// @brief Returns the coordinates of the point as Vector object.
        /// @return The vector with coordinates of the point.
        Vector AsVector() const { return point.point; }
    };

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
        MapPoint operator+(const MapPoint& p2)
        {
            if (n != p2.n)
                std::cerr << "WARNING: Adding map points with different number of electrons!\n";
            return MapPoint{n, point + p2.point, cov_mat + p2.cov_mat};
        }

        /// @brief Subtracts two map points component-wise.
        /// @param p2 The second map point.
        /// @return MapPoint The difference between the two map points.
        MapPoint operator-(const MapPoint& p2)
        {
            std::cout << "Helloooooooooooooooooooooooooooooooooooooooooooooooooooooo\n";
            if (n != p2.n)
                std::cerr << "WARNING: Subtracting map points with different number of electrons!\n";
            return MapPoint{n, point - p2.point, cov_mat - p2.cov_mat};
        }
    };

    /// @brief Scalar multiplication of a map point.
    /// @param d The scalar to multiply with.
    /// @param p The map point to multiply.
    /// @return The scaled map point.
    inline MapPoint operator*(double d, MapPoint p)
    {
        return MapPoint{p.n, d*p.point, d*p.cov_mat};
    }

    /// @brief A struct for storing the reconstructed initial points of the track. If points are reconstructed from pads, count/charge also relevant.
    struct RecoPoint
    {
        StartPoint point; // The reconstructed coordinates. [cm]
        double count;     // The number of electrons or charge.

        /// @brief Default constructor. Coordinates x,y,z set to zero, time to -1.
        RecoPoint() : point() { }

        /// @brief Constructor of RecoPoint. Sets time to -1.
        /// @param x The x-coordinate [cm].
        /// @param y The y-coordinate [cm].
        /// @param z The z-coordinate [cm].
        /// @param count The number of electrons or charge.
        RecoPoint(double x, double y, double z, int count) : point(x,y,z,-1), count(count) { }

        /// @brief Getter for the x variable.
        /// @return The x-coordinate [cm].
        double x() const { return point.x(); }

        /// @brief Getter for the y variable.
        /// @return The y-coordinate [cm].
        double y() const { return point.y(); }

        /// @brief Getter for the z variable.
        /// @return The z-coordinate [cm].
        double z() const { return point.z(); }

        /// @brief Assignment operator.
        /// @param other The object to be assigned to this object.
        void operator=(const RecoPoint& other)
        {
            point = other.point;
            count = other.count;
        }

        /// @brief Returns the coordinates of the point as Vector object.
        /// @return The vector with coordinates of the point.
        Vector AsVector() const { return point.point; }
    };

    std::vector<TMarker3DBox*> GetDataMarkers(const std::vector<RecoPoint>& data, double zbin_size = 0.3);

    /// @brief A struct for storing the simulated (later maybe also real) data from the TPC readout.
    struct DataPoint
    {
        EndPointDiscrete point; // The binned time and position information.
        double count;           // The number of electrons or charge.
    };

    /// @brief A struct for storing a point on an ionization electron driftline.
    struct DriftLinePoint
    {
        Vector point; // Current ionization electron position [cm].
        double t;     // Current time [ns].

        /// @brief Default constructor that initializes the coordinates to 0 and time to -1.
        DriftLinePoint() : point(), t(-1) { }

        /// @brief Constructor that takes individual double arguments for the coordinates and time.
        /// @param x The x-coordinate [cm].
        /// @param y The y-coordinate [cm].
        /// @param z The z-coordinate [cm].
        /// @param t The time [ns].
        DriftLinePoint(double x, double y, double z, double t) : point(x, y, z), t(t) { }

        /// @brief Getter for the x variable.
        /// @return The x-coordinate [cm].
        double x() const { return point.x; }

        /// @brief Getter for the y variable.
        /// @return The y-coordinate [cm].
        double y() const { return point.y; }

        /// @brief Getter for the z variable.
        /// @return The z-coordinate [cm].
        double z() const { return point.z; }
    };
} // namespace X17