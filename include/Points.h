#pragma once

// C++ dependencies
#include <cmath>

// ROOT dependencies
#include "TChain.h"
#include "TMarker3DBox.h"
#include "TTree.h"

// X17 dependencies
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
        EndPoint(int zero) : point(Vector(0, 0, 0)), t(0) { }

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
        EndPoint Square() { return EndPoint(point.x*point.x,point.y*point.y,point.z*point.z,t*t); }

        /// @brief Calculates the square root of the EndPoint object.
        /// @return EndPoint object with all coordinates equal to the square root of the original coordinates.
        EndPoint SquareRoot() { return EndPoint(sqrt(point.x),sqrt(point.y),sqrt(point.z),sqrt(t)); }

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
        EndPoint point; // The average final point of the ionization electrons ([cm] and [ns]).
        EndPoint dev;   // The standard deviations of individual points ([cm] and [ns]).

        /// @brief Default constructor. Coordinates x,y,z set to zero, time to -1.
        MapPoint() : point(), dev() { }

        /// @brief Constructor for initializing a MapPoint object with given values.
        /// @param p The EndPoint with (x,y,z,t) coordinates ([cm] and [ns]).
        /// @param d The EndPoint with (x,y,z,t) coordinates deviations ([cm] and [ns]).
        MapPoint(EndPoint p, EndPoint d) : point(p), dev(d) { }

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
        double xdev() const { return dev.x(); }

        /// @brief Getter for the y variable deviation.
        /// @return The y-coordinate deviation [cm].
        double ydev() const { return dev.y(); }

        /// @brief Getter for the z variable deviation.
        /// @return The z-coordinate deviation [cm].
        double zdev() const { return dev.z(); }

        /// @brief Getter for the t variable deviation.
        /// @return Time deviation [ns].
        double tdev() const { return dev.t; }

        /// @brief Assignment operator.
        /// @param p The MapPoint to be assigned to this MapPoint.
        void operator=(const MapPoint& p)
        {
            point = p.point;
            dev   = p.dev;
        }

        /// @brief Add another map point to this map point.
        /// @param[in] p The map point to add.
        void operator+=(MapPoint p)
        {
            this->point += p.point;
            this->dev += p.dev;
        }

        /// @brief Accesses the individual components of the endpoint.
        /// @param i Index of the component to access. Must be in the range [0, 3].
        /// @return The value of the component.
        /// @throws std::out_of_range if i is out of range.
        double operator[](int i) const;

        /// @brief Adds two map points component-wise and returns the result.
        /// @param p2 The second map point.
        /// @return The sum of the two map points.
        MapPoint operator+(MapPoint p2)
        {
            return MapPoint{point + p2.point, dev + p2.dev};
        }

        /// @brief Subtracts two map points component-wise.
        /// @param p2 The second map point.
        /// @return MapPoint The difference between the two map points.
        MapPoint operator-(MapPoint p2)
        {
            return MapPoint{point - p2.point, dev - p2.dev};
        }
    };

    /// @brief Scalar multiplication of a map point.
    /// @param d The scalar to multiply with.
    /// @param p The map point to multiply.
    /// @return The scaled map point.
    inline MapPoint operator*(double d, MapPoint p)
    {
        return MapPoint{d*p.point, d*p.dev};
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
} // namespace X17