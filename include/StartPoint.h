#pragma once

// ROOT dependencies
#include "Rtypes.h"

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

        ClassDefNV(StartPoint, 1);
    };
} // namespace X17