#pragma once

// ROOT dependencies
#include "Rtypes.h"

// X17 dependencies
#include "Vector.h"

namespace X17
{
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

        ClassDefNV(DriftLinePoint, 1)
    };
} // namespace X17