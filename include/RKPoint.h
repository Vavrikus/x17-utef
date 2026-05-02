#pragma once
// IWYU pragma: private, include "Points.h"

// ROOT dependencies
#include "Rtypes.h"

// X17 dependencies
#include "StartPoint.h"

namespace X17
{
  /// @brief A struct for storing a Runge-Kutta generated track point.
  struct RKPoint
  {
    StartPoint point; // Track point coordinates ([cm] and [ns]).

    /// @brief Default constructor.
    RKPoint()
      : point()
    {
    }

    /// @brief Constructor that takes individual double arguments for the coordinates and time.
    /// @param x The x-coordinate [cm].
    /// @param y The y-coordinate [cm].
    /// @param z The z-coordinate [cm].
    /// @param t The time [ns].
    RKPoint(double x, double y, double z, double t)
      : point(x, y, z, t)
    {
    }

    /// @brief Getter for the x variable.
    /// @return The x-coordinate [cm].
    [[nodiscard]] double x() const { return point.x(); }

    /// @brief Getter for the y variable.
    /// @return The y-coordinate [cm].
    [[nodiscard]] double y() const { return point.y(); }

    /// @brief Getter for the z variable.
    /// @return The z-coordinate [cm].
    [[nodiscard]] double z() const { return point.z(); }

    /// @brief Assignment operator.
    /// @param other The RKPoint to assign.
    RKPoint& operator=(const RKPoint& other) = default;

    /// @brief Returns the coordinates of the point as Vector object.
    /// @return The vector with coordinates of the point.
    [[nodiscard]] Vector AsVector() const { return point.point; }

    ClassDefNV(RKPoint, 1);
  };
} // namespace X17