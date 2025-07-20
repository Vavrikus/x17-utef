#pragma once

// ROOT dependencies
#include "TMarker3DBox.h"

// X17 dependencies
#include "StartPoint.h"

namespace X17
{
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
}