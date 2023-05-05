#pragma once

#include <vector>

#include "Points.h"

namespace X17
{
    /// @brief Struct for information about Runge-Kutta simulated track
    struct TrackRK
    {
        bool electron;                 // True if electron, false if positron.
        std::vector<RKPoint> points;   // Simulated points.
        Vector origin;                 // starting point of the track
        Vector orientation;            // normalized original direction of the track
        double kin_energy;             // kinetic energy of the track

        TrackRK() = default;

        TrackRK(const bool& e, const std::vector<RKPoint>& pts, const Vector& orig, const Vector& orient, const double& Ek)
            : electron(e), points(pts), origin(orig), orientation(orient), kin_energy(Ek) {}
    };
} // namespace X17