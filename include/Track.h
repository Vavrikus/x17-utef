#pragma once

#include "CircleFit3D.h"
#include "../VectorField.h"

/// @brief Struct for information about Runge-Kutta simulated track
struct TrackRK
{
    bool electron;                 // true if electron, false if positron
    std::vector<DataPoint> points; // simulated points
    Vector origin;                 // starting point of the track
    Vector orientation;            // normalized original direction of the track
    double kin_energy;             // kinetic energy of the track

    TrackRK() = default;

    TrackRK(const bool& e, const std::vector<DataPoint>& pts, const Vector& orig, const Vector& orient, const double& Ek)
        : electron(e), points(pts), origin(orig), orientation(orient), kin_energy(Ek) {}
};