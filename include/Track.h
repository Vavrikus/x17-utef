#pragma once

// C++ dependencies
#include <vector>

// ROOT dependencies
#include "TRandom3.h"

// X17 dependencies
#include "Points.h"
#include "Vector.h"

namespace X17
{
    /// @brief Struct for information about Runge-Kutta simulated track.
    struct TrackRK
    {
        bool electron;                 // True if electron, false if positron.
        std::vector<RKPoint> points;   // Simulated points.
        Vector origin;                 // Starting point of the track.
        Vector orientation;            // Normalized initial direction of the track.
        double kin_energy;             // Kinetic energy of the particle.

        /// @brief Default constructor.
        TrackRK() = default;

        /// @brief Constructor of TrackRK.
        /// @param e Boolean, true if the particle is electron, false if it is positron.
        /// @param pts Vector of simulated points.
        /// @param orig Starting point of the track.
        /// @param orient Normalized initial direction of the track.
        /// @param Ek Kinetic energy of the particle.
        TrackRK(bool e, const std::vector<RKPoint>& pts, Vector orig, Vector orient, double Ek)
            : electron(e), points(pts), origin(orig), orientation(orient), kin_energy(Ek) { }
    };

    /// @brief Struct for information about microscopic simulated track.
    struct TrackMicro
    {
        bool electron;                  // True if electron, false if positron.
        std::vector<MicroPoint> points; // Simulated points.
        Vector origin;                  // Starting point of the track.
        Vector orientation;             // Normalized initial direction of the track.
        double kin_energy;              // Kinetic energy of the particle.

        std::vector<std::vector<DriftLinePoint>> driftlines; // Simulated drift lines of each electron.

        /// @brief Default constructor.
        TrackMicro() = default;

        /// @brief Constructor of TrackRK.
        /// @param e Boolean, true if the particle is electron, false if it is positron.
        /// @param pts Vector of simulated points.
        /// @param orig Starting point of the track.
        /// @param orient Normalized initial direction of the track.
        /// @param Ek Kinetic energy of the particle.
        /// @param dlines Vector of drift lines.
        TrackMicro(bool e, const std::vector<MicroPoint>& pts, Vector orig, Vector orient, double Ek, const std::vector<std::vector<DriftLinePoint>>& dlines)
            : electron(e), points(pts), origin(orig), orientation(orient), kin_energy(Ek), driftlines(dlines) { }
    };

    /// @brief Generates random initial parameters for track simulation.
    /// @param rand ROOT random number generator.
    /// @param electron Boolean, true if the particle is electron, false if it is positron.
    /// @param origin Starting point of the track.
    /// @param orient Normalized initial direction of the track.
    /// @param e_kin Kinetic energy of the particle.
    void GetRandomTrackParams(TRandom3* rand, bool& electron, Vector& origin, Vector& orient, double& e_kin);
} // namespace X17