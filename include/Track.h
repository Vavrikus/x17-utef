#pragma once

// C++ dependencies
#include <vector>

// ROOT dependencies
#include "Rtypes.h"
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
        
        ClassDefNV(TrackRK, 1)
    };

    /// @brief Struct for information about microscopic simulated track.
    struct TrackMicro
    {
        bool electron;                  // True if electron, false if positron.
        std::vector<MicroPoint> points; // Simulated points.
        Vector origin;                  // Starting point of the track [cm].
        Vector orientation;             // Normalized initial direction of the track.
        double kin_energy;              // Kinetic energy of the particle [eV].

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
        
        ClassDefNV(TrackMicro, 1)
    };

    /// @brief Generates random initial parameters for track simulation.
    /// @param rand ROOT random number generator.
    /// @param electron Boolean, true if the particle is electron, false if it is positron.
    /// @param origin Starting point of the track.
    /// @param orient Normalized initial direction of the track.
    /// @param e_kin Kinetic energy of the particle.
    void GetRandomTrackParams(TRandom3* rand, bool& electron, Vector& origin, Vector& orient, double& e_kin);

    /// @brief Struct for simulated and reconstructed information about a track.
    struct TrackInfo
    {
        bool electron;     // True if electron, false if positron.
        double theta;      // Angle theta of entry into the TPC [deg].
        double phi;        // Angle phi of entry into the TPC [deg].
        double kin_energy; // Simulated kinetic energy of the particle [MeV].

        double cfit_nopads_energy_mid; // Circle fit energy (not accounting for pads) [MeV] calculated using middle field.
        double cfit_nopads_energy_avg; // Circle fit energy (not accounting for pads) [MeV] calculated using average field.
        double cfit_pads_energy_mid;   // Circle fit energy (accounting for pads) [MeV] calculated using middle field.
        double cfit_pads_energy_avg;   // Circle fit energy (accounting for pads) [MeV] calculated using average field.
        double rkfit_energy;           // Runge-Kutta fit energy (accounting for pads) [MeV].
        double rkfit_energy_err;       // Runge-Kutta fit energy error (accounting for pads) [MeV].

        /// @brief Function for calculating the relative error of the Runge-Kutta energy reconstruction.
        /// @return The relative error ΔE/E in % of the Runge-Kutta energy reconstruction.
        double GetRKRelError() { return 100 * (rkfit_energy - kin_energy) / kin_energy; }

        ClassDefNV(TrackInfo, 1)
    };
    
} // namespace X17