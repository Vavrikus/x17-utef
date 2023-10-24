#pragma once

// C++ dependencies
#include <cmath>

// X17 dependencies
#include "X17Utilities.h"

namespace X17
{
    /// @brief Struct representing a job for track simulation.
    struct TrackJob
    {
        // Main parameters of the job.
        int id;         // ID of the current job.
        int max_id;     // Maximal ID of jobs in this simulation.
        int iterations; // Number of iterations for each set of track parameters.
        bool random;    // Should the tracks be generated with random parameters?

        int angle_bins;  // Number of different theta and phi values in simulation.
        int energy_bins; // Number of different energy values in simulation.

        // Ranges for simulation.
        double theta_max = atan((constants::win_height/2)/constants::xmin); // The maximal simulated theta [rad].
        double theta_min = -theta_max;                                      // The minimal simulated theta [rad].
        double phi_max = atan((constants::win_width/2)/constants::xmin);    // The maximal simulated phi [rad].
        double phi_min = -phi_max;                                          // The minimal simulated phi [rad].
        double E_max = 13e+6;                                               // The maximal simulated energy [eV].
        double E_min = 3e+6;                                                // THe minimal simulated energy [eV].

        int n_sets;  // Total number of unique track parameter sets in all jobs.
        int min_set; // Lowest identifying number of track parameter set in this job.
        int max_set; // Highest identifying number of track parameter set in this job.


        /// @brief Parses command line arguments and sets job parameters.
        /// @param argc Number of command line arguments.
        /// @param argv Array of command line arguments.
        void SetParameters(int argc, char* argv[]);

        /// @brief Calculates initial parameters for track simulation.
        /// @param n_set The number of the current set of track parameters.
        /// @param rand ROOT random number generator.
        /// @param electron Boolean, true if the particle is electron, false if it is positron.
        /// @param origin Starting point of the track.
        /// @param orient Normalized initial direction of the track.
        /// @param e_kin Kinetic energy of the particle.
        void GetTrackParameters(int n_set, bool& electron, Vector& origin, Vector& orient, double& e_kin);
    };
    
} // namespace X17