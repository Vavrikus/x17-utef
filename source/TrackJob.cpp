// C++ dependencies
#include <cmath>
#include <iostream>

// X17 dependencies
#include "TrackJob.h"

namespace X17
{
    //// Public methods.    
    void TrackJob::SetParameters(int argc, char* argv[])
    {
        if (argc == 1)
        {
            random     = true;
            min_set    = 1;
            max_set    = 1;
            iterations = 1;
            std::cout << "\nRunning random generation of a single track, no extra parameters needed.\n\n";
        }
        else if (argc == 2)
        {
            random     = true;
            min_set    = 1;
            max_set    = 1;
            iterations = std::stoi(argv[1]);
            std::cout << "\nRunning random generation of " << iterations << " tracks.\n\n";
        }
        else
        {
            random = false;

            // Check the number of paramaters passed to main function.
            if (argc < 6)
            {
                std::cerr << "ERROR: Missing arguments in ion_electrons. Correct arguments: max_id, id, iterations, angle_bins, energy_bins.\n";
            }

            // Set parameters of this job.
            max_id      = std::stoi(argv[1]);
            id          = std::stoi(argv[2]);
            iterations  = std::stoi(argv[3]);
            angle_bins  = std::stoi(argv[4]);
            energy_bins = std::stoi(argv[5]);

            if (id <= 0)          std::cerr << "ERROR: Parameter id has to be positive.\n";
            if (id > max_id)      std::cerr << "ERROR: Parameter id cannot be bigger than max_id.\n";
            if (iterations <= 0)  std::cerr << "ERROR: Parameter iterations has to be positive.\n";
            if (angle_bins <= 1)  std::cerr << "ERROR: Parameter angle_bins has to be greater than 1.\n";
            if (energy_bins <= 1) std::cerr << "ERROR: Parameter energy_bins has to be greater than 1.\n";

            std::cout << "\nRunning with parameters:\n" << "max_id: " << max_id << ", id: " << id << ", iterations: " << iterations;
            std::cout << ", angle_bins: " << angle_bins << ", energy_bins: " << energy_bins << ".\n";

            // Integer division to avoid rounding error problems.
            n_sets  = 2 * angle_bins * angle_bins * energy_bins;
            min_set = ((n_sets * (id - 1)) / max_id) + 1;
            max_set = (n_sets * id ) / max_id;

            std::cout << "Number of unique tracks in all jobs: " << n_sets << ", in this job tracks from " << min_set << " to " << max_set;
            std::cout << " will be simulated.\n\n";
        }
    }

    void TrackJob::GetTrackParameters(int n_set, bool& electron, Vector& origin, Vector& orient, double& e_kin)
    {
        // Getting bin info from the value of n_set.
        electron     = (n_set - 1) % 2;
        int quotient = (n_set - 1) / 2;

        int phi_bin = quotient % angle_bins; // Bin number for phi between 0 and angle_bins - 1.
        quotient /= angle_bins;

        int theta_bin = quotient % angle_bins; // Bin number for theta between 0 and angle_bins - 1.
        quotient /= angle_bins;

        int E_bin = quotient; // Bin number for energy between 0 and energy_bins - 1.

        // Calculating phi, theta and energy from the bin info.
        double phi   = phi_min   + (phi_max   - phi_min)   * phi_bin   * (1.0 / (angle_bins  - 1));
        double theta = theta_min + (theta_max - theta_min) * theta_bin * (1.0 / (angle_bins  - 1));
        e_kin        = E_min     + (E_max     - E_min)     * E_bin     * (1.0 / (energy_bins - 1));

        // Setting the orientation from theta and phi.
        orient = {cos(phi)*cos(theta), sin(phi)*cos(theta), sin(theta)};

        // Setting the origin point from the orientation vector (assuming straight line motion from (0,0,0)).
        origin = orient * constants::xmin / (cos(phi)*cos(theta));
    }
} // namespace X17