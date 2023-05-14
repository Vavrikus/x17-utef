#pragma once

namespace X17
{
    /// @brief Struct representing a simulation job for simulating ionization electrons map.
    struct MapJob
    {
    public:
        // Most important parameters.
        int id;         // Id of the current job.
        int max_id;     // Maximal id of jobs for this simulation.
        int iterations; // Number of iterations to perform for each position.
        double step;    // Spacing of the grid [cm].

        // Ranges for simulation.
        // 1st sector condition (y-sqrt(3)*x<=0)&&(y+sqrt(3)*x>0) has to be satisfied later.
        double xmin =  0;  // Minimum x coordinate of simulation area.
        double xmax =  15; // Maximum x coordinate of simulation area.
        double ymin = -30; // Minimum y coordinate of simulation area.
        double ymax =  30; // Maximum y coordinate of simulation area.
        double zmin = -8;  // Minimum z coordinate of simulation area.
        double zmax =  8;  // Maximum z coordinate of simulation area.

        int nx;      // Number of steps in x direction.
        int nz;      // Number of steps in z direction.
        int nxy;     // Number of steps in x and y directions combined for each z step.
        int n_el;    // Total number of electrons in all jobs of this simulation.
        int min_el;  // Minimum index of electron to simulate in this job (inclusive).
        int max_el;  // Maximum index of electrons to simulate in this job (inclusive).

        double zsum; // Total distance to propagate all electrons (per iteration).

    public:
        /// @brief Parses command line arguments and sets job parameters.
        /// @param argc Number of command line arguments.
        /// @param argv Array of command line arguments.
        void SetParameters(int argc, char* argv[]);

        /// @brief Calculates number of electrons and total distance to propagate for current job. Sets minimal and maximal electron index.
        void SetElectronBounds();

    private:
        /// @brief Calculates number of steps in x and y directions for each z step and stores result in member variables.
        void _SetStepParameters();

        /// @brief Calculates total distance to propagate for given number of electrons.
        /// @param n_electrons Number of electrons to simulate.
        /// @return Total distance to propagate electrons.
        double _Zsum(int n_electrons) const;

        /// @brief Finds optimal maximal index of electron for given job's id.
        /// @param l_id Job id.
        /// @return Optimal maximal index of electron.
        int _Noptimal(int l_id) const;
    };
} // namespace X17