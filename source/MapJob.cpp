// C++ dependencies
#include <cmath>
#include <iostream>
#include <string>

// X17 dependencies
#include "MapJob.h"

namespace X17
{
    void MapJob::GetParameters(int argc, char* argv[])
    {
        // Check the number of paramaters passed to main function.
        if (argc < 4)
        {
            std::cerr << "ERROR: Missing arguments in ion_electrons. Correct arguments: max_id, id, step (cm).\n";
        }

        if (argc < 5) iterations = 1;
        else          iterations = std::stoi(argv[4]);

        // Set parameters of this job.
        max_id = std::stoi(argv[1]);
        id     = std::stoi(argv[2]);
        step   = std::stod(argv[3]);

        if (id <= 0)         std::cerr << "ERROR: Parameter id has to be positive.\n";
        if (id > max_id)     std::cerr << "ERROR: Parameter id cannot be bigger than max_id.\n";
        if (step <= 0.0)     std::cerr << "ERROR: Parameter step has to be positive.\n";
        if (iterations <= 0) std::cerr << "ERROR: Parameter iterations has to be positive.\n";

        std::cout << "Running with parameters:\n" << "max_id: " << max_id << ", id: " << id << ", step (cm): " << step << ", iterations: " << iterations << ".\n";
    }

    void MapJob::SetElectronBounds()
    {
        f_GetStepParameters();
        n_el = nz * nxy;
        zsum = (zmax - zmin) * nxy * nz - nxy * step * nz * (nz - 1)/2;

        std::cout << "Expecting " << n_el << " electrons across all simulations, total propagation distance (per iteration) " << zsum << ".\n";

        min_el = f_Noptimal(id-1)+1;
        max_el = f_Noptimal(id);

        std::cout << "Simulating electrons from " << min_el << " to " << max_el << ", total propagation (per iteration) " << f_Zsum(max_el)-f_Zsum(min_el-1) << ".\n";
    }

    void MapJob::f_GetStepParameters()
    {
        nx  = floor((xmax-xmin)/step)+1; // Number of x steps in each iteration.
        nz  = floor((zmax-zmin)/step)+1; // Number of z steps in each iteration.
        nxy = 0;                         // Total number of x and y steps for each z step (to be calculated below).

        // Iterate over x and add up y steps to get nxy.
        for (int i = 0; i < nx; i++) 
        {
            double x     = xmin + i*step;                   // Current x value.
            int nymin    = floor((-ymin-x*sqrt(3))/step)+1; // The minimal number of y steps to satisfy minimal sector condition.
            double ymin2 = ymin + nymin*step;               // The smallest y to satisfy minimal sector condition.

            double ny; // Number of y steps for current x iteration.

            // If maximal sector condition cannot be satisfied set ny to 0, otherwise calculate it.
            if (ymin2 > (x*sqrt(3))) ny = 0;
            else ny = floor((x*sqrt(3) - ymin2)/step) + 1;
            nxy += ny;
        }
    }

    double MapJob::f_Zsum(int n_electrons)
    {
        int full_layers = floor(n_electrons / nxy); // Full z layers.
        int remainder   = n_electrons % nxy;        // The number of electrons in the incomplete z layer on top.

        return (zmax - zmin)*nxy*full_layers - nxy*step*full_layers*(full_layers - 1)/2 + remainder*(zmax - zmin - full_layers*step);
    }

    int MapJob::f_Noptimal(int l_id)
    {
        if (l_id == 0) return 0; // Id 0 is not used but best minimal index of id 1 depends on this.

        int min_el = 0;
        int max_el = nz*nxy;
        double opt_zsum = zsum/max_id; // Optimal total propagation for one job.

        int current_el = max_el/2;
        while (min_el != current_el)
        {
            if(f_Zsum(current_el) > l_id*opt_zsum) max_el = current_el;
            else min_el = current_el;
            current_el  = (max_el + min_el)/2;
        }  

        return max_el;
    }
} // namespace X17