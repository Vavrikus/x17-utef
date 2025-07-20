#pragma once

// ROOT dependencies
#include "Rtypes.h"

// X17 dependencies
#include "EndPointDiscrete.h"

namespace X17
{
    /// @brief A struct for storing the simulated (later maybe also real) data from the TPC readout.
    struct DataPoint
    {
        EndPointDiscrete point; // The binned time and position information.
        double count;           // The number of electrons or charge.

        ClassDefNV(DataPoint, 1)
    };
} // namespace X17