#pragma once

// ROOT dependencies
#include "Rtypes.h"

namespace X17
{
    /// @brief A struct for storing the coordinates and the time of the final point of an ionization electron using pads.
    struct EndPointDiscrete
    {
        int n_pad;    // The channel (pad) number.
        int time_bin; // Number of the time bin (time bin size 100 ns).

        ClassDefNV(EndPointDiscrete, 1);
    };
} // namespace X17