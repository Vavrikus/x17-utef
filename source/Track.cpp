// ROOT dependencies
#include "TRandom3.h"

// X17 dependencies
#include "Track.h"
#include "Utilities.h"
#include "X17Utilities.h"

namespace X17
{
    void GetRandomTrackParams(TRandom3* rand, bool& electron, Vector& origin, Vector& orient, double& e_kin)
    {
        using namespace X17::constants;

        // Assuming that the target is in the YZ plane.
        constexpr double x0     =  0;              // The x coordinate of simulated origin [cm].
        constexpr double r_max  =  target_radius;  // The maximal distance from origin [cm].

        constexpr double x1     =  xmin;           // The x coordinate of simulated window point [cm].
        constexpr double y1_min = -win_width / 2;  // The minimal y coordinate of simulated window point [cm].
        constexpr double y1_max = -y1_min;         // The maximal y coordinate of simulated window point [cm].
        constexpr double z1_min = -win_height / 2; // The minimal z coordinate of simulated window point [cm].
        constexpr double z1_max = -z1_min;         // The maximal z coordinate of simulated window point [cm].

        constexpr double e_min  =  3e+6;           // The minimal simulated energy [eV].
        constexpr double e_max  =  13e+6;          // The maximal simulated energy [eV].
        

        // Simulation of the initial track parameters.
        double phi = RandomMinMax(rand,0,2*M_PI);                 // The azimuth of the initial point on the circle target.
        double r   = std::sqrt(RandomMinMax(rand,0,r_max*r_max)); // The distance of the initial point from the circle target.
        double y0  = r*cos(phi);                                  // The y-coordinate of the initial point.
        double z0  = r*sin(phi);                                  // The z-coordinate of the initial point.

        double y1 = RandomMinMax(rand,y1_min,y1_max);             // The y-coordinate of the window point (TPC entry).
        double z1 = RandomMinMax(rand,z1_min,z1_max);             // The z-coordinate of the window point (TPC entry).

        origin = {x1, y1, z1};                                    // Setting the initial point (or origin).
        orient = {x1 - x0, y1 - y0, z1 - z0};                     // Setting the initial orientation.
        orient.Normalize();

        e_kin = RandomMinMax(rand,e_min,e_max);                   // The kinetic energy of the particle.

        electron = rand->Rndm() > 0.5;                            // Choosing either electron or positron.
    }
} // namespace X17
