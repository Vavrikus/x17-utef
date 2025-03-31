#pragma once

// X17 dependencies
#include "Field.h"
#include "Vector.h"

namespace X17
{
    namespace constants
    {
        // Physical constants
            constexpr double e    = 1.602176634E-19;  // Elementary charge
            constexpr double c    = 299792458;        // Light velocity in vacuum [m/s]
            constexpr double E0   = 510998.95;        // Electron rest energy [eV]
            constexpr double m0   = e * E0 / (c * c); // Electron rest mass [kg]
            constexpr double m2cm = 100.0;            // Conversion m --> cm
            constexpr double cm2m = 0.01;             // Conversion cm --> m

        // Target parameters
            constexpr double target_radius = 0.1; // Radius of the target [cm]

        // Vacuum tube parameters
            constexpr double tube_radius = 3.5; // Outer radius of the vacuum tube [cm]

        // Magnet parameters
            constexpr double mag_width  = 10.0; // Magnet width [cm]
            constexpr double mag_height = 15.0; // Magnet height [cm]
            constexpr double mag_depth  =  1.5; // Magnet depth [cm]

            constexpr double mag_lowxx = 4.626; // Magnet lowest-x corner x-coordinate [cm]
            constexpr double mag_lowxy = 3.537; // Magnet lowest-x corner y-coordinate [cm]

        // TPC window parameters
            constexpr double win_width  = 3.8; // TPC window width [cm]
            constexpr double win_height = 4.0; // TPC window height [cm]

        // TPC first sector (contains positive x axis) boundaries
            constexpr double zmin  =  -8.00; // TPC minimal height [cm]
            constexpr double zmax  =   8.00; // TPC maximal height [cm]
            constexpr double xmin  =   6.51; // TPC minimal x [cm]
            constexpr double xmax  =  14.61; // TPC maximal x [cm]
            constexpr double ylow  =   2.75; // TPC y of positive corner closer to center [cm]
            constexpr double yhigh =   7.45; // TPC y of positive corner further from center [cm]

        // TPC first sector derived parameters
            constexpr double yxslope    = (yhigh - ylow) / (xmax - xmin); // TPC Δy/Δx slanted side (y positive)
            constexpr double yintersect = ylow - yxslope * xmin;          // TPC slanted side (y positive) y-axis intersection (extrapolated) [cm]

        // TPC assumed electric field [V/cm]
            const Vector efield = Vector(0,0,-400);

        // TPC readout pads (counting rows and columns from top right corner)
            constexpr double pad_width   =  0.60;   // Pad width (horizontal) [cm]
            constexpr double pad_height  =  0.90;   // Pad height (vertical)  [cm]
            constexpr double pad_height2 =  0.60;   // Pad height [cm] last pad in 4th and 9th column (9th and 13th diagonal row)
            constexpr double pad_height3 =  0.509;  // Pad height [cm] last pad in 3rd column (14th diagonal row)
            constexpr double pad1_offset = -0.10;   // Vertical offset [cm] of corner of first pad from top right corner of trapezoid
            constexpr double pad_offset  =  0.08;   // Shortest distance between neighboring pads [cm]
            constexpr double pad_stag    =  0.3946; // Vertical offset [cm] of topmost pads of neighboring columns

            constexpr int rows     = 15;  // Number of pad rows (diagonal)
            constexpr int columns  = 12;  // Number of pad columns
            constexpr int channels = 128; // Number of pads (channels)
    } // namespace constants
    
    /// @brief Returns true if given point is in the first sector (containing positive x-axis) of the detector (trapezoidal prism).
    /// @param x The x-coordinate [cm].
    /// @param y The y-coordinate [cm].
    /// @param z The z-coordinate [cm].
    /// @param dist Specifies the minimal distance from within the TPC walls except the bottom anode wall (makes volume smaller).
    /// @return True if in sector, false otherwise.
    bool IsInSector(double x, double y, double z, double dist = 0);

    /// @brief True if vector is in the first sector (containing positive x-axis) of the detector (trapezoidal prism).
    /// @param vec The vector representing a point in space (x,y,z) [cm].
    /// @param dist Specifies the minimal distance from within the TPC walls except the bottom anode wall (makes volume smaller).
    /// @return True if in sector, false otherwise.
    inline bool IsInSector(Vector vec, double dist = 0)
    {
        return IsInSector(vec.x,vec.y,vec.z,dist);
    }
    
    /// @brief Returns minimal and maximal (magnetic) field inside first sector TPC volume.
    /// @param magfield Magnetic field data.
    /// @param out_min Output for minimal field [T].
    /// @param out_max Output for maximal field [T].
    /// @param dist Specifies the minimal distance from within the TPC walls except the bottom anode wall (makes volume smaller).
    void GetMinMaxField(const Field<Vector>& magfield, double& out_min, double& out_max, double dist);

    /// @brief Returns minimal (magnetic) field angle to electric field (along positive y-axis) inside first sector TPC volume.
    /// @param magfield Magnetic field data.
    /// @param out_min Output for minimal angle [rad].
    /// @param out_max Output for maximal angle [rad].
    /// @param dist Specifies the minimal distance from within the TPC walls except the bottom anode wall (makes volume smaller).
    void GetMinMaxFieldAngle(const Field<Vector>& magfield, double& out_min, double& out_max, double dist);

    /// @brief Draws lines around the approximate sensitive area (walls of TPC).
    /// @param yxformat If true, the y-coordinate is drawn on the x-axis and vice versa.
    void DrawTrapezoid(bool yxformat = true);

    /// @brief Draws the tube circle.
    void DrawTube();

    /// @brief Draws the two magnets extending to the first sector.
    /// @param yxformat If true, the y-coordinate is drawn on the x-axis and vice versa.
    void DrawFirstSectorMagnets(bool yxformat = true);
} // namespace X17
