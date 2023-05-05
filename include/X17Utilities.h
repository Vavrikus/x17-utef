#pragma once

#include "TLine.h"

#include "Vector.h"
namespace X17
{
    template <typename T> class Field;
    void GetMinMaxField(const Field<Vector>& magfield, double& out_min, double& out_max, const double& dist);
    void GetMinMaxFieldAngle(const Field<Vector>& magfield, double& out_min, double& out_max, const double& dist);
}
#include "Field.h"

namespace X17
{
    namespace constants
    {
        // Physical constants
            constexpr double e    = 1.602176634E-19; // Elementary charge
            constexpr double c    = 299792458;       // Light velocity in vacuum [m/s]
            constexpr double E0   = 510998.95;       // Electron rest energy [eV]
            constexpr double m0   = e*E0/(c*c);      // Electron rest mass [kg]
            constexpr double m2cm = 100.0;           // Conversion m --> cm

        // Target parameters
            constexpr double target_radius = 0.1; // Radius of the target [cm]

        // TPC window parameters
            constexpr double win_width  = 3.8; // TPC window width [cm]
            constexpr double win_height = 4.0; // TPC window height [cm]

        // TPC first sector (contains positive x axis) boundaries
            constexpr double zmin  = -8;     // TPC minimal height [cm]
            constexpr double zmax  =  8;     // TPC maximal height [cm]
            constexpr double xmin  =  6.51;  // TPC minimal x [cm]
            constexpr double xmax  =  14.61; // TPC maximal x [cm]
            constexpr double ylow  =  2.75;  // TPC y of positive corner closer to center [cm]
            constexpr double yhigh =  7.45;  // TPC y of positive corner further from center [cm]

        // TPC first sector derived parameters
            constexpr double yxslope    = (yhigh - ylow) / (xmax - xmin); // TPC Δy/Δx slanted side (y positive)
            constexpr double yintersect = ylow - yxslope * xmin;          // TPC slanted side (y positive) y-axis intersection (extrapolated) [cm]

        // TPC assumed electric field [V/cm]
            const Vector efield = Vector(0,0,-400);

        // TPC readout pads (counting rows and columns from top right corner)
            constexpr double pad_width   =  0.6;   // Pad width (horizontal) [cm]
            constexpr double pad_height  =  0.9;   // Pad height (vertical)  [cm]
            constexpr double pad_height2 =  0.6;   // Pad height [cm] last pad in 4th and 9th column (9th and 13th diagonal row)
            constexpr double pad_height3 = 0.509;  // Pad height [cm] last pad in 3rd column (14th diagonal row)
            constexpr double pad1_offset = -0.1;   // Vertical offset [cm] of corner of first pad from top right corner of trapezoid
            constexpr double pad_offset  = 0.08;   // Shortest distance between neighboring pads [cm]
            constexpr double pad_stag    = 0.3946; // Vertical offset [cm] of topmost pads of neighboring columns

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
    bool IsInSector(const double& x, const double& y, const double& z, const double& dist = 0)
    {
        using namespace constants;

        if (z < zmin      || z > zmax-dist) return false;
        if (x < xmin+dist || x > xmax-dist) return false;

        double dy = dist / std::sqrt(1 + 1/yxslope); // change of allowed y (absolute value) for slanted surface
        if (std::abs(y) + dy > yxslope*x + yintersect) return false;
        
        return true;   
    }

    /// @brief True if vector is in the first sector (containing positive x-axis) of the detector (trapezoidal prism).
    /// @param vec The vector representing a point in space (x,y,z) [cm].
    /// @param dist Specifies the minimal distance from within the TPC walls except the bottom anode wall (makes volume smaller).
    /// @return True if in sector, false otherwise.
    bool IsInSector(const Vector& vec, const double& dist = 0)
    {
        return IsInSector(vec.x,vec.y,vec.z,dist);
    }
    
    /// @brief Returns minimal and maximal (magnetic) field inside first sector TPC volume.
    /// @param magfield Magnetic field data.
    /// @param out_min Output for minimal field [T].
    /// @param out_max Output for maximal field [T].
    /// @param dist Specifies the minimal distance from within the TPC walls except the bottom anode wall (makes volume smaller).
    void GetMinMaxField(const Field<Vector>& magfield, double& out_min, double& out_max, const double& dist)
    {
        double min_magnitude_sq = -1;
        double max_magnitude_sq = -1;

        for(int xi = 0; xi < magfield.GetXCells(); xi++)
        for(int yi = 0; yi < magfield.GetYCells(); yi++)
        for(int zi = 0; zi < magfield.GetZCells(); zi++)
        {
            const Vector& v = magfield.at(xi,yi,zi);
            const double m2cm = 100;
            double x = magfield.GetXMin() + magfield.GetStep()*xi;
            double y = magfield.GetYMin() + magfield.GetStep()*yi;
            double z = magfield.GetZMin() + magfield.GetStep()*zi;

            if (IsInSector(m2cm*x,m2cm*y,m2cm*z,dist))
            {
                if (v.SqMagnitude() < min_magnitude_sq || min_magnitude_sq == -1) 
                    min_magnitude_sq = v.SqMagnitude();
                if (v.SqMagnitude() > max_magnitude_sq)
                    max_magnitude_sq = v.SqMagnitude();
            }            
        }

        out_min = std::sqrt(min_magnitude_sq);
        out_max = std::sqrt(max_magnitude_sq);
    }

    /// @brief Returns minimal (magnetic) field angle to electric field (along positive y-axis) inside first sector TPC volume.
    /// @param magfield Magnetic field data.
    /// @param out_min Output for minimal angle [rad].
    /// @param out_max Output for maximal angle [rad].
    /// @param dist Specifies the minimal distance from within the TPC walls except the bottom anode wall (makes volume smaller).
    void GetMinMaxFieldAngle(const Field<Vector>& magfield, double& out_min, double& out_max, const double& dist)
    {
        double min_angle = -1;
        double max_angle = -1;

        for(int xi = 0; xi < magfield.GetXCells(); xi++)
        for(int yi = 0; yi < magfield.GetYCells(); yi++)
        for(int zi = 0; zi < magfield.GetZCells(); zi++)
        {
            const Vector& v = magfield.at(xi,yi,zi);
            const double m2cm = 100;
            double x = magfield.GetXMin() + magfield.GetStep()*xi;
            double y = magfield.GetYMin() + magfield.GetStep()*yi;
            double z = magfield.GetZMin() + magfield.GetStep()*zi;


            if (IsInSector(m2cm*x,m2cm*y,m2cm*z,dist))
            {
                using namespace constants;

                if (v.Angle(efield) < min_angle || min_angle == -1) 
                    min_angle = v.Angle(efield);
                if (v.Angle(efield) > max_angle)
                    max_angle = v.Angle(efield);
            }            
        }

        out_min = min_angle;
        out_max = max_angle;
    }

    /// @brief Draws lines around the approximate sensitive area (walls of TPC).
    /// @param yxformat If true, the y-coordinate is drawn on the x-axis and vice versa.
    void DrawTrapezoid(bool yxformat = true)
    {
        using namespace constants;

        TLine* l1;
        TLine* l2;
        TLine* l3;
        TLine* l4;

        if (yxformat)
        {
            l1 = new TLine(-yhigh,xmax, yhigh,xmax);
            l2 = new TLine(-ylow ,xmin, ylow ,xmin);
            l3 = new TLine(-yhigh,xmax,-ylow ,xmin);
            l4 = new TLine( yhigh,xmax, ylow ,xmin);
        }
        else
        {
            l1 = new TLine(xmax,-yhigh,xmax, yhigh);
            l2 = new TLine(xmin,-ylow ,xmin, ylow );
            l3 = new TLine(xmax,-yhigh,xmin,-ylow );
            l4 = new TLine(xmax, yhigh,xmin, ylow );
        }

        l1->Draw("same");
        l2->Draw("same");
        l3->Draw("same");
        l4->Draw("same");
    }
} // namespace X17
