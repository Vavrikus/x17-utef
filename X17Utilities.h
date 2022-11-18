#include <cmath>

#include "VectorField.h"

/// @brief Namespace for variables describing X17 detector
namespace X17
{
    // TPC first sector (contains positive x axis) boundaries
    constexpr double ymin  = -8;    // TPC minimal height [cm]
    constexpr double ymax  =  8;    // TPC maximal height [cm]
    constexpr double xmin  = 6.51;  // TPC minimal x [cm]
    constexpr double xmax  = 14.61; // TPC maximal x [cm]
    constexpr double zlow  = 2.25;  // TPC z of positive corner closer to center [cm]
    constexpr double zhigh = 7.45;  // TPC z of positive corner further from center [cm]

    // TPC first sector derived parameters
    constexpr double zxslope    = (zhigh-zlow)/(xmax-xmin); // TPC Δz/Δx slanted side (z positive)
    constexpr double zintersect = zlow-zxslope*xmin;        // TPC slanted side (z positive) z-axis intersection (extrapolated) [cm]

    // TPC assumed electric field [V/cm]
    constexpr Vector efield = {0,-400,0};

    /// @brief Returns true if given point is in the first sector (containing positive x-axis) of the detector (trapezoidal prism).
    /// @param x x coordinate [cm]
    /// @param y y coordinate [cm]
    /// @param z z coordinate [cm]
    /// @param dist specifies minimal distance from within the TPC walls (makes volume smaller)
    /// @return 
    bool IsInSector(const double& x, const double& y, const double& z, const double& dist = 0)
    {
        if (y < X17::ymin+dist || y > X17::ymax-dist) return false;
        if (x < X17::xmin+dist || x > X17::xmax-dist) return false;
        double dz = dist/std::sqrt(1+1/zxslope); // change of allowed z (absolute value) for slanted surface
        if (std::abs(z)+dz > X17::zxslope*x+X17::zintersect) return false;
        
        return true;   
    }

    /// @brief Returns minimal and maximal (magnetic) field inside first sector TPC volume.
    /// @param bfield Magnetic field
    /// @param out_min Output for minimal field [T]
    /// @param out_max Output for maximal field [T]
    /// @param dist specifies minimal distance from within the TPC walls (makes volume smaller)
    /// @return 
    void GetMinMaxField(const VectorField& bfield, double& out_min, double& out_max, const double& dist = 0)
    {
        double min_magnitude_sq = -1;
        double max_magnitude_sq = -1;

        for(int xi = 0; xi < bfield.ximax; xi++)
        for(int yi = 0; yi < bfield.yimax; yi++)
        for(int zi = 0; zi < bfield.zimax; zi++)
        {
            const Vector& v = bfield.field[xi][yi][zi];
            const double m2cm = 100;
            double x = bfield.xmin + bfield.step*xi;
            double y = bfield.ymin + bfield.step*yi;
            double z = bfield.zmin + bfield.step*zi;

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
    /// @param bfield Magnetic field
    /// @param out_min Output for minimal angle [rad]
    /// @param out_max Output for maximal angle [rad]
    /// @param dist specifies minimal distance from within the TPC walls (makes volume smaller)
    /// @return 
    void GetMinMaxFieldAngle(const VectorField& bfield, double& out_min, double& out_max, const double& dist = 0)
    {
        double min_angle = -1;
        double max_angle = -1;

        for(int xi = 0; xi < bfield.ximax; xi++)
        for(int yi = 0; yi < bfield.yimax; yi++)
        for(int zi = 0; zi < bfield.zimax; zi++)
        {
            const Vector& v = bfield.field[xi][yi][zi];
            const double m2cm = 100;
            double x = bfield.xmin + bfield.step*xi;
            double y = bfield.ymin + bfield.step*yi;
            double z = bfield.zmin + bfield.step*zi;

            if (IsInSector(m2cm*x,m2cm*y,m2cm*z,dist))
            {
                if (v.Angle(X17::efield) < min_angle || min_angle == -1) 
                    min_angle = v.Angle(X17::efield);
                if (v.Angle(X17::efield) > max_angle)
                    max_angle = v.Angle(X17::efield);
            }            
        }

        out_min = min_angle;
        out_max = max_angle;
    }
} // namespace X17