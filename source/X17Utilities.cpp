// ROOT dependencies
#include "TLine.h"

// X17 dependencies
#include "X17Utilities.h"

namespace X17
{    
    bool IsInSector(const double& x, const double& y, const double& z, const double& dist)
    {
        using namespace constants;

        if (z < zmin      || z > zmax-dist) return false;
        if (x < xmin+dist || x > xmax-dist) return false;

        double dy = dist / std::sqrt(1 + 1/yxslope); // change of allowed y (absolute value) for slanted surface
        if (std::abs(y) + dy > yxslope*x + yintersect) return false;
        
        return true;   
    }
    
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

    void DrawTrapezoid(bool yxformat)
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
