// C++ dependencies
#include <cmath>
#include <stdexcept>

// ROOT dependencies
#include "TBox.h"
#include "TEllipse.h"
#include "TLine.h"

// X17 dependencies
#include "Field.h"
#include "Vector.h"
#include "X17Utilities.h"

namespace X17
{
    //// Useful functions for X17.

    bool IsInSector(double x, double y, double z, double dist)
    {
        using namespace constants;

        if (z < zmin + dist || z > zmax) return false;
        if (x < xmin + dist || x > xmax - dist) return false;

        double dy = dist / std::sqrt(1 + 1 / yxslope); // Change of allowed y (absolute value) for slanted surface.
        if (std::abs(y) + dy > yxslope * x + yintersect) return false;
        
        return true;   
    }
    
    void GetMinMaxField(const Field<Vector>& magfield, double& out_min, double& out_max, double dist)
    {
        double min_magnitude_sq = -1;
        double max_magnitude_sq = -1;

        for(int xi = 0; xi < magfield.GetXCells(); xi++)
        for(int yi = 0; yi < magfield.GetYCells(); yi++)
        for(int zi = 0; zi < magfield.GetZCells(); zi++)
        {
            const Vector& v = magfield.at(xi,yi,zi);
            
            double x = magfield.GetXMin() + magfield.GetStep()*xi;
            double y = magfield.GetYMin() + magfield.GetStep()*yi;
            double z = magfield.GetZMin() + magfield.GetStep()*zi;

            if (IsInSector(x,y,z,dist))
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

    void GetMinMaxFieldAngle(const Field<Vector>& magfield, double& out_min, double& out_max, double dist)
    {
        double min_angle = -1;
        double max_angle = -1;

        for(int xi = 0; xi < magfield.GetXCells(); xi++)
        for(int yi = 0; yi < magfield.GetYCells(); yi++)
        for(int zi = 0; zi < magfield.GetZCells(); zi++)
        {
            const Vector& v = magfield.at(xi,yi,zi);
            
            double x = magfield.GetXMin() + magfield.GetStep()*xi;
            double y = magfield.GetYMin() + magfield.GetStep()*yi;
            double z = magfield.GetZMin() + magfield.GetStep()*zi;


            if (IsInSector(x,y,z,dist))
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
    void DrawTube()
    {
        using namespace constants;

        TEllipse* circle = new TEllipse(0, 0, tube_radius, tube_radius);
        circle->SetFillStyle(0); // Set fill style to transparent
        circle->Draw("same");

    }
    void DrawFirstSectorMagnets(bool yxformat)
    {
        using namespace constants;

        // Magnets are at a 30-degree angle from the x-axis
        double angle = 30.0 * M_PI / 180.0;

        // Top magnet corners
        double x1 = mag_lowxx;
        double y1 = mag_lowxy;
        double x2 = x1 + mag_depth * sin(angle);
        double y2 = y1 - mag_depth * cos(angle);
        double x3 = x2 + mag_width * cos(angle);
        double y3 = y2 + mag_width * sin(angle);
        double x4 = x3 - mag_depth * sin(angle);
        double y4 = y3 + mag_depth * cos(angle);

        // Top magnet
        TLine *t1, *t2, *t3, *t4;

        // Bottom magnet
        TLine *b1, *b2, *b3, *b4;

        if (yxformat)
        {
            t1 = new TLine(y1, x1, y2, x2);
            t2 = new TLine(y2, x2, y3, x3);
            t3 = new TLine(y3, x3, y4, x4);
            t4 = new TLine(y4, x4, y1, x1);

            b1 = new TLine(-y1, x1, -y2, x2);
            b2 = new TLine(-y2, x2, -y3, x3);
            b3 = new TLine(-y3, x3, -y4, x4);
            b4 = new TLine(-y4, x4, -y1, x1);
        }
        else
        {
            t1 = new TLine(x1, y1, x2, y2);
            t2 = new TLine(x2, y2, x3, y3);
            t3 = new TLine(x3, y3, x4, y4);
            t4 = new TLine(x4, y4, x1, y1);

            b1 = new TLine(x1, -y1, x2, -y2);
            b2 = new TLine(x2, -y2, x3, -y3);
            b3 = new TLine(x3, -y3, x4, -y4);
            b4 = new TLine(x4, -y4, x1, -y1);
        }

        for (auto l : {t1, t2, t3, t4, b1, b2, b3, b4})
        {
            l->SetLineColor(kRed);
            l->SetLineWidth(2);
            l->Draw("same");
        }
    }
} // namespace X17
