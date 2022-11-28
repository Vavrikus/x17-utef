#include <cmath>
#include <iostream>
#include <string.h>

#include "TCanvas.h"
#include "TLine.h"
#include "TText.h"

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

    /// @brief Returns true if given point is in the first sector (containing positive x-axis) of the detector (trapezoidal prism).
    /// @param x x coordinate [cm]
    /// @param y y coordinate [cm]
    /// @param z z coordinate [cm]
    /// @param dist specifies minimal distance from within the TPC walls except the bottom anode wall (makes volume smaller)
    /// @return 
    bool IsInSector(const double& x, const double& y, const double& z, const double& dist = 0)
    {
        if (y < X17::ymin      || y > X17::ymax-dist) return false;
        if (x < X17::xmin+dist || x > X17::xmax-dist) return false;
        double dz = dist/std::sqrt(1+1/zxslope); // change of allowed z (absolute value) for slanted surface
        if (std::abs(z)+dz > X17::zxslope*x+X17::zintersect) return false;
        
        return true;   
    }

    /// @brief Returns minimal and maximal (magnetic) field inside first sector TPC volume.
    /// @param bfield Magnetic field
    /// @param out_min Output for minimal field [T]
    /// @param out_max Output for maximal field [T]
    /// @param dist specifies minimal distance from within the TPC walls except the bottom anode wall (makes volume smaller)
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
    /// @param dist specifies minimal distance from within the TPC walls except the bottom anode wall (makes volume smaller)
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

    /// @brief Triangle number
    /// @param n Side of triangle
    /// @return n-th triangle number
    int triangle(int n) {return n*(n+1)/2;}

    /// @brief Triangle with base as the first row
    /// @param columns The total number of columns in the triangle
    /// @param index The index of element
    /// @return The number of row of element with given index 
    int triangle_row(int columns,int index){return ceil((2*columns+1-sqrt(pow((2*columns+1),2)-8*index))/2);}

    /// @brief Returns coordinates of top right and bottom left corners of ith pad
    /// @param i Number of pad (channel) between 1 and 128
    /// @param out_xlow Coordinate x of bottom left corner
    /// @param out_zlow Coordinate z of bottom left corner
    /// @param out_xhigh Coordinate x of top right corner
    /// @param out_zhigh Coordinate z of top right corner
    void GetPadCorners(const int& i, double& out_xlow, double& out_zlow, double& out_xhigh, double& out_zhigh)
    {
        if(i < 1 || i > channels) std::cerr << "ERROR: Invalid channel number " << i << " (must be between 1 and " << channels << ").";
        
        int column; // Column number of given pad (0 corresponds to 12)
        int row;    // Row (diagonal) number of given pad
        
        int triangle_row1   = 7;  // first row with smaller number of columns
        int triangle_row2   = 10; // first row breaking the first triangle
        int part1_channels  = (triangle_row1-1)*columns;
        int part12_channels = part1_channels+triangle(columns-1)-triangle(columns-1-triangle_row2+triangle_row1);

        // Part with same length rows (first six rows)
        if (i <= part1_channels)  {row = (i-1)/columns+1; column = i%columns;}

        // First triangular part (rows 7-9)
        else if (i <= part12_channels)
        {
            row = triangle_row1-1+triangle_row(columns-1,i-part1_channels); // Maximal number in row r is diference of c-th and (c-r)th triangulae number
            column = (i-part1_channels) - triangle(columns-1) + triangle(columns-1+triangle_row1-row); // Column number is the offset of channel from difference of c-th and r-th Tn
        }

        // Second triangular part (rows 10-15, missing number in last row does not matter)
        else
        {
            int max_cols_p3 = columns-2-triangle_row2+triangle_row1;
            row = triangle_row2-1+triangle_row(max_cols_p3,i-part12_channels);
            column = (i-part12_channels) - triangle(columns-2-triangle_row2+triangle_row1) + triangle(max_cols_p3+triangle_row2-row);
        }

        if (column == 0) column = 12; // Number 0 corresponds to 12th column

        double i_pad_height = pad_height;
        
        // Special cases for height
        if (i == 102) i_pad_height = pad_height2;
        if (i == 124) i_pad_height = pad_height2;
        if (i == 127) i_pad_height = pad_height3;

        out_xhigh = xmax - (column-1)*(pad_width+pad_offset);
        out_zhigh = zhigh - (row-1)*(pad_height+pad_offset)-(column-1)*pad_stag;
        out_xlow  = out_xhigh - pad_width;
        out_zlow  = out_zhigh - i_pad_height;
    }

    /// @brief Control function for drawing the pads
    void DrawPads()
    {
        TCanvas* c = new TCanvas("c_pads","GEM readout pads",600,600*2*(zhigh+1)/(xmax-xmin+2));
        vector<TLine*> pad_lines;
        vector<TText*> pad_numbers;

        for (int i = 1; i <= channels; i++)
        {
            double x1,z1,x2,z2;
            GetPadCorners(i,x1,z1,x2,z2);

            pad_lines.push_back(new TLine(x1,z1,x1,z2)); // left
            pad_lines.push_back(new TLine(x2,z1,x2,z2)); // right
            pad_lines.push_back(new TLine(x1,z1,x2,z1)); // bottom
            pad_lines.push_back(new TLine(x1,z2,x2,z2)); // top

            TText* pad_number = new TText((x1+x2)/2,(z1+z2)/2,to_string(i).c_str());
            pad_number->SetTextAlign(22);
            pad_number->SetTextSize(0.035);
            pad_numbers.push_back(pad_number);
        }
        
        c->Range(xmin-1,-zhigh-1,xmax+1,zhigh+1);
        for(auto l : pad_lines)   l->Draw("AL");
        for(auto t : pad_numbers) t->Draw("same");
    }
} // namespace X17