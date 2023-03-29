#pragma once

#include <cmath>
#include <iostream>
#include <string.h>

#include "TCanvas.h"
#include "TFile.h"
#include "TLine.h"
#include "TText.h"

#include "VectorField.h"

/// @brief Namespace for variables describing X17 detector and related functions
namespace X17
{
    // Physical constants
        constexpr double e  = 1.602176634E-19; // Elementary charge
        constexpr double c  = 299792458;       // Light velocity in vacuum [m/s]
        constexpr double E0 = 510998.95;       // Electron rest energy [eV]
        constexpr double m0 = e*E0/(c*c);      // Electron rest mass [kg]

    // TPC first sector (contains positive x axis) boundaries
        constexpr double zmin  = -8;    // TPC minimal height [cm]
        constexpr double zmax  =  8;    // TPC maximal height [cm]
        constexpr double xmin  = 6.51;  // TPC minimal x [cm]
        constexpr double xmax  = 14.61; // TPC maximal x [cm]
        constexpr double ylow  = 2.75;  // TPC y of positive corner closer to center [cm]
        constexpr double yhigh = 7.45;  // TPC y of positive corner further from center [cm]

    // TPC first sector derived parameters
        constexpr double yxslope    = (yhigh-ylow)/(xmax-xmin); // TPC Δy/Δx slanted side (y positive)
        constexpr double yintersect = ylow-yxslope*xmin;        // TPC slanted side (y positive) y-axis intersection (extrapolated) [cm]

    // TPC assumed electric field [V/cm]
        constexpr Vector efield = {0,0,-400};

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
    /// @return true/false
    bool IsInSector(const double& x, const double& y, const double& z, const double& dist = 0)
    {
        if (z < X17::zmin      || z > X17::zmax-dist) return false;
        if (x < X17::xmin+dist || x > X17::xmax-dist) return false;
        double dy = dist/std::sqrt(1+1/yxslope); // change of allowed y (absolute value) for slanted surface
        if (std::abs(y)+dy > X17::yxslope*x+X17::yintersect) return false;
        
        return true;   
    }

    /// @brief True if vector is in the first sector (containing positive x-axis) of the detector (trapezoidal prism).
    /// @param vec point in space (x,y,z) [cm]
    /// @param dist specifies minimal distance from within the TPC walls except the bottom anode wall (makes volume smaller)
    /// @return true/false
    bool IsInSector(const Vector& vec, const double& dist = 0)
    {
        return IsInSector(vec.vx,vec.vy,vec.vz,dist);
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
    int triangle_row(int columns,int index) {return ceil((2*columns+1-sqrt(pow((2*columns+1),2)-8*index))/2);}

    /// @brief Returns height of i-th pad.
    /// @param i Number of pad (channel) between 1 and 128
    /// @param effective If true, gap width will be added to the height
    /// @return Height of i-th pad 
    double GetPadHeight(const int& i, bool effective = false)
    {
        double height = pad_height;
        
        // Special cases for height
        if (i == 102) height = pad_height2;
        if (i == 124) height = pad_height2;
        if (i == 127) height = pad_height3;

        return height;
    }

    /// @brief Returns column, row and height of i-th pad (using input parameters)
    /// @param i Number of pad (channel) between 1 and 128
    /// @param out_column Column number of given pad
    /// @param out_row Row (diagonal) number of given pad
    /// @param out_height Height of given pad
    void GetPadInfo(const int& i, int& out_column, int& out_row, double& out_height)
    {
        if(i < 1 || i > channels) std::cerr << "ERROR: Invalid channel number " << i << " (must be between 1 and " << channels << ").";

        int triangle_row1   = 7;  // first row with smaller number of columns
        int triangle_row2   = 10; // first row breaking the first triangle
        int part1_channels  = (triangle_row1-1)*columns;
        int part12_channels = part1_channels+triangle(columns-1)-triangle(columns-1-triangle_row2+triangle_row1);

        // Part with same length rows (first six rows)
        if (i <= part1_channels)  {out_row = (i-1)/columns+1; out_column = i%columns;}

        // First triangular part (rows 7-9)
        else if (i <= part12_channels)
        {
            out_row = triangle_row1-1+triangle_row(columns-1,i-part1_channels); // Maximal number in row r is diference of c-th and (c-r)th triangular number
            out_column = (i-part1_channels) - triangle(columns-1) + triangle(columns-1+triangle_row1-out_row); // Column number is the offset of channel from difference of c-th and r-th Tn
        }

        // Second triangular part (rows 10-15, missing number in last row does not matter)
        else
        {
            int max_cols_p3 = columns-2-triangle_row2+triangle_row1;
            out_row = triangle_row2-1+triangle_row(max_cols_p3,i-part12_channels);
            out_column = (i-part12_channels) - triangle(columns-2-triangle_row2+triangle_row1) + triangle(max_cols_p3+triangle_row2-out_row);
        }

        if (out_column == 0) out_column = 12; // Number 0 corresponds to 12th column

        out_height = GetPadHeight(i);
    }

    /// @brief Returns coordinates of the center of i-th pad
    /// @param i Number of pad (channel) between 1 and 128
    /// @param out_x Coordinate x of center
    /// @param out_y Coordinate y of center
    void GetPadCenter(const int& i, double& out_x, double& out_y)
    {
        int column; // Column number of given pad (0 corresponds to 12)
        int row;    // Row (diagonal) number of given pad

        double i_pad_height; // Height of the i-th pad
        
        GetPadInfo(i,column,row,i_pad_height);

        out_x = xmax - (column-1)*(pad_width+pad_offset) - pad_width/2;
        out_y = -(yhigh + pad1_offset - (row-1)*(pad_height+pad_offset)-(column-1)*pad_stag - i_pad_height/2);
    }

    /// @brief Returns coordinates of top right and bottom left corners of i-th pad (using input parameters)
    /// @param i Number of pad (channel) between 1 and 128
    /// @param out_xlow Coordinate x of bottom left corner
    /// @param out_ylow Coordinate y of bottom left corner
    /// @param out_xhigh Coordinate x of top right corner
    /// @param out_yhigh Coordinate y of top right corner
    /// @param nogaps If true, gaps between pads are evenly divided between neighbouring pads
    void GetPadCorners(const int& i, double& out_xlow, double& out_ylow, double& out_xhigh, double& out_yhigh, bool nogaps = false)
    {        
        double gaps_correction = 0; // removes gaps if nogaps is true
        if(nogaps) gaps_correction = pad_offset/2.0;

        int column; // Column number of given pad (0 corresponds to 12)
        int row;    // Row (diagonal) number of given pad

        double i_pad_height; // Height of the i-th pad
        
        GetPadInfo(i,column,row,i_pad_height);

        out_xhigh = xmax - (column-1)*(pad_width+pad_offset) + gaps_correction;
        out_ylow = -(yhigh + pad1_offset - (row-1)*(pad_height+pad_offset)-(column-1)*pad_stag + gaps_correction);
        out_xlow  = out_xhigh - pad_width - 2*gaps_correction;
        out_yhigh  = out_ylow - (-i_pad_height - 2*gaps_correction);
    }

    /// @brief Function for retrieving the number of a pad (channel) effective area hit for given position (currently slow approach, may need optimization in future)
    /// @param x x coordinate [cm]
    /// @param y y coordinate [cm]
    /// @return Number of pad containing the given point (-1 if no pad contains this point)
    int GetPad(const double& x, const double& y)
    {
        for (int i = 1; i <= channels; i++)
        {
            double xlow,ylow,xhigh,yhigh;
            GetPadCorners(i,xlow,ylow,xhigh,yhigh,true);
            if(x >= xlow && x < xhigh && y >= ylow && y < yhigh) return i;
        }

        return -1;        
    }

    /// @brief Draws lines around approximate sensitive area (walls of TPC)
    /// @param yxformat If true, y coordinate is drawn on x axis and vice versa
    void DrawTrapezoid(bool yxformat = true)
    {
        TLine* l1;
        TLine* l2;
        TLine* l3;
        TLine* l4;

        if (yxformat)
        {
            l1 = new TLine(-X17::yhigh,X17::xmax, X17::yhigh,X17::xmax);
            l2 = new TLine(-X17::ylow ,X17::xmin, X17::ylow ,X17::xmin);
            l3 = new TLine(-X17::yhigh,X17::xmax,-X17::ylow ,X17::xmin);
            l4 = new TLine( X17::yhigh,X17::xmax, X17::ylow ,X17::xmin);
        }
        else
        {
            l1 = new TLine(X17::xmax,-X17::yhigh,X17::xmax, X17::yhigh);
            l2 = new TLine(X17::xmin,-X17::ylow ,X17::xmin, X17::ylow );
            l3 = new TLine(X17::xmax,-X17::yhigh,X17::xmin,-X17::ylow );
            l4 = new TLine(X17::xmax, X17::yhigh,X17::xmin, X17::ylow );
        }

        l1->Draw("same");
        l2->Draw("same");
        l3->Draw("same");
        l4->Draw("same");
    }

    /// @brief Test function for drawing the pads with their channel numbers
    /// @param nogaps If true, gaps between pads are evenly divided between neighbouring pads
    void DrawPads(bool nogaps = false, TCanvas* c = nullptr)
    {
        if(c == nullptr) c = new TCanvas("c_pads","GEM readout pads",600,600*2*(yhigh+1)/(xmax-xmin+2));
        vector<TLine*> pad_lines;
        vector<TText*> pad_numbers;

        for (int i = 1; i <= channels; i++)
        {
            double x1,y1,x2,y2;
            GetPadCorners(i,x1,y1,x2,y2,nogaps);

            pad_lines.push_back(new TLine(x1,y1,x1,y2)); // left
            pad_lines.push_back(new TLine(x2,y1,x2,y2)); // right
            pad_lines.push_back(new TLine(x1,y1,x2,y1)); // bottom
            pad_lines.push_back(new TLine(x1,y2,x2,y2)); // top

            // double xc,yc;
            // GetPadCenter(i,xc,yc);
            // pad_lines.push_back(new TLine(x1,y1,xc,yc));
            // pad_lines.push_back(new TLine(x1,y2,xc,yc));
            // pad_lines.push_back(new TLine(x2,y1,xc,yc));
            // pad_lines.push_back(new TLine(x2,y2,xc,yc));

            TText* pad_number = new TText((x1+x2)/2,(y1+y2)/2,to_string(i).c_str());
            pad_number->SetTextAlign(22);
            pad_number->SetTextSize(0.035);
            pad_numbers.push_back(pad_number);
        }
        
        c->Range(xmin-1,-yhigh-1,xmax+1,yhigh+1);
        for(auto l : pad_lines)   l->Draw("AL");
        for(auto t : pad_numbers) t->Draw("same");

        DrawTrapezoid(false);
    }


    /// @brief Draw pads using coordinates of electrons ending up in the corners
    /// @param time time to propagate [ns]
    void DrawPadsDistortion(const double& time, TCanvas* c = nullptr,Field<SensorData>* map = nullptr)
    {
        // Get map from initial to final electron positions from file
        if (map == nullptr)
        {
            TFile* inFile2 = new TFile("map.root");
            map = (Field<SensorData>*)inFile2->Get("map");
        }

        string c_title = "GEM readout pads distorted " + to_string(time) + " ns";

        if (c == nullptr) 
        {
            c = new TCanvas("c_pads_distorted",c_title.c_str(),600*2*(yhigh+1)/(xmax-xmin+2),600);
            c->Range(-yhigh-1,xmin-1,yhigh+1,xmax+1);
        }
        // else
        // {
        //     c->GetSelectedPad()->Range(-yhigh-1,xmin-1,yhigh+1,xmax+1);
        // }
        vector<TLine*> pad_lines;

        for (int i = 1; i <= channels; i++)
        {
            double x1,y1,x2,y2;
            GetPadCorners(i,x1,y1,x2,y2);

            SensorData bottomleft  = map->Invert(x1,y1,time);//RecoPoint(map,x1,y1,5000,0.01);
            SensorData bottomright = map->Invert(x1,y2,time);//RecoPoint(map,x2,y1,5000,0.01);
            SensorData topleft     = map->Invert(x2,y1,time);//RecoPoint(map,x1,y2,5000,0.01);
            SensorData topright    = map->Invert(x2,y2,time);//RecoPoint(map,x2,y2,5000,0.01);

            pad_lines.push_back(new TLine(-bottomleft.y1,  bottomleft.x1,  -topleft.y1,     topleft.x1));     // left
            pad_lines.push_back(new TLine(-bottomright.y1, bottomright.x1, -topright.y1,    topright.x1));    // right
            pad_lines.push_back(new TLine(-bottomleft.y1,  bottomleft.x1,  -bottomright.y1, bottomright.x1)); // bottom
            pad_lines.push_back(new TLine(-topleft.y1,     topleft.x1,     -topright.y1,    topright.x1));    // top
        }
        
        for(auto l : pad_lines) l->Draw("same");
        DrawTrapezoid();
    }
} // namespace X17