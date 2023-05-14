// C++ dependencies
#include <iostream>
#include <vector>

// ROOT dependencies
#include "TFile.h"
#include "TLine.h"
#include "TPolyLine3D.h"
#include "TText.h"

// X17 dependencies
#include "Field.h"
#include "PadLayout.h"
#include "Points.h"
#include "Reconstruction.h"
#include "X17Utilities.h"

namespace X17
{
    //// Public DefaultLayout methods.

    double DefaultLayout::GetPadHeight(int i, bool effective) const
    {
        using namespace constants;

        double height = pad_height;
        
        // Special cases for height.
        if (i == 102) height = pad_height2;
        if (i == 124) height = pad_height2;
        if (i == 127) height = pad_height3;

        return height;
    }

    void DefaultLayout::GetPadInfo(int i, int& out_column, int& out_row, double& out_height) const
    {
        using namespace constants;
        if(i < 1 || i > channels) std::cerr << "WARNING: Invalid channel number " << i << " (must be between 1 and " << channels << ").\n";

        int triangle_row1   = 7;  // First row with smaller number of columns.
        int triangle_row2   = 10; // First row breaking the first triangle.
        int part1_channels  = (triangle_row1 - 1) * columns;
        int part12_channels = part1_channels + triangle(columns-1) - triangle(columns-1-triangle_row2+triangle_row1);

        // Part with same length rows (first six rows).
        if (i <= part1_channels)  {out_row = (i - 1) / columns + 1; out_column = i % columns;}

        // First triangular part (rows 7-9).
        else if (i <= part12_channels)
        {
            // Maximal number in row r is diference of c-th and (c-r)th triangular number.
            out_row    = triangle_row1 - 1 + triangle_row(columns - 1, i - part1_channels);

            // Column number is the offset of channel from difference of c-th and r-th triangle number.
            out_column = (i - part1_channels) - triangle(columns - 1) + triangle(columns - 1 + triangle_row1 - out_row);
        }

        // Second triangular part (rows 10-15, missing number in last row does not matter).
        else
        {
            int max_cols_p3 = columns - 2 - triangle_row2 + triangle_row1;
            out_row    = triangle_row2 - 1 + triangle_row(max_cols_p3, i - part12_channels);
            out_column = (i - part12_channels) - triangle(columns - 2 - triangle_row2 + triangle_row1) + triangle(max_cols_p3 + triangle_row2 - out_row);
        }

        if (out_column == 0) out_column = 12; // Number 0 corresponds to 12th column.

        out_height = GetPadHeight(i);
    }
    void DefaultLayout::GetPadCenter(int i, double& out_x, double& out_y) const
    {
        using namespace constants;

        int column; // Column number of given pad (0 corresponds to 12).
        int row;    // Row (diagonal) number of given pad.

        double i_pad_height; // Height of the i-th pad.
        
        GetPadInfo(i,column,row,i_pad_height);

        out_x = xmax - (column - 1) * (pad_width + pad_offset) - pad_width / 2;
        out_y = -(yhigh + pad1_offset - (row - 1) * (pad_height + pad_offset) - (column - 1) * pad_stag - i_pad_height / 2);
    }

    void DefaultLayout::GetPadCorners(int i, double& out_xlow, double& out_ylow, double& out_xhigh, double& out_yhigh, bool nogaps) const
    {   
        using namespace constants;

        double gaps_correction = 0; // Removes gaps if nogaps is true.
        if (nogaps) gaps_correction = pad_offset / 2.0;

        int column; // Column number of given pad (0 corresponds to 12).
        int row;    // Row (diagonal) number of given pad.

        double i_pad_height; // Height of the i-th pad.
        
        GetPadInfo(i,column,row,i_pad_height);

        out_xhigh = xmax - (column - 1) * (pad_width + pad_offset) + gaps_correction;
        out_ylow  = -(yhigh + pad1_offset - (row - 1) * (pad_height + pad_offset) - (column - 1) * pad_stag + gaps_correction);
        out_xlow  = out_xhigh - pad_width - 2 * gaps_correction;
        out_yhigh = out_ylow - (-i_pad_height - 2 * gaps_correction);
    }

    int DefaultLayout::GetPad(double x, double y) const
    {
        using namespace constants;

        for (int i = 1; i <= channels; i++)
        {
            double xlow,ylow,xhigh,yhigh;
            GetPadCorners(i,xlow,ylow,xhigh,yhigh,true);
            if(x >= xlow && x < xhigh && y >= ylow && y < yhigh) return i;
        }

        return -1;        
    }

    void DefaultLayout::DrawPads(bool nogaps, TCanvas* c) const
    {
        using namespace constants;

        if(c == nullptr) c = new TCanvas("c_pads","GEM readout pads",600,600*2*(yhigh+1)/(xmax-xmin+2));
        std::vector<TLine*> pad_lines;
        std::vector<TText*> pad_numbers;

        for (int i = 1; i <= channels; i++)
        {
            double x1,y1,x2,y2;
            GetPadCorners(i,x1,y1,x2,y2,nogaps);

            pad_lines.push_back(new TLine(x1,y1,x1,y2)); // Left line.
            pad_lines.push_back(new TLine(x2,y1,x2,y2)); // Right line.
            pad_lines.push_back(new TLine(x1,y1,x2,y1)); // Bottom line.
            pad_lines.push_back(new TLine(x1,y2,x2,y2)); // Top line.

            // double xc,yc;
            // GetPadCenter(i,xc,yc);
            // pad_lines.push_back(new TLine(x1,y1,xc,yc));
            // pad_lines.push_back(new TLine(x1,y2,xc,yc));
            // pad_lines.push_back(new TLine(x2,y1,xc,yc));
            // pad_lines.push_back(new TLine(x2,y2,xc,yc));

            TText* pad_number = new TText((x1+x2)/2,(y1+y2)/2,std::to_string(i).c_str());
            pad_number->SetTextAlign(22);
            pad_number->SetTextSize(0.035);
            pad_numbers.push_back(pad_number);
        }
        
        c->Range(xmin-1,-yhigh-1,xmax+1,yhigh+1);
        for(auto l : pad_lines)   l->Draw("AL");
        for(auto t : pad_numbers) t->Draw("same");

        DrawTrapezoid(false);
    }

    void DefaultLayout::DrawPads3D(const double height) const
    {
        using namespace constants;

        std::vector<TPolyLine3D*> pad_lines;

        for (int i = 1; i <= channels; i++)
        {
            double x1,y1,x2,y2;
            GetPadCorners(i,x1,y1,x2,y2,true);

            TPolyLine3D* pad = new TPolyLine3D(5);
            
            pad->SetPoint(0,x1,y1,height);
            pad->SetPoint(1,x1,y2,height);
            pad->SetPoint(2,x2,y2,height);
            pad->SetPoint(3,x2,y1,height);
            pad->SetPoint(4,x1,y1,height);

            pad_lines.push_back(pad);
        }
        
        for(auto l : pad_lines)   l->Draw("AL");
    }

    void DefaultLayout::DrawPadsDistortion(double time, TCanvas* c, Field<MapPoint>* map) const
    {
        using namespace constants;

        // Get map from initial to final electron positions from file.
        if (map == nullptr)
        {
            TFile* inFile2 = new TFile("/home/vavrik/work/X17/data/ion_map/map.root");
            map = (Field<MapPoint>*)inFile2->Get("map");
        }

        std::string c_title = "GEM readout pads distorted " + std::to_string(time) + " ns";

        if (c == nullptr) 
        {
            c = new TCanvas("c_pads_distorted",c_title.c_str(),600*2*(yhigh+1)/(xmax-xmin+2),600);
            c->Range(-yhigh-1,xmin-1,yhigh+1,xmax+1);
        }
        // else
        // {
        //     c->GetSelectedPad()->Range(-yhigh-1,xmin-1,yhigh+1,xmax+1);
        // }
        std::vector<TLine*> pad_lines;

        for (int i = 1; i <= channels; i++)
        {
            double x1,y1,x2,y2;
            GetPadCorners(i,x1,y1,x2,y2);

            RecoPoint bottomleft  = Reconstruct(*map,EndPoint{x1,y1,0,time});
            RecoPoint bottomright = Reconstruct(*map,EndPoint{x1,y2,0,time});
            RecoPoint topleft     = Reconstruct(*map,EndPoint{x2,y1,0,time});
            RecoPoint topright    = Reconstruct(*map,EndPoint{x2,y2,0,time});

            pad_lines.push_back(new TLine(-bottomleft.y,  bottomleft.x,  -topleft.y,     topleft.x));     // Left line.
            pad_lines.push_back(new TLine(-bottomright.y, bottomright.x, -topright.y,    topright.x));    // Right line.
            pad_lines.push_back(new TLine(-bottomleft.y,  bottomleft.x,  -bottomright.y, bottomright.x)); // Bottom line.
            pad_lines.push_back(new TLine(-topleft.y,     topleft.x,     -topright.y,    topright.x));    // Top line.
        }
        
        for(auto l : pad_lines) l->Draw("same");
        DrawTrapezoid();
    }
} // namespace X17
