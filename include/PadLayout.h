#pragma once

#include <cmath>
#include <string>

#include "TCanvas.h"
#include "TFile.h"
#include "TPolyLine3D.h"
#include "TText.h"

#include "Points.h"
#include "Reconstruction.h"
#include "X17Utilities.h"

namespace X17
{
    /// @brief Abstract class representing a pad layout for the TPC detector.
    class PadLayout
    {
    public:
        /// @brief Returns the height of i-th pad.
        /// @param i The number of the pad (channel) between one and the number of pads.
        /// @param effective If true, gap width will be added to the height.
        /// @return The height of the i-th pad.
        virtual double GetPadHeight(const int& i, const bool& effective = false) = 0;

        /// @brief Function for retrieving column, row and height of the i-th pad.
        /// @param i The number of the pad (channel) between one and the number of pads.
        /// @param out_column Will be set to the column number of the pad.
        /// @param out_row Will be set to the row (diagonal) number of the pad.
        /// @param out_height Will be set to the height of the pad.
        virtual void GetPadInfo(const int& i, int& out_column, int& out_row, double& out_height) = 0;

        /// @brief Returns coordinates of the center of the i-th pad.
        /// @param i The number of the pad (channel) between one and the number of pads.
        /// @param out_x Will be set to the x-coordinate of the center.
        /// @param out_y Will be set to the y-coordinate of the center.
        virtual void GetPadCenter(const int& i, double& out_x, double& out_y) = 0;
        
        /// @brief Returns coordinates of the top right and the bottom left corner of i-th pad.
        /// @param i The number of the pad (channel) between one and the number of pads.
        /// @param out_xlow Will be set to the x-coordinate of the bottom left corner.
        /// @param out_ylow Will be set to the y-coordinate of the bottom left corner.
        /// @param out_xhigh Will be set to the x-coordinate of the top right corner.
        /// @param out_yhigh Will be set to the y-coordinate of the top right corner.
        /// @param nogaps If true, the gaps are evenly divided between neighbouring pads.
        virtual void GetPadCorners(const int& i, double& out_xlow, double& out_ylow, double& out_xhigh, double& out_yhigh, bool nogaps = false) = 0;

        /// @brief Function for retrieving the number of a pad (channel) effective area hit for given position (currently slow approach, may need optimization in future).
        /// @param x The x-coordinate [cm].
        /// @param y The y-coordinate [cm].
        /// @return Number of pad containing the given point (-1 if no pad contains this point).
        virtual int GetPad(const double& x, const double& y) = 0;

        /// @brief Test function for drawing of the pads with their channel numbers.
        /// @param nogaps If true, the gaps are evenly divided between neighbouring pads.
        virtual void DrawPads(bool nogaps = false, TCanvas* c = nullptr) = 0;

        /// @brief Function for drawing of the pads in 3D with certain height.
        /// @param height Specifies how high should the pads be drawn.
        virtual void DrawPads3D(const double height = -8) = 0;

        /// @brief Draw the pads using the coordinates of the electrons ending up in the corners.
        /// @param time Time to propagate backwards [ns].
        virtual void DrawPadsDistortion(const double& time, TCanvas* c = nullptr, Field<MapPoint>* map = nullptr) = 0;
    };

    /// @brief A function for the n-th triangle number calculation.
    /// @param n Side of the triangle (i.e. the number of columns).
    /// @return The n-th triangle number.
    int triangle(int n) { return n * (n + 1) / 2; }

    /// @brief A function that returns the row of an element in inversed triangle (starts with base).
    /// @param columns The total (maximal) number of columns in the triangle.
    /// @param index The index of an element.
    /// @return The number of row of element with given index.
    int triangle_row(int columns,int index) { return ceil((2 * columns + 1 - sqrt(pow((2 * columns + 1), 2) - 8 * index)) / 2); }

    /// @brief Singleton class for the default pad layout of the TPC detector (the one that is expected to be used).
    class DefaultLayout : public PadLayout
    {
    private:
        /// @brief Default constructor.
        DefaultLayout() = default;

        /// @brief Deleted copy constructor. 
        DefaultLayout(const DefaultLayout&) = delete;

        /// @brief Deleted assignment operator. 
        DefaultLayout& operator=(const DefaultLayout&) = delete;

    public:
        /// @brief Function for retrieving the singleton instance
        /// @return The singleton instance
        static DefaultLayout& GetDefaultLayout()
        {
            static DefaultLayout instance;
            return instance;
        }

        /// @brief Returns the height of i-th pad.
        /// @param i The number of the pad (channel) between 1 and 128.
        /// @param effective If true, gap width will be added to the height.
        /// @return The height of the i-th pad.
        double GetPadHeight(const int& i, const bool& effective = false)
        {
            using namespace constants;

            double height = pad_height;
            
            // Special cases for height
            if (i == 102) height = pad_height2;
            if (i == 124) height = pad_height2;
            if (i == 127) height = pad_height3;

            return height;
        }

        /// @brief Function for retrieving column, row and height of the i-th pad.
        /// @param i The number of the pad (channel) between one and the number of pads.
        /// @param out_column Will be set to the column number of the pad.
        /// @param out_row Will be set to the row (diagonal) number of the pad.
        /// @param out_height Will be set to the height of the pad.
        void GetPadInfo(const int& i, int& out_column, int& out_row, double& out_height)
        {
            using namespace constants;
            if(i < 1 || i > channels) std::cerr << "WARNING: Invalid channel number " << i << " (must be between 1 and " << channels << ").\n";

            int triangle_row1   = 7;  // first row with smaller number of columns
            int triangle_row2   = 10; // first row breaking the first triangle
            int part1_channels  = (triangle_row1 - 1) * columns;
            int part12_channels = part1_channels + triangle(columns-1) - triangle(columns-1-triangle_row2+triangle_row1);

            // Part with same length rows (first six rows)
            if (i <= part1_channels)  {out_row = (i-1)/columns+1; out_column = i%columns;}

            // First triangular part (rows 7-9)
            else if (i <= part12_channels)
            {
                // Maximal number in row r is diference of c-th and (c-r)th triangular number
                out_row    = triangle_row1-1 + triangle_row(columns-1,i-part1_channels);

                // Column number is the offset of channel from difference of c-th and r-th triangle number
                out_column = (i-part1_channels) - triangle(columns-1) + triangle(columns-1+triangle_row1-out_row);
            }

            // Second triangular part (rows 10-15, missing number in last row does not matter)
            else
            {
                int max_cols_p3 = columns-2-triangle_row2+triangle_row1;
                out_row    = triangle_row2-1 + triangle_row(max_cols_p3,i-part12_channels);
                out_column = (i-part12_channels) - triangle(columns-2-triangle_row2+triangle_row1) + triangle(max_cols_p3+triangle_row2-out_row);
            }

            if (out_column == 0) out_column = 12; // Number 0 corresponds to 12th column

            out_height = GetPadHeight(i);
        }

        /// @brief Returns coordinates of the center of the i-th pad.
        /// @param i The number of the pad (channel) between 1 and the 128.
        /// @param out_x Will be set to the x-coordinate of the center.
        /// @param out_y Will be set to the y-coordinate of the center.
        void GetPadCenter(const int& i, double& out_x, double& out_y)
        {
            using namespace constants;

            int column; // Column number of given pad (0 corresponds to 12)
            int row;    // Row (diagonal) number of given pad

            double i_pad_height; // Height of the i-th pad
            
            GetPadInfo(i,column,row,i_pad_height);

            out_x = xmax - (column-1)*(pad_width+pad_offset) - pad_width/2;
            out_y = -(yhigh + pad1_offset - (row-1)*(pad_height + pad_offset) - (column-1)*pad_stag - i_pad_height/2);
        }

        /// @brief Returns coordinates of the top right and the bottom left corner of i-th pad.
        /// @param i The number of the pad (channel) between 1 and the 128.
        /// @param out_xlow Will be set to the x-coordinate of the bottom left corner.
        /// @param out_ylow Will be set to the y-coordinate of the bottom left corner.
        /// @param out_xhigh Will be set to the x-coordinate of the top right corner.
        /// @param out_yhigh Will be set to the y-coordinate of the top right corner.
        /// @param nogaps If true, the gaps are evenly divided between neighbouring pads.
        void GetPadCorners(const int& i, double& out_xlow, double& out_ylow, double& out_xhigh, double& out_yhigh, bool nogaps = false)
        {   
            using namespace constants;

            double gaps_correction = 0; // removes gaps if nogaps is true
            if (nogaps) gaps_correction = pad_offset / 2.0;

            int column; // Column number of given pad (0 corresponds to 12)
            int row;    // Row (diagonal) number of given pad

            double i_pad_height; // Height of the i-th pad
            
            GetPadInfo(i,column,row,i_pad_height);

            out_xhigh = xmax - (column-1)*(pad_width+pad_offset) + gaps_correction;
            out_ylow  = -(yhigh + pad1_offset - (row-1)*(pad_height+pad_offset) - (column-1)*pad_stag + gaps_correction);
            out_xlow  = out_xhigh - pad_width - 2*gaps_correction;
            out_yhigh = out_ylow - (-i_pad_height - 2*gaps_correction);
        }

        /// @brief Function for retrieving the number of a pad (channel) effective area hit for given position (currently slow approach, may need optimization in future).
        /// @param x The x-coordinate [cm].
        /// @param y The y-coordinate [cm].
        /// @return Number of pad containing the given point (-1 if no pad contains this point).
        int GetPad(const double& x, const double& y)
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

        /// @brief Test function for drawing of the pads with their channel numbers.
        /// @param nogaps If true, the gaps are evenly divided between neighbouring pads.
        void DrawPads(bool nogaps = false, TCanvas* c = nullptr)
        {
            using namespace constants;

            if(c == nullptr) c = new TCanvas("c_pads","GEM readout pads",600,600*2*(yhigh+1)/(xmax-xmin+2));
            std::vector<TLine*> pad_lines;
            std::vector<TText*> pad_numbers;

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

        /// @brief Function for drawing of the pads in 3D with certain height.
        /// @param height Specifies how high should the pads be drawn.
        void DrawPads3D(const double height = -8)
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

        /// @brief Draw the pads using the coordinates of the electrons ending up in the corners.
        /// @param time Time to propagate backwards [ns].
        void DrawPadsDistortion(const double& time, TCanvas* c = nullptr, Field<MapPoint>* map = nullptr)
        {
            using namespace constants;

            // Get map from initial to final electron positions from file
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

                pad_lines.push_back(new TLine(-bottomleft.y,  bottomleft.x,  -topleft.y,     topleft.x));     // left
                pad_lines.push_back(new TLine(-bottomright.y, bottomright.x, -topright.y,    topright.x));    // right
                pad_lines.push_back(new TLine(-bottomleft.y,  bottomleft.x,  -bottomright.y, bottomright.x)); // bottom
                pad_lines.push_back(new TLine(-topleft.y,     topleft.x,     -topright.y,    topright.x));    // top
            }
            
            for(auto l : pad_lines) l->Draw("same");
            DrawTrapezoid();
        }
    };
} // namespace X17
