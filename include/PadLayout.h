#pragma once

// C++ dependencies
#include <cmath>
#include <string>

// ROOT dependencies
#include "TCanvas.h"
#include "TFile.h"
#include "TPolyLine3D.h"
#include "TText.h"

// X17 dependencies
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
        double GetPadHeight(const int& i, const bool& effective = false);

        /// @brief Function for retrieving column, row and height of the i-th pad.
        /// @param i The number of the pad (channel) between one and the number of pads.
        /// @param out_column Will be set to the column number of the pad.
        /// @param out_row Will be set to the row (diagonal) number of the pad.
        /// @param out_height Will be set to the height of the pad.
        void GetPadInfo(const int& i, int& out_column, int& out_row, double& out_height);

        /// @brief Returns coordinates of the center of the i-th pad.
        /// @param i The number of the pad (channel) between 1 and the 128.
        /// @param out_x Will be set to the x-coordinate of the center.
        /// @param out_y Will be set to the y-coordinate of the center.
        void GetPadCenter(const int& i, double& out_x, double& out_y);

        /// @brief Returns coordinates of the top right and the bottom left corner of i-th pad.
        /// @param i The number of the pad (channel) between 1 and the 128.
        /// @param out_xlow Will be set to the x-coordinate of the bottom left corner.
        /// @param out_ylow Will be set to the y-coordinate of the bottom left corner.
        /// @param out_xhigh Will be set to the x-coordinate of the top right corner.
        /// @param out_yhigh Will be set to the y-coordinate of the top right corner.
        /// @param nogaps If true, the gaps are evenly divided between neighbouring pads.
        void GetPadCorners(const int& i, double& out_xlow, double& out_ylow, double& out_xhigh, double& out_yhigh, bool nogaps = false);

        /// @brief Function for retrieving the number of a pad (channel) effective area hit for given position (currently slow approach, may need optimization in future).
        /// @param x The x-coordinate [cm].
        /// @param y The y-coordinate [cm].
        /// @return Number of pad containing the given point (-1 if no pad contains this point).
        int GetPad(const double& x, const double& y);

        /// @brief Test function for drawing of the pads with their channel numbers.
        /// @param nogaps If true, the gaps are evenly divided between neighbouring pads.
        void DrawPads(bool nogaps = false, TCanvas* c = nullptr);

        /// @brief Function for drawing of the pads in 3D with certain height.
        /// @param height Specifies how high should the pads be drawn.
        void DrawPads3D(const double height = -8);

        /// @brief Draw the pads using the coordinates of the electrons ending up in the corners.
        /// @param time Time to propagate backwards [ns].
        void DrawPadsDistortion(const double& time, TCanvas* c = nullptr, Field<MapPoint>* map = nullptr);
    };
} // namespace X17
