// C++ dependencies
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

// X17 dependencies
#include "Field.h"
#include "Matrix.h"
#include "Points.h"
#include "Reconstruction.h"

namespace X17
{
    //// Functions related to reconstruction.

    template<typename T>
    std::array<T,3> CubeCorner(const T (&indices)[6], int corner)
    {
        T x,y,z;

        if(((int)corner/4)     == 0) x = indices[0];
        else                         x = indices[1];

        if(((int)corner/2) % 2 == 0) y = indices[2];
        else                         y = indices[3];

        if(corner % 2 == 0)          z = indices[4];
        else                         z = indices[5];

        return {x,y,z};
    }

    std::array<std::array<double,8>,3> GetInterpolCoef(const Field<MapPoint>& map, const int (&indices)[6])
    {
        double bounds[6]; // Array {xmin, xmax, ymin, ymax, zmin, zmax}.

        bounds[0] = indices[0] * map.GetStep() + map.GetXMin();
        bounds[1] = indices[1] * map.GetStep() + map.GetXMin();

        bounds[2] = indices[2] * map.GetStep() + map.GetYMin();
        bounds[3] = indices[3] * map.GetStep() + map.GetYMin();

        bounds[4] = indices[4] * map.GetStep() + map.GetZMin();
        bounds[5] = indices[5] * map.GetStep() + map.GetZMin();

        const int coor_ids[] = {0, 1, 3};             // The x, y and t coordinate indices in SensorData (operator[]).
        std::array<std::array<double,8>,3> xytvalues; // An array of size 3, containing arrays of the x, y, and t values at each corner of the pseudocube used for interpolation.
        std::array<std::array<double,8>,3> output;

        // Filling xytvalues.
        for (int j = 0; j < 3; j++)
        {
            for (int i = 0; i < 8; i++)
            {
                const auto [xi,yi,zi] = CubeCorner(indices,i);
                xytvalues[j][i] = map.at(xi,yi,zi)[coor_ids[j]];
            }
        }

        // Solving system of linear equation for each coordinate to get the interpolation coefficients.
        for (int i = 0; i < 3; i++)
        {
            double arr[8*9];

            // Adding rows to the matrix for each cube corner a + b*x + c*y + d*t + e*xy + f*xt + g*yt + h*xyt = (x0,y0,z0).
            for (int j = 0; j < 8; j++)
            {
                arr[9*j]   = 1;
                arr[9*j+1] = xytvalues[0][j];
                arr[9*j+2] = xytvalues[1][j];
                arr[9*j+3] = xytvalues[2][j];
                arr[9*j+4] = xytvalues[0][j] * xytvalues[1][j];
                arr[9*j+5] = xytvalues[0][j] * xytvalues[2][j];
                arr[9*j+6] = xytvalues[1][j] * xytvalues[2][j];
                arr[9*j+7] = xytvalues[0][j] * xytvalues[1][j] * xytvalues[2][j];

                const auto corner = CubeCorner(bounds,j);
                arr[9*j+8] = corner[i];
            }
            
            Matrix<8,9> matrix = Matrix<8,9>(arr);
            matrix.Reduce();
            output[i] = matrix.GetColumn(8);
        }

        return output;
    }

    void FindMapMinMaxIndex(const Field<MapPoint>& map, int& i_min, int& i_max, int (&i_mid)[3], int var, double value)
    {
        // bool i_mid_is_max = false; // Does midpoint index contain higher value?

        // The z-index is connected to variable t.
        int var2 = (var == 2) ? 3 : var;

        // Z search is inverted so the stopping condition has to be inverted
        auto should_continue = [&i_min, &i_max, &var](){
            if (var == 2) return (i_max + 1) < i_min;
            else return (i_min + 1) < i_max;
        };

    #ifdef DEBUG
        // Indices for low and high values.
        int i_low[3], i_high[3];
        std::copy(std::begin(i_mid), std::end(i_mid), std::begin(i_low));
        std::copy(std::begin(i_mid), std::end(i_mid), std::begin(i_high));
        if(!should_continue()) std::cerr << "WARNING: Binary search (" << var << ") skipped.\n";
    #endif

        while (should_continue())
        {
            i_mid[var] = (i_min + i_max) / 2;

            double mid = map.at(i_mid)[var2];

            if(value <= mid) { i_max = i_mid[var]; }//i_mid_is_max = true; }
            if(value >  mid) { i_min = i_mid[var]; }//i_mid_is_max = false; }

        #ifdef DEBUG
            i_low[var]  = i_min;
            i_high[var] = i_max;
            double low  = map.at(i_low)[var2];
            double high = map.at(i_high)[var2];
            if(mid < low || mid > high) std::cerr << "WARNING: Ordering mistake in binary search.\n";
        #endif
        }

        // Set up actual min/max variables.
        // if(i_mid_is_max) { i_max = i_mid[var]; i_min = i_max - 1; }
        // else             { i_min = i_mid[var]; i_max = i_min + 1; }
    }
    
    void PrintCube(const Field<MapPoint>& map, const int (&indices)[6], EndPoint end_point)
    {
        std::cout << "\n-------------------------------------------------------------------------------\n";
        std::cout << "Cube for (x1,y1,t1) = (" << end_point.x() << ", " << end_point.y() << ", " << end_point.t << "): \n";
        for (int i = 0; i < 8; i++)
        {
            const auto [ci,cj,ck] = CubeCorner({0, 1, 0, 1, 0, 1}, i);
            const auto [xi,yi,zi] = CubeCorner(indices, i);
            std::cout << ci << cj << ck << ": [" << xi << "][" << yi << "][" << zi << "],";
            std::cout << " (x,y,t) = (" << map.at(xi,yi,zi).point.x() << ", " << map.at(xi,yi,zi).point.y() << ", " << map.at(xi,yi,zi).point.t << ")\n";
        }

        std::cout << "\n";

        int ximin = indices[0], ximax = indices[1];
        int yimin = indices[2], yimax = indices[3];
        int zimin = indices[4], zimax = indices[5];

        // Bound checks.
        for (int xi : {ximin, ximax})
        for (int yi : {yimin, yimax})
        for (int zi : {zimin, zimax})
        {
            const EndPoint& mappoint = map.at(xi, yi, zi).point;

            if (xi == ximin && end_point.x() < mappoint.x())
                std::cerr << "INFO: Minimal x bound not minimal.\n";
            if (xi == ximax && end_point.x() > mappoint.x())
                std::cerr << "INFO: Maximal x bound not maximal.\n";

            if (yi == yimin && end_point.y() < mappoint.y())
                std::cerr << "INFO: Minimal y bound not minimal.\n";
            if (yi == yimax && end_point.y() > mappoint.y())
                std::cerr << "INFO: Maximal y bound not maximal.\n";

            if (zi == zimin && end_point.t < mappoint.t)
                std::cerr << "INFO: Minimal z bound not minimal.\n";
            if (zi == zimax && end_point.t > mappoint.t)
                std::cerr << "INFO: Maximal z bound not maximal.\n";
        }
        std::cout << "-------------------------------------------------------------------------------\n\n";
    }

    RecoPoint Reconstruct(const Field<MapPoint>& map, EndPoint end_point)
    {
        const double x1 = end_point.x();
        const double y1 = end_point.y();
        const double t1 = end_point.t;

        // Find 8 closest points using binary search (assuming ordering).
        int ximin = 0;                   // Starting minimal search x index.
        int ximax = map.GetXCells() - 1; // Starting maximal search x index.
        int yimin = 0;                   // Starting minimal search y index.
        int yimax = map.GetYCells() - 1; // Starting maximal search y index.
        int zimin = map.GetZCells() - 1; // Starting minimal search z index (inverted - time is maximal for z minimal).
        int zimax = 0;                   // Starting maximal search z index (inverted - time is maximal for z minimal).

        int i_mid[3] = { (ximin + ximax) / 2, (yimin + yimax) / 2, (zimin + zimax) / 2}; // Midpoint array.

        for (int i = 0; i < 2; i++)
        {
            ximin = 0;
            ximax = map.GetXCells() - 1;

            FindMapMinMaxIndex(map, ximin, ximax, i_mid, 0, x1);

            // Some part of the field is initialized with zeros and isn't therefore ordered.
            if(i_mid[0] != 0)
            {
                while(map.at(i_mid[0], yimin, i_mid[2]).point.y() == 0) yimin++;
                while(map.at(i_mid[0], yimax, i_mid[2]).point.y() == 0) yimax--;
            }

            FindMapMinMaxIndex(map, yimin, yimax, i_mid, 1, y1);
            FindMapMinMaxIndex(map, zimin, zimax, i_mid, 2, t1);
        }

        // Sanity check.
        // PrintCube(map, { ximin, ximax, yimin, yimax, zimin, zimax }, end_point);

        // Interpolation from map.
        std::array<std::array<double,8>,3> coef = GetInterpolCoef(map, {ximin, ximax, yimin, yimax, zimin, zimax});

        double xout = coef[0][0] + coef[0][1]*x1 + coef[0][2]*y1 + coef[0][3]*t1 + coef[0][4]*x1*y1 + coef[0][5]*x1*t1 + coef[0][6]*y1*t1 + coef[0][7]*x1*y1*t1;
        double yout = coef[1][0] + coef[1][1]*x1 + coef[1][2]*y1 + coef[1][3]*t1 + coef[1][4]*x1*y1 + coef[1][5]*x1*t1 + coef[1][6]*y1*t1 + coef[1][7]*x1*y1*t1;
        double zout = coef[2][0] + coef[2][1]*x1 + coef[2][2]*y1 + coef[2][3]*t1 + coef[2][4]*x1*y1 + coef[2][5]*x1*t1 + coef[2][6]*y1*t1 + coef[2][7]*x1*y1*t1;

        return RecoPoint(xout, yout, zout, 1);
    }

    double Offset(MapPoint p, double x1, double y1, double t1)
    {
        constexpr double tfact = 0.00327; // Time is measured at different scale, it needs weight.
        return sqrt(pow(x1 - p.point.x(), 2) + pow(y1 - p.point.y(), 2) + pow(tfact * (t1 - p.point.t), 2));
    }

    RecoPoint ReconstructOld(const Field<MapPoint>& map, double x1, double y1, double t1, double max_err)
    {
        // Start looking at the same position.
        double x = x1;
        double y = y1;
        double z = (map.GetZMax() + map.GetZMin()) / 2;
        double step = map.GetStep() / 10;

        double offset;       // Metric of distance between points.
        int iterations = 0;  // Number of iterations should not exceed 100.
        double damp = 0.005; // Damping coefficient.

        // Loop for offset minimization.
        do
        {
            // Calculate offset gradient.
            MapPoint xa = map.GetField(x+step,y,z);
            MapPoint xb = map.GetField(x-step,y,z);
            MapPoint ya = map.GetField(x,y+step,z);
            MapPoint yb = map.GetField(x,y-step,z);
            MapPoint za = map.GetField(x,y,z+step);
            MapPoint zb = map.GetField(x,y,z-step);

            double oxa = Offset(xa,x1,y1,t1);
            double oxb = Offset(xb,x1,y1,t1);
            double oya = Offset(ya,x1,y1,t1);
            double oyb = Offset(yb,x1,y1,t1);
            double oza = Offset(za,x1,y1,t1);
            double ozb = Offset(zb,x1,y1,t1);

            double gradx = (oxa-oxb)/(2*step);
            double grady = (oya-oyb)/(2*step);
            double gradz = (oza-ozb)/(2*step);

            // Adjust current guess by minus gradient.
            x -= damp*gradx; y -= damp*grady; z -= damp*gradz;

            // Check bounds.
            if (x < map.GetXMin()) x = map.GetXMin();
            if (x > map.GetXMax()) x = map.GetXMax();
            if (y < map.GetYMin()) y = map.GetYMin();
            if (y > map.GetYMax()) y = map.GetYMax();
            if (z < map.GetZMin()) z = map.GetZMin();
            if (z > map.GetZMax()) z = map.GetZMax();
            // cout << "gradx: " << x << " grady: " << y << " gradz: " << z << "\n";

            // Calculate values at current position.
            MapPoint cur = map.GetField(x,y,z);
            offset = Offset(cur,x1,y1,t1);

            // Make sure step isn't too high.
            if(offset < 10 * step) step /= 10;

            iterations++;
            // cout << "iter: " << iterations << "\n";        
            if (iterations == 1000) std::cout << "1000 iterations.\n";
        }
        while ((offset > max_err) && (iterations < 1000));

        return RecoPoint(x,y,z,1);
    }
} // namespace X17