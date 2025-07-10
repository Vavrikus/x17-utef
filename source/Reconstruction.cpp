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

// #define DEBUG

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
        double bounds[6]; // Array {xmin, xmax, ymin, ymax, z_tmin, z_tmax}.

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
        // (i_min for var == 2 is zi_tmin <--> zimax and vice versa).
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

    bool CellContainsReco(const Field<MapPoint> &map, const int (&indices)[6], EndPoint end_point, int& out_status)
    {
        int ximin = indices[0], ximax = indices[1];
        int yimin = indices[2], yimax = indices[3];
        int zi_tmin = indices[4], zi_tmax = indices[5];
        
        bool contains_smaller[3] = { false, false, false };
        bool contains_bigger[3]  = { false, false, false };

        for (int xi : {ximin, ximax})
        for (int yi : {yimin, yimax})
        for (int zi : {zi_tmin, zi_tmax})
        {
            const EndPoint& mappoint = map.at(xi, yi, zi).point;

            if (end_point.x() > mappoint.x())
                contains_smaller[0] = true;
            if (end_point.x() < mappoint.x())
                contains_bigger[0] = true;

            if (end_point.y() > mappoint.y())
                contains_smaller[1] = true;
            if (end_point.y() < mappoint.y())
                contains_bigger[1] = true;

            if (end_point.t > mappoint.t)
            contains_smaller[2] = true;
            if (end_point.t < mappoint.t)
            contains_bigger[2] = true;
        }
        
        bool is_ok = true;
        out_status = -1;
        for (int i = 0; i < 3; i++)
        {
            if (!contains_smaller[i])
            {
                is_ok = false;
                out_status = 2*i;
            }
            else if (!contains_bigger[i])
            {
                is_ok = false;
                out_status = 2*i+1;
            }
        }

        return is_ok;
    }

    void PrintCube(const Field<MapPoint>& map, const int (&indices)[6], EndPoint end_point)
    {
        int status;
        if (CellContainsReco(map,indices,end_point,status)) return;

        std::cout << "\n-------------------------------------------------------------------------------\n";
        std::cerr << "ERROR: Found pseudocell does not contain the endpoint in its circumscribed cube (status " << status << ").\n\n";

        std::cout << "Cube for (x1,y1,t1) = (" << end_point.x() << ", " << end_point.y() << ", " << end_point.t << "): \n";
        for (int i = 0; i < 8; i++)
        {
            const auto [ci,cj,ck] = CubeCorner({0, 1, 0, 1, 0, 1}, i);
            const auto [xi,yi,zi] = CubeCorner(indices, i);
            std::cout << ci << cj << ck << ": [" << xi << "][" << yi << "][" << zi << "],";
            std::cout << " (x,y,t) = (" << map.at(xi,yi,zi).point.x() << ", " << map.at(xi,yi,zi).point.y() << ", " << map.at(xi,yi,zi).point.t << ")\n";
        }
        std::cout << "-------------------------------------------------------------------------------\n\n";
    }

    RecoPoint Reconstruct(const Field<MapPoint>& map, EndPoint end_point, TGraph2D* g_map_pts)
    {
        const double x1 = end_point.x();
        const double y1 = end_point.y();
        const double t1 = end_point.t;

        // Find 8 closest points using binary search (assuming ordering).
        int boundaries[6] = { 0, map.GetXCells() - 1, 0, map.GetYCells() - 1, map.GetZCells() - 1, 0 };
        int& ximin   = boundaries[0]; // Starting minimal search x index.
        int& ximax   = boundaries[1]; // Starting maximal search x index.
        int& yimin   = boundaries[2]; // Starting minimal search y index.
        int& yimax   = boundaries[3]; // Starting maximal search y index.
        int& zi_tmin = boundaries[4]; // Starting minimal t search z index (inverted - time is maximal for z minimal).
        int& zi_tmax = boundaries[5]; // Starting maximal t search z index (inverted - time is maximal for z minimal).

        int cells[3] = { map.GetXCells(), map.GetYCells(), map.GetZCells() };

        int i_mid[3] = { (ximin + ximax) / 2, (yimin + yimax) / 2, (zi_tmin + zi_tmax) / 2}; // Midpoint array.

        int i = 0;
        int status;
        do
        {
            if (i > 3)
            {
                // std::cerr << "WARNING: Three iterations did not find the pseudocell.\n";
                break;
            }

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
            FindMapMinMaxIndex(map, zi_tmin, zi_tmax, i_mid, 2, t1);

            i++;
        }
        while (!CellContainsReco(map, boundaries, end_point, status));

        // Fix off-by-one-cell error
        // PrintCube(map, boundaries, end_point);
        bool fix_failed = false;
        while (status != -1)
        {
            // std::cout << "Fixing status " << status << "\n";

            int var = status / 2;
            bool err_max = status % 2;

            bool var_is_z = var == 2;
            if (var_is_z) err_max = !err_max;
            
            if (err_max && boundaries[2*var+(int)!var_is_z] < cells[var]-1)
            {
                // std::cout << "Boundaries increased\n";
                boundaries[2*var]++;
                boundaries[2*var+1]++;
            }
            else if (!err_max && boundaries[2*var+(int)var_is_z] > 0)
            {
                // std::cout << "Boundaries decreased\n";
                boundaries[2*var]--;
                boundaries[2*var+1]--;
            }
            else
            {
                fix_failed = true;
                // std::cerr << "ERROR: Fixing failed\n";
                break;
            }

            CellContainsReco(map, boundaries, end_point, status);
        }

        // Sanity check.
        if(!fix_failed) PrintCube(map, boundaries, end_point);

        // Plotting map points used for reconstruction.
        if (g_map_pts)
        for (int i = 0; i < 8; i++)
        {
            const auto [xi,yi,zi] = CubeCorner(boundaries,i);
            // X17::MapPoint map_pt = map.at(xi,yi,zi);

            double map_step = map.GetStep();
            g_map_pts->AddPoint(map.GetXMin() + xi*map_step, map.GetYMin() + yi*map_step, map.GetZMin() + zi*map_step);
        }

        // Interpolation from map.
        std::array<std::array<double,8>,3> coef = GetInterpolCoef(map, boundaries);

        double xout = coef[0][0] + coef[0][1]*x1 + coef[0][2]*y1 + coef[0][3]*t1 + coef[0][4]*x1*y1 + coef[0][5]*x1*t1 + coef[0][6]*y1*t1 + coef[0][7]*x1*y1*t1;
        double yout = coef[1][0] + coef[1][1]*x1 + coef[1][2]*y1 + coef[1][3]*t1 + coef[1][4]*x1*y1 + coef[1][5]*x1*t1 + coef[1][6]*y1*t1 + coef[1][7]*x1*y1*t1;
        double zout = coef[2][0] + coef[2][1]*x1 + coef[2][2]*y1 + coef[2][3]*t1 + coef[2][4]*x1*y1 + coef[2][5]*x1*t1 + coef[2][6]*y1*t1 + coef[2][7]*x1*y1*t1;

        double x,y,z;
        // std::cout << "Min corner:          (" << map.GetXMin()+ximin*map.GetStep() << ", " << map.GetYMin()+yimin*map.GetStep() << ", " << map.GetZMin()+zi_tmax*map.GetStep() << ")\n";
        // std::cout << "Reconstructed point: (" << xout << ", " << yout << ", " << zout << ")\n";
        return RecoPoint(xout, yout, zout, 1);
    }

    double Offset(MapPoint p, double x1, double y1, double t1, double tfact)
    {
        return std::sqrt(pow(x1 - p.point.x(), 2) + pow(y1 - p.point.y(), 2) + pow(tfact * (t1 - p.point.t), 2));
    }

    RecoPoint ReconstructOld(const Field<MapPoint>& map, double x1, double y1, double t1, double max_err, bool gas9010)
    {
        // Start looking at the same position.
        double x = x1;
        double y = y1;
        double z = (map.GetZMax() + map.GetZMin()) / 2;
        double step = map.GetStep() / 1000;

        double offset;      // Metric of distance between points.
        int iterations = 0; // Number of iterations should not exceed 100.
        double damp = 0.05; // Damping coefficient.

        double tfact = gas9010 ? 0.00327 : 0.000939;

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

            double oxa = Offset(xa,x1,y1,t1,tfact);
            double oxb = Offset(xb,x1,y1,t1,tfact);
            double oya = Offset(ya,x1,y1,t1,tfact);
            double oyb = Offset(yb,x1,y1,t1,tfact);
            double oza = Offset(za,x1,y1,t1,tfact);
            double ozb = Offset(zb,x1,y1,t1,tfact);

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
            offset = Offset(cur,x1,y1,t1,tfact);

            // Make sure step isn't too high.
            if (offset < 10 * step) step /= 10;

            double damped_grad_size = damp*std::sqrt(gradx*gradx + grady*grady + gradz*gradz);
            if (damped_grad_size > 0.1 * offset) damp /= 10; 

            iterations++;
            
            if (iterations == 1000) 
                std::cout << "1000 iterations.\n";
        }
        while ((offset > max_err) && (iterations < 1000));

        return RecoPoint(x,y,z,1);
    }
} // namespace X17