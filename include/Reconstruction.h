#pragma once

// C++ dependencies
#include <array>
#include <iostream>
#include <vector>

// X17 dependencies
#include "Field.h"
#include "Matrix.h"
#include "Points.h"

namespace X17
{
    /// @brief Templated function for cube corners selection.
    /// @tparam T Any number type.
    /// @param indices Array {xmin, xmax, ymin, ymax, zmin, zmax} in this order.
    /// @param corner  Index of the corner (0-7 <----> xyz: 000 --> 111).
    /// @return Array with x,y,z of corner in this order.
    template<typename T>
    std::array<T,3> CubeCorner(const T (&indices)[6], int corner);

    /// @brief Calculates the interpolation coefficients for a polynomial in three dimensions,
    ///        using the values of the MapPoint member variables in the given indices of the Field.
    /// @param map A reference to the Field of MapPoint values.
    /// @param indices An array of 6 integers representing the indices of the MapPoint values to use for interpolation.
    ///                The array should contain indices for {xmin, xmax, ymin, ymax, zmin, zmax}.
    /// @return A 2D vector containing the interpolation coefficients for each of the x, y, and t coordinates of the MapPoint values.
    ///         The first index of the vector selects the x (0), y (1), z (2), or t(3) coordinate,
    ///         while the second index selects the polynomial coefficient.
    std::vector<std::vector<double>> GetInterpolCoef(const Field<MapPoint>& map, const int (&indices)[6]);

    /// @brief Finds the indices of the elements closest to a given value in a 3D map along a given variable.
    /// @param map The 3D map to search through.
    /// @param i_min A reference to an integer that will hold the index of the element with lower value.
    /// @param i_max A reference to an integer that will hold the index of the element with higher value.
    /// @param i_mid An integer array that will hold the indices of the midpoint elements along each variable.
    /// @param var The variable to search along -- x (0), y (1), or z (2).
    /// @param value The value to search for.
    void FindMapMinMaxIndex(const Field<MapPoint>& map, int& i_min, int& i_max, int (&i_mid)[3], int var, double value);

    /// @brief Prints the cube defined by a set of indices and checks if an endpoint coordinates are within the bounds.
    /// @param map The 3D map of ionization electron drift.
    /// @param indices An array of 6 integers representing the indices of the MapPoint values to use for interpolation.
    ///                The array should contain indices for {xmin, xmax, ymin, ymax, zmin, zmax}.
    /// @param end_point The electron endpoint to be checked.
    /// @return void.
    void PrintCube(const Field<MapPoint>& map, const int (&indices)[6], const EndPoint& end_point);

    /// @brief Reconstructs an EndPoint using interpolation from the ionization electron drift map.
    /// @param map The 3D map of ionization electron drift used for interpolation.
    /// @param end_point An EndPoint object that represents the point in 3D space to be reconstructed. (x1,y1,t1)
    /// @return A RecoPoint object that represents the reconstructed point in 3D space. (x0,y0,z0)
    RecoPoint Reconstruct(const Field<MapPoint>& map, const EndPoint& end_point);

    /// @brief Reconstructs an EndPoint using interpolation from the ionization electron drift map.
    /// @param map The 3D map of ionization electron drift used for interpolation.
    /// @param point An MicroPoint object that represents the simulated point in 3D space to be reconstructed. (x1,y1,t1)
    /// @return A RecoPoint object that represents the reconstructed point in 3D space. (x0,y0,z0)
    RecoPoint Reconstruct(const Field<MapPoint>& map, const MicroPoint& point) { return Reconstruct(map, point.end); }

    /// @brief Estimates the distance between two points, with given coordinates and time values in a SensorData object. The time value is weighted by a factor.
    /// @param p The MapPoint from which we want to calculate distance.
    /// @param x1 The final x-coordinate.
    /// @param y1 The final y-coordinate.
    /// @param t1 The final time.
    /// @return The estimated distance between points.
    double Offset(const MapPoint& p, const double& x1, const double& y1, const double& t1);

    /// @brief Reconstructs a point with given final coordinates using the map of ion electron drift. Uses gradient descend. Old function.
    /// @param map The field map to use for reconstruction.
    /// @param x1 The x-coordinate of the initial point.
    /// @param y1 The y-coordinate of the initial point.
    /// @param t1 The time of the initial point.
    /// @param max_err The maximum allowed error between reconstructed point and actual point.
    /// @return The reconstructed point.
    RecoPoint ReconstructOld(const Field<MapPoint>& map, const double& x1, const double& y1, const double& t1, const double& max_err);
} // namespace X17