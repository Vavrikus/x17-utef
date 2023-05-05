#pragma once

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "Vector.h"
#include "X17Utilities.h"

namespace X17
{
    /// @brief A rectangular field with a regular grid of values. The field is defined by its minimum
    ///        and maximum coordinates along each axis, the grid spacing, and a default value
    ///        for all points. The field is represented as a 1D array of values stored in a vector.
    /// @tparam T The type of the values stored in the field.
    template<typename T>
    class Field
    {
    private:
        double m_xmin;      // The minimum x coordinate of the field (meters).
        double m_ymin;      // The minimum y coordinate of the field (meters).
        double m_zmin;      // The minimum z coordinate of the field (meters).
        double m_xmax;      // The maximum x coordinate of the field (meters).
        double m_ymax;      // The maximum y coordinate of the field (meters).
        double m_zmax;      // The maximum z coordinate of the field (meters).
        double m_step_size; // The spacing between grid points (meters).
        T m_def_value;      // The default value for all points in the field.

        int m_num_x_cells; // The number of grid points along the x axis.
        int m_num_y_cells; // The number of grid points along the y axis.
        int m_num_z_cells; // The number of grid points along the z axis.

        std::vector<T> m_data; // The array of values representing the field. Filled in the order x,y,z.

        /// @brief Maps given 3D coordinates to the corresponding grid cell indices.
        /// @param x The x-coordinate of the point to be mapped (meters).
        /// @param y The y-coordinate of the point to be mapped (meters).
        /// @param z The z-coordinate of the point to be mapped (meters).
        /// @param xi Output parameter that will be set to the x-index of the corresponding grid cell.
        /// @param yi Output parameter that will be set to the y-index of the corresponding grid cell.
        /// @param zi Output parameter that will be set to the z-index of the corresponding grid cell.
        /// @throws std::out_of_range if the specified coordinates are out of bounds.
        void GetPointIndices(const double& x, const double& y, const double& z, int& xi, int& yi, int& zi) const
        {
            // Check if the point is out of bounds.
            if(x < m_xmin || x > m_xmax || y < m_ymin || y > m_ymax || z < m_zmin || z > m_zmax)
            {
                throw std::out_of_range("Cannot read field out of bounds.");
            }

            // Calculate the indices of the corresponding grid cell.
            xi = round((x - m_xmin) / m_step_size);
            yi = round((y - m_ymin) / m_step_size);
            zi = round((z - m_zmin) / m_step_size);
        }

    public:
        /// @brief Deleted copy constructor. 
        Field(const Field&) = delete;

        /// @brief Constructs a Field object with the given dimensions and step size, and sets the default value for all points.
        /// @param min_corner The minimum corner of the field as a Vector object (meters).
        /// @param max_corner The maximum corner of the field as a Vector object (meters).
        /// @param step The spacing between grid points.
        /// @param def_value The default value for all points in the field.
        Field(const Vector& min_corner, const Vector& max_corner, const double& step, const T& def_value)
            : m_xmin(min_corner.x), m_xmax(max_corner.x), m_ymin(min_corner.y), m_ymax(max_corner.y),
              m_zmin(min_corner.z), m_zmax(max_corner.z), m_step_size(step), m_def_value(def_value)
        {
            // Determine the number of grid points along each axis.
            this->m_num_x_cells = round((max_corner.x - min_corner.x) / step) + 1;
            this->m_num_y_cells = round((max_corner.y - min_corner.y) / step) + 1;
            this->m_num_z_cells = round((max_corner.z - min_corner.z) / step) + 1;

            // Resize the data array and initialize it with default value.
            this->m_data.resize(GetNCells(), def_value);
        }

        /// @brief Get the value of the data at a given index (xi,yi,zi).
        /// @param xi The x index of the data.
        /// @param yi The y index of the data.
        /// @param zi The z index of the data.
        /// @return The value of the data at the given index.
        T& at(const int& xi, const int& yi, const int& zi)
        {
            return m_data[zi * m_num_y_cells * m_num_x_cells + yi * m_num_x_cells + xi];
        }

        /// @brief Get the value of the data at a given index (xi,yi,zi).
        /// @param xi The x index of the data.
        /// @param yi The y index of the data.
        /// @param zi The z index of the data.
        /// @return The value of the data at the given index.
        const T& at(const int& xi, const int& yi, const int& zi) const
        {
            return m_data[zi * m_num_y_cells * m_num_x_cells + yi * m_num_x_cells + xi];
        }

        /// @brief Get the value of the data at a given index (xi,yi,zi).
        /// @param i_arr The (xi,yi,zi) index array of the data.
        /// @return The value of the data at the given index.
        const T& at(const int (&i_arr)[3]) const
        {
            return this->at(i_arr[0], i_arr[1], i_arr[2]);
        }

        /// @brief Get the minimum x coordinate of the field.
        /// @return double The minimum x coordinate.
        double GetXMin() const { return m_xmin; }

        /// @brief Get the maximum x coordinate of the field.
        /// @return double The maximum x coordinate.
        double GetXMax() const { return m_xmax; }

        /// @brief Get the minimum y coordinate of the field.
        /// @return double The minimum y coordinate.
        double GetYMin() const { return m_ymin; }

        /// @brief Get the maximum y coordinate of the field.
        /// @return double The maximum y coordinate.
        double GetYMax() const { return m_ymax; }

        /// @brief Get the minimum z coordinate of the field.
        /// @return double The minimum z coordinate.
        double GetZMin() const { return m_zmin; }

        /// @brief Get the maximum z coordinate of the field.
        /// @return double The maximum z coordinate.
        double GetZMax() const { return m_zmax; }

        /// @brief Get the step size of the field.
        /// @return double The step size.
        double GetStep() const { return m_step_size; }

        /// @brief Get the number of x cells.
        /// @return double The number of x cells.
        double GetXCells() const { return m_num_x_cells; }

        /// @brief Get the number of y cells.
        /// @return double The number of y cells.
        double GetYCells() const { return m_num_y_cells; }

        /// @brief Get the number of z cells.
        /// @return double The number of z cells.
        double GetZCells() const { return m_num_z_cells; }

        /// @brief Returns the number of grid cells.
        /// @return Total number of grid cells.
        int GetNCells() const
        {
            return this->m_num_x_cells * this->m_num_y_cells * this->m_num_z_cells;
        }

        /// @brief Returns a pointer to the field value at the specified coordinates.
        /// @param x The x-coordinate of the point to retrieve (meters).
        /// @param y The y-coordinate of the point to retrieve (meters).
        /// @param z The z-coordinate of the point to retrieve (meters).
        /// @return A pointer to the field value at the specified coordinates.
        /// @throws std::out_of_range if the specified coordinates are out of bounds.
        T* GetPoint(const double& x, const double& y, const double& z)
        {
            int xi, yi, zi;
            GetPointIndices(x, y, z, xi, yi, zi);

            return &this->at(xi,yi,zi);
        }

        /// @brief Sets the value of the field at the specified point.
        /// @param x The x coordinate of the point (meters).
        /// @param y The y coordinate of the point (meters).
        /// @param z The z coordinate of the point (meters).
        /// @param new_value The new value to set at the point.
        /// @throws std::out_of_range if the specified coordinates are out of bounds.
        void SetPoint(const double& x, const double& y, const double& z, const T& new_value)
        {
            *(this->GetPoint(x,y,z)) = new_value;
        }

        /// @brief Sets the value of the field at the specified point.
        /// @param position The vector with position of the point (meters).
        /// @param new_value The new value to set at the point.
        /// @throws std::out_of_range if the specified coordinates are out of bounds.
        void SetPoint(const Vector& position, const T& new_value)
        {
            *(this->GetPoint(position.x,position.y,position.z)) = new_value;
        }

        /// @brief Get the field value at the specified location using trilinear interpolation.
        /// @param x The x-coordinate of the location to sample (meters).
        /// @param y The y-coordinate of the location to sample (meters).
        /// @param z The z-coordinate of the location to sample (meters).
        /// @return The interpolated field value at the specified location.
        T GetField(double x, double y, double z) const
        {
            int xi, yi, zi;
            GetPointIndices(x, y, z, xi, yi, zi);

            // Determine the grid coordinates of given point.
            int x_grid = m_xmin + m_step_size * xi;
            int y_grid = m_ymin + m_step_size * yi;
            int z_grid = m_zmin + m_step_size * zi;

            // Determine the indices of the eight surrounding points.
            int xi2, yi2, zi2;
            if(x - x_grid < 0) xi2 = xi-1; else xi2 = xi + 1;
            if(y - y_grid < 0) yi2 = yi-1; else yi2 = yi + 1;
            if(z - z_grid < 0) zi2 = zi-1; else zi2 = zi + 1;

            // Make sure the indices are in bounds.
            std::clamp(xi2, 0, m_num_x_cells - 1);
            std::clamp(yi2, 0, m_num_y_cells - 1);
            std::clamp(zi2, 0, m_num_z_cells - 1);

            // Compute the trilinear interpolation coefficients.
            double dx = abs((x - x_grid) / m_step_size);
            double dy = abs((y - y_grid) / m_step_size);
            double dz = abs((z - z_grid) / m_step_size);

            // Get references to the eight surrounding points.
            T& c000 = this->at(xi,  yi,  zi);
            T& c001 = this->at(xi,  yi,  zi2);
            T& c010 = this->at(xi,  yi2, zi);
            T& c011 = this->at(xi,  yi2, zi2);
            T& c100 = this->at(xi2, yi,  zi);
            T& c101 = this->at(xi2, yi,  zi2);
            T& c110 = this->at(xi2, yi2, zi);
            T& c111 = this->at(xi2, yi2, zi2);

            // Compute the interpolated field value.
            T c00 = (1 - dx) * c000 + dx * c100;
            T c01 = (1 - dx) * c001 + dx * c101;
            T c10 = (1 - dx) * c010 + dx * c110;
            T c11 = (1 - dx) * c011 + dx * c111;

            T c0 = (1 - dy) * c00 + dy * c10;
            T c1 = (1 - dy) * c01 + dy * c11;

            return (1 - dz) * c0 + dz * c1;
        }

        /// @brief Returns the trilinearly interpolated value of the field at the specified position.
        /// @param vec The Vector object representing the position to interpolate the field at (meters).
        /// @return The interpolated value of the field at the specified position.
        T GetField(const Vector& vec) const
        {
            return GetField(vec.x,vec.y,vec.z);
        }
    };

    /// @brief Loads field data from a file and stores it in a Field<Vector>.
    ///        The function reads a file containing field data in the format X Y Z VX VY VZ, where X, Y, Z
    ///        represent the spatial coordinates and VX, VY, VZ represent the vector components of the field.
    ///        The function stores the field data in a Field<Vector> object.
    /// @param filename The name of the file to load the field data from.
    /// @param min_corner The minimum corner of the field as a Vector object (meters).
    /// @param max_corner The maximum corner of the field as a Vector object (meters).
    /// @param step The spacing between grid points.
    /// @param printInfo If true, prints information about maximal and minimal field magnitude and angle values.
    Field<Vector>* LoadField(const char* filename, const Vector& min_corner, const Vector& max_corner, const double& step, bool printInfo = false)
    {
        Field<Vector>* field = new Field<Vector>(min_corner,max_corner,step,{0,0,0});

        std::ifstream inf {filename};
        std::cout << "\nLoading field from " << filename << "\n";

        int lines_read = 0;
        int lines_processed = 0;
        int lines_expected = field->GetNCells();

        while (inf)
        {
            std::string X,Y,Z,VX,VY,VZ;
            inf >> X; inf >> Y; inf >> Z; inf >> VX; inf >> VY; inf >> VZ;
            lines_read++;

            try
            {                
                if (X != "")
                {
                    double x,y,z,vx,vy,vz;
                    x = stod(X); y = stod(Y); z = stod(Z);
                    vx = stod(VX); vy = stod(VY); vz = stod(VZ);

                    *(field->GetPoint(x,y,z)) = Vector{vx,vy,vz};
                    lines_processed++;
                }
            }

            catch (const std::exception& e)
            {
                std::cerr << e.what() << '\n';
                std::cout << "X: " << X << " Y: " << Y << " Z: " << Z << " VX: " << VX << " VY: " << VY << " VZ: " << VZ << "\n";
            }
        }

        std::cout << "Lines read: " << lines_read << " processed: " << lines_processed << " expected: " << lines_expected << "\n\n";

        if (printInfo)
        {
            double minfield,maxfield,minangle,maxangle;
            GetMinMaxField(*field,minfield,maxfield,0);
            GetMinMaxFieldAngle(*field,minangle,maxangle,0);
            std::cout << "At least 0.0 cm from TPC walls: minimal magnetic field: " << minfield << " maximal: " << maxfield;
            std::cout << " minimal angle to electric field: " << minangle << " maximal: " << maxangle << "\n";
            GetMinMaxField(*field,minfield,maxfield,0.5);
            GetMinMaxFieldAngle(*field,minangle,maxangle,0.5);
            std::cout << "At least 0.5 cm from TPC walls: Minimal magnetic field: " << minfield << " maximal: " << maxfield;
            std::cout << " minimal angle to electric field: " << minangle << " maximal: " << maxangle << "\n";
            GetMinMaxField(*field,minfield,maxfield,1);
            GetMinMaxFieldAngle(*field,minangle,maxangle,1);
            std::cout << "At least 1.0 cm from TPC walls: Minimal magnetic field: " << minfield << " maximal: " << maxfield;
            std::cout << " minimal angle to electric field: " << minangle << " maximal: " << maxangle << "\n";
        }

        return field;
    }
} // namespace X17