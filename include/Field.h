#pragma once

// C++ dependencies
#include <vector>

// X17 dependencies
#include "Vector.h"

namespace X17
{
    /// @brief Rectangular field with a regular grid of values. The field is defined by its minimum
    ///        and maximum coordinates along each axis, the grid spacing, and a default value
    ///        for all points. The field is represented as a 1D array of values stored in a vector.
    /// @tparam T The type of the values stored in the field.
    template<typename T>
    class Field
    {
    private:
        double m_xmin;      // Minimum x coordinate of the field [cm].
        double m_ymin;      // Minimum y coordinate of the field [cm].
        double m_zmin;      // Minimum z coordinate of the field [cm].
        double m_xmax;      // Maximum x coordinate of the field [cm].
        double m_ymax;      // Maximum y coordinate of the field [cm].
        double m_zmax;      // Maximum z coordinate of the field [cm].
        double m_step_size; // Spacing between grid points [cm].
        T m_def_value;      // Default value for all points in the field.

        int m_num_x_cells; // Number of grid points along the x axis.
        int m_num_y_cells; // Number of grid points along the y axis.
        int m_num_z_cells; // Number of grid points along the z axis.

        std::vector<T> m_data; // Array of values representing the field. Filled in the order x,y,z.

    public:
        /// @brief Deleted copy constructor. 
        Field(const Field&) = delete;

        /// @brief Constructs a Field object with the given dimensions and step size, and sets the default value for all points.
        /// @param min_corner The minimum corner of the field as a Vector object [cm].
        /// @param max_corner The maximum corner of the field as a Vector object [cm].
        /// @param step The spacing between grid points.
        /// @param def_value The default value for all points in the field.
        Field(Vector min_corner, Vector max_corner, double step, const T& def_value);

        /// @brief Get the value of the data at a given index (xi,yi,zi).
        /// @param xi The x index of the data.
        /// @param yi The y index of the data.
        /// @param zi The z index of the data.
        /// @return Value of the data at the given index.
        T& at(int xi, int yi, int zi)
        {
            return m_data[zi * m_num_y_cells * m_num_x_cells + yi * m_num_x_cells + xi];
        }

        /// @brief Get the value of the data at a given index (xi,yi,zi).
        /// @param xi The x index of the data.
        /// @param yi The y index of the data.
        /// @param zi The z index of the data.
        /// @return Value of the data at the given index.
        const T& at(int xi, int yi, int zi) const
        {
            return m_data[zi * m_num_y_cells * m_num_x_cells + yi * m_num_x_cells + xi];
        }

        /// @brief Get the value of the data at a given index (xi,yi,zi).
        /// @param i_arr The (xi,yi,zi) index array of the data.
        /// @return Value of the data at the given index.
        const T& at(const int (&i_arr)[3]) const
        {
            return this->at(i_arr[0], i_arr[1], i_arr[2]);
        }

        /// @brief Get the minimum x coordinate of the field.
        /// @return Minimum x coordinate.
        double GetXMin() const { return m_xmin; }

        /// @brief Get the maximum x coordinate of the field.
        /// @return Maximum x coordinate.
        double GetXMax() const { return m_xmax; }

        /// @brief Get the minimum y coordinate of the field.
        /// @return Minimum y coordinate.
        double GetYMin() const { return m_ymin; }

        /// @brief Get the maximum y coordinate of the field.
        /// @return Maximum y coordinate.
        double GetYMax() const { return m_ymax; }

        /// @brief Get the minimum z coordinate of the field.
        /// @return Minimum z coordinate.
        double GetZMin() const { return m_zmin; }

        /// @brief Get the maximum z coordinate of the field.
        /// @return Maximum z coordinate.
        double GetZMax() const { return m_zmax; }

        /// @brief Get the step size of the field.
        /// @return Step size.
        double GetStep() const { return m_step_size; }

        /// @brief Get the number of x cells.
        /// @return Number of x cells.
        double GetXCells() const { return m_num_x_cells; }

        /// @brief Get the number of y cells.
        /// @return Number of y cells.
        double GetYCells() const { return m_num_y_cells; }

        /// @brief Get the number of z cells.
        /// @return Number of z cells.
        double GetZCells() const { return m_num_z_cells; }

        /// @brief Returns the number of grid cells.
        /// @return Total number of grid cells.
        int GetNCells() const
        {
            return this->m_num_x_cells * this->m_num_y_cells * this->m_num_z_cells;
        }

        /// @brief Returns a pointer to the field value at the specified coordinates.
        /// @param x The x-coordinate of the point to retrieve [cm].
        /// @param y The y-coordinate of the point to retrieve [cm].
        /// @param z The z-coordinate of the point to retrieve [cm].
        /// @return Pointer to the field value at the specified coordinates.
        /// @throws std::out_of_range if the specified coordinates are out of bounds.
        T* GetPoint(double x, double y, double z)
        {
            int xi, yi, zi;
            _GetPointIndices(x, y, z, xi, yi, zi);

            return &this->at(xi,yi,zi);
        }

        /// @brief Sets the value of the field at the specified point.
        /// @param x The x coordinate of the point [cm].
        /// @param y The y coordinate of the point [cm].
        /// @param z The z coordinate of the point [cm].
        /// @param new_value New value to set at the point.
        /// @throws std::out_of_range if the specified coordinates are out of bounds.
        void SetPoint(double& x, double& y, double& z, const T& new_value)
        {
            *(this->GetPoint(x,y,z)) = new_value;
        }

        /// @brief Sets the value of the field at the specified point.
        /// @param position Vector with position of the point [cm].
        /// @param new_value New value to set at the point.
        /// @throws std::out_of_range if the specified coordinates are out of bounds.
        void SetPoint(Vector position, const T& new_value)
        {
            *(this->GetPoint(position.x,position.y,position.z)) = new_value;
        }

        /// @brief Get the field value at the specified location using trilinear interpolation.
        /// @param x The x-coordinate of the location to sample [cm].
        /// @param y The y-coordinate of the location to sample [cm].
        /// @param z The z-coordinate of the location to sample [cm].
        /// @return Interpolated field value at the specified location.
        T GetField(double x, double y, double z) const;

        /// @brief Returns the trilinearly interpolated value of the field at the specified position.
        /// @param vec Vector object representing the position to interpolate the field at [cm].
        /// @return Interpolated value of the field at the specified position.
        T GetField(Vector vec) const
        {
            return GetField(vec.x,vec.y,vec.z);
        }

    private:
        /// @brief Maps given 3D coordinates to the corresponding grid cell indices.
        /// @param x The x-coordinate of the point to be mapped [cm].
        /// @param y The y-coordinate of the point to be mapped [cm].
        /// @param z The z-coordinate of the point to be mapped [cm].
        /// @param xi Output parameter that will be set to the x-index of the corresponding grid cell.
        /// @param yi Output parameter that will be set to the y-index of the corresponding grid cell.
        /// @param zi Output parameter that will be set to the z-index of the corresponding grid cell.
        /// @throws std::out_of_range if the specified coordinates are out of bounds.
        void _GetPointIndices(double x, double y, double z, int& xi, int& yi, int& zi) const;
    };

    /// @brief Loads field data from a file and stores it in a Field<Vector>.
    ///        The function reads a file containing field data in the format X Y Z VX VY VZ, where X, Y, Z
    ///        represent the spatial coordinates and VX, VY, VZ represent the vector components of the field.
    ///        The function stores the field data in a Field<Vector> object.
    /// @param filename Name of the file to load the field data from.
    /// @param min_corner Minimum corner of the field as a Vector object [cm].
    /// @param max_corner Maximum corner of the field as a Vector object [cm].
    /// @param step Spacing between grid points.
    /// @param printInfo If true, prints information about maximal and minimal field magnitude and angle values.
    Field<Vector>* LoadField(const char* filename, Vector min_corner, Vector max_corner, double step, bool printInfo = false);
} // namespace X17

// Templated function definitions.
#include "Field.inl"