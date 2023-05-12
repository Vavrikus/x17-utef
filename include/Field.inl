// X17 dependencies
#include "Field.h"

namespace X17
{
    template <typename T>
    void Field<T>::GetPointIndices(const double& x, const double& y, const double& z, int& xi, int& yi, int& zi) const
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

    template <typename T>
    Field<T>::Field(const Vector& min_corner, const Vector& max_corner, const double& step, const T& def_value)
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

    template <typename T>
    T Field<T>::GetField(double x, double y, double z) const
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
        const T& c000 = this->at(xi,  yi,  zi);
        const T& c001 = this->at(xi,  yi,  zi2);
        const T& c010 = this->at(xi,  yi2, zi);
        const T& c011 = this->at(xi,  yi2, zi2);
        const T& c100 = this->at(xi2, yi,  zi);
        const T& c101 = this->at(xi2, yi,  zi2);
        const T& c110 = this->at(xi2, yi2, zi);
        const T& c111 = this->at(xi2, yi2, zi2);

        // Compute the interpolated field value.
        T c00 = (1 - dx) * c000 + dx * c100;
        T c01 = (1 - dx) * c001 + dx * c101;
        T c10 = (1 - dx) * c010 + dx * c110;
        T c11 = (1 - dx) * c011 + dx * c111;

        T c0 = (1 - dy) * c00 + dy * c10;
        T c1 = (1 - dy) * c01 + dy * c11;

        return (1 - dz) * c0 + dz * c1;
    }
} // namespace X17
