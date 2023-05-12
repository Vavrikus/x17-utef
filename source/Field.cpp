// X17 dependencies
#include "Field.h"
#include "X17Utilities.h"

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

    Field<Vector>* LoadField(const char* filename, const Vector& min_corner, const Vector& max_corner, const double& step, bool printInfo)
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