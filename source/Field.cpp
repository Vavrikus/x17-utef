// X17 dependencies
#include "Field.h"
#include "X17Utilities.h"

namespace X17
{
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
                    using namespace constants;

                    double x,y,z;    // Coordinates of the point in space [cm].
                    double vx,vy,vz; // Components of the vector in given point.
                    
                    // Magnetic and electric field data files contain coordinates in meters.
                    x = m2cm*stod(X); y = m2cm*stod(Y); z = m2cm*stod(Z);
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