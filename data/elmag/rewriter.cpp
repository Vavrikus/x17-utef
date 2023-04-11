#include <fstream>
#include <iostream>

void rewriter()
{
    std::ifstream inf{"/home/vavrik/work/X17/mag_data/VecE.txt"};
    std::ofstream outf{"/home/vavrik/work/X17/mag_data/VecE2.txt"};

    while (inf)
    {
        std::string X,Y,Z,VX,VY,VZ;
        inf >> X; inf >> Y; inf >> Z; inf >> VX; inf >> VY; inf >> VZ;
        try
        {                
            if (X != "")
            {
                // changing coordinate system from magnetic simulation (x,y,z) --> (y,z,x)
                outf << Z << " " << X << " " << Y << " " << VZ << " " << VX << " " << VY << "\n";
            }
            else std::cout << "Line skipped.\n";
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
            std::cout << "X: " << X << " Y: " << Y << " Z: " << Z << " VX: " << VX << " VY: " << VY << " VZ: " << VZ << "\n";
        }
    }

}

int main()
{
    rewriter();
    return 0;
}