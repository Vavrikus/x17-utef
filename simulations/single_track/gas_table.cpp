#include <iostream>

#include "Garfield/MediumMagboltz.hh"

using namespace Garfield;
using namespace std;

int main(int argc, char * argv[])
{
    // Default pressure 760 torr = 1 atm
    // Default temperature 293.15 K
    // Default nColl = 10 (E+7)

    MediumMagboltz gas;
    gas.SetComposition("ar", 70., "co2", 30.);
    gas.SetFieldGrid(400,400,1,false,0,1,2,0,Pi,3);
    // gas.EnableThermalMotion();
    gas.GenerateGasTable();
    gas.WriteGasFile("Ar_70_CO2_30_M2A3_noTM.gas");

    return 0;
}