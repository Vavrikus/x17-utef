// C++ dependencies
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

// ROOT dependencies
#include "TApplication.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

// Garfield++ dependencies
#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentGrid.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewSignal.hh"

// X17 dependencies
#include "MapJob.h"
#include "Points.h"
#include "Utilities.h"
#include "Vector.h"

using namespace Garfield;

// Arguments from console: max_id, id, step (cm).
int main(int argc, char* argv[])
{
  TApplication app("app", &argc, argv);

  X17::MapJob job;
  job.SetParameters(argc, argv);
  job.SetElectronBounds();

  // Set the output file, name dependend on the job's id.
  std::string filename = "ion" + std::to_string(job.id) + ".root";

  TFile outFile(filename.c_str(), "RECREATE", "Electrons from ionization track");
  TTree electrons("electrons", "Tree of initial and final points of secondary electrons");

  // Setting up the output TTree, contains initial and final time, position and energy of electrons.
  X17::MicroPoint point;
  point.MakeTTreeBranches(&electrons);

  // Set the gas mixture.
  MediumMagboltz gas;
  gas.SetComposition("ar", 70., "co2",
                     30.); // std::cout << "BAD GAS COMPOSITION!!!!\n"; // Change to 70/30, set temperature, etc.

  // Add magnetic and electric field.
  ComponentGrid grid;
  const double m2cm = 100.;
  grid.LoadMagneticField("../../../data/elmag/VecB2.txt", "xyz", m2cm);
  grid.LoadElectricField("../../../data/elmag/VecE2.txt", "xyz", false, false, m2cm);
  grid.SetMedium(&gas);

  // Assemble a sensor.
  Sensor sensor;
  sensor.AddComponent(&grid);
  constexpr double space = 3; // Extra space on the sensor so electrons always end on zmax.
  sensor.SetArea(job.xmin - space, job.ymin - space, job.zmin - space, job.xmax + space, job.ymax + space, job.zmax);

  // We use microscopic tracking for the electron avalanche simulation.
  AvalancheMicroscopic aval;
  aval.SetSensor(&sensor);
  // Switch on signal calculation.
  aval.EnableSignalCalculation();
  aval.EnableMagneticField();

  // The actual simulation.
  int i_el = 0; // Current index of the electron.

  // Pre-calculate the number of steps for each axis
  const int nX = iround((job.xmax - job.xmin) / job.step);
  const int nY = iround((job.ymax - job.ymin) / job.step);
  const int nZ = iround((job.zmax - job.zmin) / job.step);

  // Loop using integers
  for (int iz = 0; iz <= nZ; ++iz)
  {
    double z = job.zmin + iz * job.step;

    for (int iy = 0; iy <= nY; ++iy)
    {
      double y = job.ymin + iy * job.step;

      for (int ix = 0; ix <= nX; ++ix)
      {
        double x = job.xmin + ix * job.step;

        // Only inside of the first sector.
        if ((job.SectorLineDist(x, y, false) > 0) || (job.SectorLineDist(x, y, true) < 0))
          continue;
        i_el++;

        // Check if this electron is supposed to be simulated by this job.
        if ((i_el < job.min_el) || (i_el > job.max_el))
          continue;

        std::cout << "Progress: " << i_el - job.min_el + 1 << "/" << job.max_el - job.min_el + 1;
        std::cout << "   z: " << std::setw(4) << z << " cm   y: " << std::setw(4) << y << " cm   x: " << std::setw(4)
                  << x << " cm\n";

        // Actual for loop for electron generation.
        for (int j = 0; j < job.iterations; j++)
        {
          int status;
          aval.AvalancheElectron(x, y, z, 0, 0.1, 0, 0, 0);
          aval.GetElectronEndpoint(0, point.start.point.x, point.start.point.y, point.start.point.z, point.start.t,
                                   point.e0, point.end.point.x, point.end.point.y, point.end.point.z, point.end.t,
                                   point.e1, status);
          electrons.Fill();
        }
      }
    }
  }

  outFile.Write();
  outFile.Close();

  return 0;
}