// C++ dependencies
#include <iostream>
#include <cmath>

// ROOT dependencies
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <TGraph2D.h>
#include <TTree.h>

// Garfield++ dependencies
#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/ComponentGrid.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewSignal.hh"

// X17 dependencies
#include "MapJob.h"
#include "Points.h"

using namespace Garfield;

// Arguments from console: max_id, id, step (cm).
int main(int argc, char *argv[])
{
    TApplication app("app", &argc, argv);

    X17::MapJob job;
    job.SetParameters(argc,argv);
    job.SetElectronBounds();

    // Set the output file, name dependend on the job's id.
    std::string filename = "ion" + std::to_string(job.id) + ".root";

    TFile outFile(filename.c_str(),"RECREATE","Electrons from ionization track");
    TTree electrons("electrons","Tree of initial and final points of secondary electrons");

    // Setting up the output TTree, contains initial and final time, position and energy of electrons.
    X17::MicroPoint point;
    point.MakeTTreeBranches(&electrons);

    // Set the gas mixture.
    MediumMagboltz gas;
    gas.SetComposition("ar", 70., "co2", 30.); // std::cout << "BAD GAS COMPOSITION!!!!\n"; // Change to 70/30, set temperature, etc.

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
    sensor.SetArea(job.xmin-space, job.ymin-space, job.zmin-space, job.xmax+space, job.ymax+space, job.zmax);

    // We use microscopic tracking for the electron avalanche simulation.
    AvalancheMicroscopic aval;
    aval.SetSensor(&sensor); 
    // Switch on signal calculation. 
    aval.EnableSignalCalculation(); 
    aval.EnableMagneticField();

    // The actual simulation.
    int i_el = 0; // Current index of the electron.
    for (double z = job.zmin; z <= job.zmax; z += job.step)
    for (double y = job.ymin; y <= job.ymax; y += job.step)
    for (double x = job.xmin; x <= job.xmax; x += job.step)
    {
        // Only inside of the first sector.
        if(!((job.SectorLineDist(x,y,false) <= 0) && (job.SectorLineDist(x,y,true) >= 0))) continue;
        i_el++;

        // Check if this electron is supposed to be simulated by job with this id.
        if (((i_el < job.min_el)||(i_el > job.max_el))) continue;

        // Actual for loop for electron generation.
        for (int j = 0; j < job.iterations; j++)
        {                        
            int status;
            aval.AvalancheElectron(x, y, z, 0, 0.1, 0, 0, 0);
            aval.GetElectronEndpoint(0,point.x0,point.y0,point.z0,point.t0,point.e0,point.x1,point.y1,point.z1,point.t1,point.e1,status);
            electrons.Fill();
            // std::cout << "X: " << x << " Y: " << y << " Z: " << z << "\n";
        }
    }

    outFile.Write();
    outFile.Close();

    return 0;
}