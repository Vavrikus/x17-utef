// C++ dependencies
#include <iostream>

// ROOT dependencies
#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>

// Garfield++ dependencies
#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/ComponentGrid.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewSignal.hh"

// X17 dependencies
#include "Points.h"
#include "X17Utilities.h"

using namespace Garfield;
using namespace X17::constants;

int main(int argc, char *argv[])
{
    TApplication app("app", &argc, argv);

    // Set the gas mixture.
    MediumMagboltz gas;
    gas.SetComposition("ar", 90., "co2", 10.); std::cout << "BAD GAS COMPOSITION!!!!\n"; // Change to 70/30, set temperature, etc.

    // Set the output file.
    TFile outFile("../../../data/single_track/electrons.root","RECREATE","Electrons from ionization track");
    TTree electrons("electrons","Tree of initial and final points of electrons");

    X17::MicroPoint point;
    electrons.Branch("point",&point);

    // Add magnetic and electric field.
    ComponentGrid grid;
    const double m2cm = 100.;
    grid.LoadMagneticField("../../../data/elmag/VecB2.txt", "xyz", m2cm); 
    grid.LoadElectricField("../../../data/elmag/VecE2.txt", "xyz", false, false, m2cm);
    grid.SetMedium(&gas);

    // Assemble a sensor.
    Sensor sensor;
    sensor.AddComponent(&grid);
    constexpr double space = 1; // Extra space on the sensor so electrons always end on zmax.
    sensor.SetArea(-space, -space, zmin-space, xmax+space, space, zmax);

    // We use microscopic tracking for the electron avalanche simulation.
    AvalancheMicroscopic aval;
    aval.SetSensor(&sensor);
    aval.EnableSignalCalculation(); 
    aval.EnableMagneticField();

    // Simulate an ionizing particle (negative pion) using Heed.
    TrackHeed track;
    track.SetParticle("electron");    // Set the particle type.
    constexpr double momentum = 8.e6; // Set the particle momentum [eV / c].
    track.SetMomentum(momentum);
    track.SetSensor(&sensor);
    track.EnableMagneticField();
    track.EnableElectricField();
    track.DisableDeltaElectronTransport();  // This will disable secondary electrons in the track.
    track.EnablePhotonReabsorption(false);  // Enable/disable fluorescence reabsorption.

    // Construct a viewer to visualise the drift lines.
    ViewDrift driftView;
    track.EnablePlotting(&driftView);
    aval.EnablePlotting(&driftView);
    aval.EnableExcitationMarkers(false);
    aval.EnableIonisationMarkers(false);
    
    // Get the default parameters.
    double maxrange = 0., rforstraight = 0., stepstraight = 0., stepcurved = 0.;
    track.GetSteppingLimits(maxrange, rforstraight, stepstraight, stepcurved);

    // Reduce the step size [rad].
    stepcurved = 0.04;
    maxrange = 0.2;
    track.SetSteppingLimits(maxrange, rforstraight, stepstraight, stepcurved);

    // Set the starting point and momentum vector of the particle. Randomize in the future.
    double xt = 0.0; // [cm]
    double yt = 0.0; // [cm]
    double zt = 0.0; // [cm]
    double ti = 0;   // [ns]
    double px = 1;
    double py = 0;
    double pz = 0;
    int k = 0;

    // Now simulate a track, with p0 = x̂.
    track.NewTrack(xt, yt, zt, ti, px, py, pz);

    // Loop over the clusters.
    double xc, yc, zc, tc, ec, extra;
    int nc;
    while (track.GetCluster(xc, yc, zc, tc, nc, ec, extra)) 
    {
        std::cout << "Distance to origin: " << sqrt(xc*xc+yc*yc+zc*zc) << "  time " << tc << "  number " << k << "\n";
        for (int j = 0; j < nc; ++j) 
        {
            double xe, ye, ze, te, ee, dxe, dye, dze;
            track.GetElectron(j, xe, ye, ze, te, ee, dxe, dye, dze);
            // Simulate the drift/avalanche of this electron.
            aval.AvalancheElectron(xe, ye, ze, te, 0.1, dxe, dye, dze);
            // Move electrons that hit the mesh plane into the amplification gap.
            int status;
            aval.GetElectronEndpoint(0, point.start.point.x, point.start.point.y, point.start.point.z, point.start.t, point.e0, point.end.point.x, point.end.point.y, point.end.point.z, point.end.t, point.e1, status);
            electrons.Fill();
            k++;
        }
    }

    std::cout << "Total electron generated: "<< k << std::endl;

    // Plotting.
    // These canvases will be used to display the drift lines and the field.
    if (bool make_plots = true)
    {
        gStyle->SetPadRightMargin(0.15);
        TCanvas* c1 = new TCanvas("c1", "", 800, 800);
        driftView.SetPlane(1,0, 0, 0, 0, 0);
        driftView.SetArea(zmin, -space, zmax, space);
        driftView.SetCanvas(c1);
        constexpr bool twod = true;
        driftView.Plot(twod);
        c1->SaveAs("plot.pdf");

        TCanvas* c2 = new TCanvas("c2", "", 800, 800);
        driftView.SetPlane(0,1,0, 0, 0, 0);
        driftView.SetArea(0, zmin, xmax + space, zmax);
        driftView.SetCanvas(c2);
        constexpr bool twod2 = true;
        driftView.Plot(twod2);
        c2->SaveAs("plot2.pdf");  

        TCanvas* c3 = new TCanvas("c3", "", 800, 800);
        driftView.SetPlane(0,0,1, 0, 0, 0);
        driftView.SetArea(0, -space, xmax + space, space);
        driftView.SetCanvas(c3);
        constexpr bool twod3 = true;
        driftView.Plot(twod3);
        c3->SaveAs("plot3.pdf");
    }

    electrons.Print();
    outFile.Write();
    outFile.Close();

    app.Run(true);

    return 0;
}