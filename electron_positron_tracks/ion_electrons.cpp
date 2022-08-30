#include <iostream>
#include <cmath>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <TGraph2D.h>
#include <TTree.h>

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

using namespace Garfield;
using namespace std;

int main(int argc, char * argv[])
{
    TApplication app("app", &argc, argv);

    // Set the gas mixture.
    MediumMagboltz gas;
    gas.SetComposition("ar", 90., "co2", 10.);

    // ofstream myfile;
    // myfile.open ("electrons.txt"); //name of the file to save data
    TFile outFile("ion.root","RECREATE","Electrons from ionization track");
    TTree electrons("electrons","Tree of initial and final points of secondary electrons");

    double x0, y0, z0, t0, e0;
    double x1, y1, z1, t1, e1;
    electrons.Branch("x0",&x0);
    electrons.Branch("y0",&y0);
    electrons.Branch("z0",&z0);
    electrons.Branch("t0",&t0);
    electrons.Branch("e0",&e0);
    electrons.Branch("x1",&x1);
    electrons.Branch("y1",&y1);
    electrons.Branch("z1",&z1);
    electrons.Branch("t1",&t1);
    electrons.Branch("e1",&e1);


    ComponentGrid grid;
    const double m2cm = 100.;
    grid.LoadMagneticField("VecB.txt", "xyz", m2cm); 
    grid.LoadElectricField("VecE.txt", "xyz",false,false, m2cm);
    grid.SetMedium(&gas);


    // 8 cm drift gap.
    constexpr double yDrift = 8;


    // Assemble a sensor.
    Sensor sensor;
    sensor.AddComponent(&grid); 
    sensor.SetArea(-25, -yDrift, -15, 25, yDrift, 15.);


    //~ // We use microscopic tracking for simulating the electron avalanche.
    // AvalancheMC aval;
    AvalancheMicroscopic aval;
    aval.SetSensor(&sensor); 
    // Switch on signal calculation. 
    aval.EnableSignalCalculation(); 
    aval.EnableMagneticField();

    //~ // Construct a viewer to visualise the drift lines.
    ViewDrift driftView;
    aval.EnablePlotting(&driftView);
    aval.EnableExcitationMarkers();
    aval.EnableIonisationMarkers();

    
    double xmin,xmax,ymin,ymax,zmin,zmax;
    xmin = -100; xmax = 100; ymin = -7; ymax = 8; zmin = 0; zmax = 15;
    double step = 3;

    for (double x = xmin; x <= xmax; x+=step)
    {
        for (double y = ymin; y <= ymax; y+=step)
        {
            for (double z = zmin; z <= zmax; z+=step)
            {
                //only inside 1 sector
                if((x-sqrt(3)*z<=0)&&(x+sqrt(3)*z>0))
                {
                    cout << "xyz: " << x << " " << y << " " << z << "\n";

                    int status;
                    // aval.AvalancheElectron(x, y, z, 0);
                    // aval.GetElectronEndpoint(0,x0,y0,z0,t0,x1,y1,z1,t1,status);
                    aval.AvalancheElectron(x, y, z, 0, 0.1, 0, 0, 0);
                    aval.GetElectronEndpoint(0,x0,y0,z0,t0,e0,x1,y1,z1,t1,e1,status);
                    electrons.Fill();
                }
            }
        }
    }

    gStyle->SetPadRightMargin(0.15);
    TCanvas* c2 = new TCanvas("c2", "", 800, 800);
    driftView.SetPlane(1,0, 0, 0, 0, 0);
    driftView.SetArea(-15, -8, 15, 8);
    driftView.SetCanvas(c2);
    constexpr bool twod = true;
    driftView.Plot(twod);

    electrons.Print();
    outFile.Write();
    outFile.Close();

    TCanvas* c1 = new TCanvas("c1", "", 800, 800);
    electrons.Draw("z0:t1","","ap");
    app.Run(true);

    return 0;
}
