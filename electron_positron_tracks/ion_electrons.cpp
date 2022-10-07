//C++ dependencies
#include <iostream>
#include <cmath>

//ROOT dependencies
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <TGraph2D.h>
#include <TTree.h>

//Garfield++ dependencies
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

//arguments: max_id, id, step (cm)
int main(int argc, char * argv[])
{
    TApplication app("app", &argc, argv);

    int max_id,id;
    double step;

    //check number of paramaters passed to main function
    if (argc < 4)
    {
        cerr << "ERROR: Missing arguments in ion_electrons. Correct arguments: max_id, id, step (cm).\n";
        return 1;
    }

    //set parameters of program
    max_id = stoi(argv[1]);
    id     = stoi(argv[2]);
    step   = stod(argv[3]);

    if (id <= 0) {cerr << "ERROR: Parameter id has to be positive.\n";return 2;}
    if (id > max_id) {cerr << "ERROR: Parameter id cannot be bigger than max_id.\n"; return 3;}
    cout << "Running with parameters:\n" << "max_id: " << max_id << ", id: " << id << ", step (cm): " << step << "\n";

    //range for simulating electrons, 1st sector condition (x-sqrt(3)*z<=0)&&(x+sqrt(3)*z>0) has to be satisfied later
    double xmin,xmax,ymin,ymax,zmin,zmax;
    xmin = -30; xmax = 30; ymin = -8; ymax = 8; zmin = 0; zmax = 15;

    //count the total number of electrons
    int n_electrons = 0;
    for (double x = xmin; x <= xmax; x+=step)
    {
        for (double y = ymin; y <= ymax; y+=step)
        {
            for (double z = zmin; z <= zmax; z+=step)
            {
                if((x-sqrt(3)*z<=0)&&(x+sqrt(3)*z>0)) n_electrons++;
            }
        }
    }

    int min_electron = ceil(((id-1)/(max_id*1.0))*n_electrons);  //minimal (exclusive) index of electron simulated in job with this id
    int max_electron = ceil((id/(max_id*1.0))*n_electrons);      //maximal (inclusive) index of electron simulated in job with this id

    cout << "Expecting " << n_electrons << " electrons across all simulations.\n";
    cout << "Simulating electrons from " << min_electron << " to " << max_electron << ".\n";

    //set output file, name dependent on id
    string filename = "ion" + to_string(id) + ".root";

    TFile outFile(filename.c_str(),"RECREATE","Electrons from ionization track");
    TTree electrons("electrons","Tree of initial and final points of secondary electrons");

    //setting up output TTree, contains initial and final time, position and energy of electrons
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


    // Set the gas mixture.
    MediumMagboltz gas;
    gas.SetComposition("ar", 90., "co2", 10.);

    ComponentGrid grid;
    const double m2cm = 100.;
    grid.LoadMagneticField("VecB.txt", "xyz", m2cm); 
    grid.LoadElectricField("VecE.txt", "xyz",false,false, m2cm);
    grid.SetMedium(&gas);


    // 8 cm drift gap. //constexpr double yDrift = 8; (already set up with ymax/ymin)
    //extra space so electrons end on ymax -> sensor limits
    double space = 3;

    // Assemble a sensor.
    Sensor sensor;
    sensor.AddComponent(&grid); 
    sensor.SetArea(xmin-space, ymin-space, xmin-space, xmax+space, ymax, zmax+space);


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

    int i = 0;
    for (double x = xmin; x <= xmax; x+=step)
    {
        for (double y = ymin; y <= ymax; y+=step)
        {
            for (double z = zmin; z <= zmax; z+=step)
            {
                //only inside 1 sector
                if((x-sqrt(3)*z<=0)&&(x+sqrt(3)*z>0))
                {
                    i++;
                    //if this electron is supposed to be simulated by job with this id
                    if ((i > min_electron)&&(i<=max_electron))
                    {
                        int status;
                        aval.AvalancheElectron(x, y, z, 0, 0.1, 0, 0, 0);
                        aval.GetElectronEndpoint(0,x0,y0,z0,t0,e0,x1,y1,z1,t1,e1,status);
                        electrons.Fill();
                    }
                }
            }
        }
    }

    outFile.Write();
    outFile.Close();

    return 0;
}
//plotting
    // gStyle->SetPadRightMargin(0.15);
    // TCanvas* c2 = new TCanvas("c2", "", 800, 800);
    // driftView.SetPlane(1,0, 0, 0, 0, 0);
    // driftView.SetArea(-15, -8, 15, 8);
    // driftView.SetCanvas(c2);
    // constexpr bool twod = true;
    // driftView.Plot(twod);

    // electrons.Print();

    // TCanvas* c1 = new TCanvas("c1", "", 800, 800);
    // electrons.Draw("z0:t1","","ap");

    // app.Run(true);