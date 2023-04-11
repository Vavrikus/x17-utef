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

struct Job
{
    // most important parameters are jobs id, maximal id in given run,
    // step in centimeters and number of iterations per position
    int id, max_id, iterations;
    double step;

    // ranges for simulation
    // 1st sector condition (y-sqrt(3)*x<=0)&&(y+sqrt(3)*x>0) has to be satisfied later
    double xmin = 0;
    double xmax = 15;
    double ymin = -30;
    double ymax = 30;
    double zmin = -8;
    double zmax = 8;

    int nx,nz,nxy;     // variables for number of steps in x, z, x and y combined
    int n_el;          // total number of electrons in all simulations
    int min_el,max_el; // min/max indexes of electrons simulated in this simulation (inclusive)

    double zsum;       // total distance to propagate all electrons

    void GetParameters(int argc, char* argv[])
    {
        // check number of paramaters passed to main function
        if (argc < 4)
        {
            cerr << "ERROR: Missing arguments in ion_electrons. Correct arguments: max_id, id, step (cm).\n";
        }

        if (argc < 5) iterations = 1;
        else          iterations = stoi(argv[4]);

        // set parameters of program
        max_id = stoi(argv[1]);
        id     = stoi(argv[2]);
        step   = stod(argv[3]);

        if (id <= 0)         cerr << "ERROR: Parameter id has to be positive.\n";
        if (id > max_id)     cerr << "ERROR: Parameter id cannot be bigger than max_id.\n";
        if (step <= 0.0)     cerr << "ERROR: Parameter step has to be positive.\n";
        if (iterations <= 0) cerr << "ERROR: Parameter iterations has to be positive.\n";

        cout << "Running with parameters:\n" << "max_id: " << max_id << ", id: " << id << ", step (cm): " << step << ", iterations: " << iterations << ".\n";
    }

    void SetElectronBounds()
    {
        f_GetStepParameters();
        n_el = nz*nxy;
        zsum = (zmax-zmin)*nxy*nz-nxy*step*nz*(nz-1)/2;

        cout << "Expecting " << n_el << " electrons across all simulations, total propagation distance (per iteration) " << zsum << ".\n";

        min_el = f_Noptimal(id-1)+1;
        max_el = f_Noptimal(id);

        cout << "Simulating electrons from " << min_el << " to " << max_el << ", total propagation (per iteration) " << f_Zsum(max_el)-f_Zsum(min_el-1) << ".\n";
    }

private:
    void f_GetStepParameters()
    {
        nx = floor((xmax-xmin)/step)+1; // number of x steps in each iteration
        nz = floor((zmax-zmin)/step)+1; // number of z steps in each iteration
        nxy = 0;                        // total number of x and y steps for each z step

        // iterate over x and add up y steps to get nxy
        for (int i = 0; i < nx; i++) 
        {
            double x     = xmin + i*step;                   // current x
            int nymin    = floor((-ymin-x*sqrt(3))/step)+1; // minimal number of y steps to satisfy minimal sector condition
            double ymin2 = ymin + nymin*step;               // smallest y to satisfy minimal sector condition

            double ny; // number of y steps for current x iteration
            if (ymin2 > (x*sqrt(3))) ny = 0; // if maximal sector condition cannot be satisfied
            else ny = floor((x*sqrt(3)-ymin2)/step)+1;
            nxy += ny;
        }
    }

    // calculates total z to propagate for given number of electrons
    double f_Zsum(int n_electrons)
    {
        int full_layers = floor(n_electrons/nxy); // full z layers
        int remainder   = n_electrons%nxy;        // electrons in incomplete z layer on top
        return (zmax-zmin)*nxy*full_layers - nxy*step*full_layers*(full_layers-1)/2 + remainder*(zmax-zmin-full_layers*step);
    }

    // finds best maximal index of electron for this id
    int f_Noptimal(int l_id)
    {
        if (l_id == 0) return 0; // id 0 is not used but best minimal index of id 1 depends on this

        int min_el = 0;
        int max_el = nz*nxy;
        double opt_zsum = zsum/max_id; // optimal total propagation for one job

        int current_el = max_el/2;
        while (min_el != current_el)
        {
            if(f_Zsum(current_el) > l_id*opt_zsum) max_el = current_el;
            else min_el = current_el;
            current_el = (max_el+min_el)/2; 
        }  

        return max_el;
    }
};

//arguments: max_id, id, step (cm)
int main(int argc, char* argv[])
{
    TApplication app("app", &argc, argv);

    Job job;
    job.GetParameters(argc,argv);
    job.SetElectronBounds();

    //set output file, name dependent on id
    string filename = "ion" + to_string(job.id) + ".root";

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
    gas.SetComposition("ar", 90., "co2", 10.); cout << "BAD GAS COMPOSITION!!!!\n";

    ComponentGrid grid;
    const double m2cm = 100.;
    grid.LoadMagneticField("../../data/elmag/VecB2.txt", "xyz", m2cm); 
    grid.LoadElectricField("../../data/elmag/VecE2.txt", "xyz",false,false, m2cm);
    grid.SetMedium(&gas);


    // 8 cm drift gap. //constexpr double zDrift = 8; (already set up with zmax/zmin)
    //extra space so electrons end on zmax -> sensor limits
    double space = 3;

    // Assemble a sensor.
    Sensor sensor;
    sensor.AddComponent(&grid); 
    sensor.SetArea(job.xmin-space, job.ymin-space, job.zmin-space, job.xmax+space, job.ymax+space, job.zmax);


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

    int i_el = 0; // index of electron
    for (double z = job.zmin; z <= job.zmax; z+=job.step)
    {
        for (double y = job.ymin; y <= job.ymax; y+=job.step)
        {
            for (double x = job.xmin; x <= job.xmax; x+=job.step)
            {
                //only inside 1 sector
                if((y-sqrt(3)*x<=0)&&(y+sqrt(3)*x>0))
                {
                    i_el++;

                    //if this electron is supposed to be simulated by job with this id
                    if (!((i_el < job.min_el)||(i_el > job.max_el)))
                    {
                        for (int j = 0; j < job.iterations; j++)
                        {                        
                            int status;
                            aval.AvalancheElectron(x, y, z, 0, 0.1, 0, 0, 0);
                            aval.GetElectronEndpoint(0,x0,y0,z0,t0,e0,x1,y1,z1,t1,e1,status);
                            electrons.Fill();

                            //cout << "x: " << x << " y: " << y << " z: " << z << " i_el: " << i_el << "\n";
                        }
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