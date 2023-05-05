// ROOT dependencies
#include "TCanvas.h"
#include "TH3F.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TTree.h"

// X17 dependencies
#include "Field.h"
#include "PadLayout.h"
#include "RK4.h"
#include "Track.h"
#include "Utilities.h"
#include "X17Utilities.h"

using namespace X17::constants;

int rk_tracks()
{
    constexpr int n_tracks  =  100000;        // The number of tracks to be simulated by Runge-Kutta.
    constexpr double step   =  1E-13;         // The step of Runge-Kutta [s].

    // Assuming that the target is in the YZ plane.
    constexpr double x0     =  0;             // The x coordinate of simulated origin [cm].
    constexpr double r_max  =  target_radius; // The maximal distance from origin [cm].

    constexpr double x1     =  xmin;          // The x coordinate of simulated window point [cm].
    constexpr double y1_min = -win_width/2;   // The minimal y coordinate of simulated window point [cm].
    constexpr double y1_max = -y1_min;        // The maximal y coordinate of simulated window point [cm].
    constexpr double z1_min = -win_height/2;  // The minimal z coordinate of simulated window point [cm].
    constexpr double z1_max = -z1_min;        // The maximal z coordinate of simulated window point [cm].

    constexpr double e_min  =  4e+6;          // The minimal simulated energy [eV].
    constexpr double e_max  =  12e+6;         // The maximal simulated energy [eV].

    // Loading the magnetic field.
    X17::Field<X17::Vector>& magfield = X17::LoadField("../../data/elmag/VecB2.txt",{-0.2,-0.3,-0.3},{0.2,0.3,0.3},0.005);

    // Some necessary variables for simulating and saving track parameters.
    TTree* simulated_tracks = new TTree("rk_tracks","Runge-Kutta simulated tracks");
    TRandom3* rand = new TRandom3(0);

    X17::TrackRK current_track;
    X17::Vector origin,orientation;
    
    simulated_tracks->Branch("track",&current_track);

    std::vector<TGraph2D*> track_graphs;

    // The loop for the track simulation.
    for (int i = 0; i < n_tracks; i++)
    {
        if ((100 * i) % n_tracks == 0) std::cout << "Progress: " << 100*i/n_tracks << " \%\n";

        // Simulation of the initial track parameters.
        double phi = RandomMinMax(rand,0,2*M_PI);            // The azimuth of the initial point on the circle target.
        double r   = sqrt(RandomMinMax(rand,0,r_max*r_max)); // The distance of the initial point from the circle target.
        double y0  = r*cos(phi);                             // The y-coordinate of the initial point.
        double z0  = r*sin(phi);                             // The z-coordinate of the initial point.

        double y1 = RandomMinMax(rand,y1_min,y1_max);        // The y-coordinate of the window point (TPC entry).
        double z1 = RandomMinMax(rand,z1_min,z1_max);        // The z-coordinate of the window point (TPC entry).

        origin = {x1, y1, z1};                               // Setting the initial point (or origin).
        orientation = {x1 - x0, y1 - y0, z1 - z0};           // Setting the initial orientation.
        orientation.Normalize();

        double kin_en = RandomMinMax(rand,e_min,e_max);      // The kinetic energy of the particle.

        bool electron = rand->Rndm() > 0.5;                  // Choosing either electron or positron.

        // The actual track simulation.
        X17::RK4<8>* track = GetTrackRK(magfield,electron,step,kin_en,origin,orientation);
        track->Run();

        std::vector<X17::Matrix<8,1>> results = track->GetResults();
        std::vector<X17::RKPoint> points;

        using namespace X17::constants;
        for (auto r : results) points.emplace_back(m2cm*r.at(1,0),m2cm*r.at(2,0),m2cm*r.at(3,0),1e+9*r.at(0,0));

        current_track = X17::TrackRK(electron,points,origin,orientation,kin_en);
        simulated_tracks->Fill();

        if((100*i)%n_tracks == 0) track_graphs.push_back(GetGraphRK(track));
    }

    // Plotting some of the tracks.
    double height = 8;
    using namespace X17::constants;
    TCanvas* c_tracks = new TCanvas("c_tracks","Example tracks");
    // A histogram for scalling of the axes.
    TH3F* scale = new TH3F("scale","Example tracks;x [cm];y [cm];z [cm]",1,xmin,xmax,1,-yhigh,yhigh,1,-height,height);
    scale->Draw("");
    gStyle->SetOptStat(0);
    scale->GetXaxis()->SetTitleOffset(1.5);
    scale->GetZaxis()->SetTitleOffset(1.5);
    for(auto g : track_graphs) 
    {
        g->SetLineColor(kRed);
        g->Draw("LINE same");
    }
    
    X17::DefaultLayout::GetDefaultLayout().DrawPads3D(height);
    
    TFile* outfile = new TFile("../../data/rk_tracks/rk_tracks.root","RECREATE");
    simulated_tracks->Write();

    return 0;
}

int main()
{
    return rk_tracks();
}