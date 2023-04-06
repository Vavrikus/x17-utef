#include "TH3F.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TTree.h"

#include "RK4.h"
#include "Track.h"
#include "../X17Utilities.h"

double RandomMinMax(TRandom3* rand, double min, double max) {return min + (max-min)*rand->Rndm();}

int rk_tracks()
{
    constexpr int n_tracks  =  100000;             // number of tracks to be simulated by Runge-Kutta
    constexpr double step   =  1E-13;              // step Runge-Kutta [ns]

    // assuming target in YZ plane
    constexpr double x0     =  0;                  // x coordinate of simulated window point [cm]
    constexpr double r_max  =  X17::target_radius; // maximal distance from origin [cm]

    constexpr double x1     =  X17::xmin;          // x coordinate of simulated window point [cm]
    constexpr double y1_min = -X17::win_width/2;   // minimal y coordinate of simulated window point [cm]
    constexpr double y1_max = -y1_min;             // maximal y coordinate of simulated window point [cm]
    constexpr double z1_min = -X17::win_height/2;  // minimal z coordinate of simulated window point [cm]
    constexpr double z1_max = -z1_min;             // maximal z coordinate of simulated window point [cm]

    constexpr double e_min  =  4e+6;               // minimal simulated energy [eV]
    constexpr double e_max  =  12e+6;              // maximal simulated energy [eV]


    VectorField* magfield = new VectorField(-0.2,0.2,-0.3,0.3,-0.3,0.3,0.005);
    magfield->LoadField("../mag_data/VecB2.txt");

    TTree* simulated_tracks = new TTree("rk_tracks","Runge-Kutta simulated tracks");
    TRandom3* rand = new TRandom3(0);

    TrackRK current_track;
    Vector origin,orientation;
    
    simulated_tracks->Branch("track",&current_track);
    // simulated_tracks->Branch("origin",&origin);
    // simulated_tracks->Branch("orientation",&orientation);

    std::vector<TGraph2D*> track_graphs;

    for (int i = 0; i < n_tracks; i++)
    {
        if((100*i)%n_tracks == 0) std::cout << "Progress: " << 100*i/n_tracks << " \%\n";

        double phi = RandomMinMax(rand,0,2*M_PI);
        double r   = sqrt(RandomMinMax(rand,0,r_max*r_max));
        double y0  = r*cos(phi);
        double z0  = r*sin(phi);

        double y1 = RandomMinMax(rand,y1_min,y1_max);
        double z1 = RandomMinMax(rand,z1_min,z1_max);

        origin = {x1,y1,z1};
        orientation = {x1-x0,y1-y0,z1-z0};
        orientation.Normalize();

        double kin_en = RandomMinMax(rand,e_min,e_max);

        bool electron = rand->Rndm() > 0.5;

        RK4<8>* track = GetTrackRK(magfield,electron,step,kin_en,origin,orientation);
        track->Run();

        std::vector<Matrix<8,1>> results = track->GetResults();
        std::vector<DataPoint> points;
        for (auto r : results) points.emplace_back(X17::m2cm*r.at(1,0),X17::m2cm*r.at(2,0),X17::m2cm*r.at(3,0),1);

        current_track = TrackRK(electron,points,origin,orientation,kin_en);
        simulated_tracks->Fill();

        if((100*i)%n_tracks == 0) track_graphs.push_back(GetGraphRK(track));
    }

    double height = 8;
    TCanvas* c_tracks = new TCanvas("c_tracks","Example tracks");
    // histogram for scalling axes
    TH3F* scale = new TH3F("scale","Example tracks;x [cm];y [cm];z [cm]",1,X17::xmin,X17::xmax,1,-X17::yhigh,X17::yhigh,1,-height,height);
    scale->Draw("");
    gStyle->SetOptStat(0);
    scale->GetXaxis()->SetTitleOffset(1.5);
    scale->GetZaxis()->SetTitleOffset(1.5);
    for(auto g : track_graphs) 
    {
        g->SetLineColor(kRed);
        g->Draw("LINE same");
    }
    
    X17::DrawPads3D(height);
    
    TFile* outfile = new TFile("rk_tracks.root","RECREATE");
    simulated_tracks->Write();

    return 0;
}

int main()
{
    return rk_tracks();
}