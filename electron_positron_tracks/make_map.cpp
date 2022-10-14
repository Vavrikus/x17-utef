// my dependencies
#include "../VectorField.h"

// C++ dependencies
#include <iostream>
#include <string>

// ROOT dependencies
#include "TChain.h"
#include "TFile.h"

using namespace std;

// loads data (named ion(id).root) from folder
TChain* LoadData(int max_id, string folder = "/home/vavrik/work/X17/electron_positron_tracks/data/")
{
    TChain* map_data = new TChain("map_data","Data from ionization electrons simulation.");

    for (int i = 1; i <= max_id; i++)
    {
        string filepath = folder + "ion" + to_string(i) + ".root?#electrons";
        map_data->Add(filepath.c_str());
    }

    return map_data;
}

int make_map()
{
    // load data from all files (results of individual jobs)
    TChain* map_data_in = LoadData(200);

    //map_data_in->Print();
    cout << "\n\nEntries: " << map_data_in->GetEntries() << "\n";

    // set branches for tchain containing data
    double x0, y0, z0, t0, e0;
    double x1, y1, z1, t1, e1;
    map_data_in->SetBranchAddress("x0",&x0);
    map_data_in->SetBranchAddress("y0",&y0);
    map_data_in->SetBranchAddress("z0",&z0);
    map_data_in->SetBranchAddress("t0",&t0);
    map_data_in->SetBranchAddress("e0",&e0);
    map_data_in->SetBranchAddress("x1",&x1);
    map_data_in->SetBranchAddress("y1",&y1);
    map_data_in->SetBranchAddress("z1",&z1);
    map_data_in->SetBranchAddress("t1",&t1);
    map_data_in->SetBranchAddress("e1",&e1);

    // prepare field for output
    Field<SensorData> map;
    map.SetDefault({0,0,0,0});
    map.InitField(-30,30,-8,8,0,15,1);

    // variables for checking with the previous position
    double x0_prev = 0; double y0_prev = 0; double z0_prev = 0;
    int same_prev = 1;
    
    // variables for averages
    double x1_avg = 0; double y1_avg = 0; double z1_avg = 0; double t1_avg = 0;
    
    // calculate averages
    for (int i = 0; i < map_data_in->GetEntries(); i++)
    {
        map_data_in->GetEntry(i);
        cout << "i: " << i << "\n";

        if ((x0 == x0_prev) && (y0 == y0_prev) && (z0 == z0_prev) && (i != map_data_in->GetEntries()-1))
        {
            x1_avg += x1; y1_avg += y1; z1_avg += z1; t1_avg = t1;
            same_prev++;
        }

        else
        {
            x1_avg /= same_prev; y1_avg /= same_prev; z1_avg /= same_prev;
            if(i != 0) 
                map.SetPoint(x0_prev,y0_prev,z0_prev,{x1,y1,z1,t1});
            // cout << "same_prev: " << same_prev << " x0_prev: " << x0_prev << " y0_prev: " << y0_prev  << " z0_prev: " << z0_prev << "\n";
            // cout << "x1_avg: " << x1_avg << " y1_avg: " << y1_avg  << " z1_avg: " << z1_avg << " t1_avg: " << t1_avg << "\n";

            x0_prev = x0; y0_prev = y0; z0_prev = z0;
            same_prev = 1;
        }
        
        // cout << "x0: " << x0 << " y0: " << y0  << " z0: " << z0  << " t0: " << t0  << " e0: " << e0 << "\n"; 
        // cout << "x1: " << x1 << " y1: " << y1  << " z1: " << z1  << " t1: " << t1  << " e1: " << e1 << "\n";
    }

    // output
    TFile* outfile = new TFile("map.root","RECREATE");
    outfile->WriteObject(&map,"map");
    outfile->Close();

    return 0;
}

int main(int argc, char const *argv[])
{
    return make_map();
}