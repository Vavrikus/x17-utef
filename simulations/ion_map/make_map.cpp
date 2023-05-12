// C++ dependencies
#include <iostream>
#include <string>

// ROOT dependencies
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TStyle.h"

// X17 dependencies
#include "Field.h"
#include "MapTask.h"
#include "PadLayout.h"
#include "Points.h"
#include "Utilities.h"
#include "X17Utilities.h"

#include "MapTasks.h"

/// @brief Loads ionization electron data from files (named ion(id).root) in a specified folder.
/// @param max_id The maximum ID number of the files to load.
/// @param folder The folder containing the files to load. Default is "../../data/ion_map/sample_1.0/".
/// @return A TChain pointer containing the loaded data.
TChain* LoadData(int max_id, std::string folder = "../../data/ion_map/sample_1.0/")
{
    TChain* map_data = new TChain("map_data","Data from ionization electrons simulation.");

    for (int i = 1; i <= max_id; i++)
    {
        std::string filepath = folder + "ion" + std::to_string(i) + ".root?#electrons";
        map_data->Add(filepath.c_str());
    }

    return map_data;
}

int make_map()
{    
    // Load data from all files (results of individual jobs).
    TChain* map_data_in = LoadData(200);
    std::cout << "Number of simulated electrons: " << map_data_in->GetEntries() << "\n";

    // Set branches for TChain containing data.
    X17::MicroPoint point; // An object that will hold the information about current ionization electron loaded.
    point.SetTChainBranches(map_data_in);

    // Prepare the field that will hold the final values.
    X17::Field<X17::MapPoint> map({0,-30,-8},{15,30,8},1,X17::MapPoint());

    // Variables for checking with the previous position.
    X17::Vector v_prev(0,0,0);        // Vector for the previous position. Used for comparison with current position.
    int same_prev = 1;                // The number of entries with the same position including the current entry.
    X17::EndPoint p_avg(0,0,0,0);     // The object used for calculating of the averages.
    std::vector<X17::EndPoint> p_vec; // Stores the endpoints with the same initial position.

    // Calculate the averages and standard deviations.
    for (int i = 0; i < map_data_in->GetEntries(); i++)
    {
        map_data_in->GetEntry(i);
        
        double percent_complete = 100 * i * 1.0 / map_data_in->GetEntries();
        if (floor(percent_complete) - floor(percent_complete * (i - 1.0) / i) == 1) std::cout << "Progress: " << floor(percent_complete) << "\%\n";

        if (point.GetInitPos() == v_prev && (i != map_data_in->GetEntries()-1))
        {
            p_avg += point.end;
            p_vec.push_back(point.end);
            same_prev++;
        }

        else
        {
            if(i == map_data_in->GetEntries()-1) 
            {
                same_prev++;
                p_vec.push_back(point.end);
                p_avg = point.end; // CHECK IF CORRECT!
            }

            p_avg /= same_prev;

            if(i != 0) 
            {
                map.SetPoint(v_prev,X17::MapPoint(p_avg,stdev(p_vec,p_avg)));
                p_vec.clear();
            }

            p_vec.push_back(point.end);
            p_avg  = point.end;
            v_prev = point.GetInitPos();
            same_prev = 1;
        }
    }

    // Save the compiled map.
    TFile* outfile = new TFile("../../data/ion_map/map.root","RECREATE");
    outfile->WriteObject(&map,"map");

    // Plotting.
    bool MakePlots = true;

    if(MakePlots)
    {
        gStyle->SetOptStat(0);

        // Plotting limits for some of the plots.
        double ymin = -10; // Minimal plotted y-coordinate.
        double ymax = 10;  // Maximal plotted y-coordinate.
        double xmin = 5;   // Minimal plotted x-coordinate.
        double xmax = 16;  // Maximal plotted x-coordinate.
        
        std::vector<MapTask*> plot_tasks; // The vector containing all plotting tasks to be plotted.
        plot_tasks.push_back(new Hist_YX_DX(&map));
        plot_tasks.push_back(new Hist_YX_DY(&map));
        plot_tasks.push_back(new Hist_YX_T1(&map,xmin,xmax,ymin,ymax));
        plot_tasks.push_back(new Graph_YX(&map,xmin,xmax,ymin,ymax));
        plot_tasks.push_back(new Graph_ZT(&map));
        plot_tasks.push_back(new Graph_XZ(&map));
        plot_tasks.push_back(new Graph_XT(&map));
        plot_tasks.push_back(new Hist_XZ_T1(&map));

        for(MapTask* t : plot_tasks) t->PreLoop();

        for (int z_xyi = 0; z_xyi <= map.GetZCells(); z_xyi++)
        {
            double z = map.GetZMin() + z_xyi * map.GetStep();
            for(MapTask* t : plot_tasks) t->Z_Loop_Start(z);

            for (int xi = 0; xi <= map.GetXCells(); xi++)
            {
                double x = map.GetXMin()+xi*map.GetStep();

                for (int yi = 0; yi <= map.GetYCells(); yi++)
                {
                    double y = map.GetYMin() + yi * map.GetStep();
                    X17::MapPoint& current = map.at(xi,yi,z_xyi);

                    for(MapTask* t : plot_tasks) t->XYZ_Loop(x,y,z,current);
                }
            }
            for(MapTask* t : plot_tasks) t->Z_Loop_End();
        }

        double y_xz = 0;                                       // The y-coordinate of xz plots.
        int y_xzi   = round(y_xz-map.GetYMin())/map.GetStep(); // The y-coordinate index of xz plots.

        for (int xi = 0; xi <= map.GetXCells(); xi++)
        {
            double x = map.GetXMin() + xi * map.GetStep();
            for (int zi = 0; zi <= map.GetZCells(); zi++)
            {
                double z = map.GetZMin() + zi * map.GetStep();
                X17::MapPoint& current = map.at(xi,y_xzi,zi);

                for(MapTask* t : plot_tasks) t->ZX_Loop(z,x,current);
            }
        }

        for(MapTask* t : plot_tasks) t->PostLoop();

        outfile->Close();

        // Drawing the distortion of the pads.
        TFile* mapfile = new TFile("../../data/ion_map/map.root");
        X17::Field<X17::MapPoint>* map = (X17::Field<X17::MapPoint>*)mapfile->Get("map");
        TCanvas* c_pads = new TCanvas("c_pads", "Pads distortion for different times.");
        c_pads->Divide(4,4);

        for (int i = 0; i < 16; i++)
        {
            using namespace X17::constants;
            
            c_pads->cd(i+1);
            TGraph* gr = new TGraph();
            gr->AddPoint(-yhigh,xmin);
            gr->AddPoint(yhigh,xmax);
            gr->Draw("AP");

            X17::DefaultLayout& pads = X17::DefaultLayout::GetDefaultLayout();
            pads.DrawPadsDistortion((i + 1) * 5000 / 16, c_pads, map);
            // X17::DrawTrapezoid();
        }
    }

    return 0;
}

int main(int argc, char const *argv[])
{
    return make_map();
}