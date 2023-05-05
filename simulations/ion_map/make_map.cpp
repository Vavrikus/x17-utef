// C++ dependencies
#include <iostream>
#include <string>

// ROOT dependencies
#include "TArrow.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH2F.h"
#include "TStyle.h"

// X17 dependencies
#include "Field.h"
#include "MapTask.h"
#include "PadLayout.h"
#include "Points.h"
#include "Utilities.h"
#include "X17Utilities.h"

int main(int argc, char const *argv[])
{
    return make_map();
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

    // variables for checking with the previous position
    X17::Vector v_prev(0,0,0);        // Vector for the previous position. Used for comparison with current position.
    int same_prev = 1;                // The number of entries with the same position including the current entry.
    X17::EndPoint p_avg(0,0,0,0);     // The object used for calculating of the averages.
    std::vector<X17::EndPoint> p_vec; // Stores the endpoints with the same initial position.

    // Calculate the averages and standard deviations.
    for (int i = 0; i < map_data_in->GetEntries(); i++)
    {
        map_data_in->GetEntry(i);
        
        double percent_complete = 100*i*1.0/map_data_in->GetEntries();
        if (floor(percent_complete) - floor(percent_complete*(i-1.0)/i) == 1) std::cout << "Progress: " << floor(percent_complete) << "\%\n";

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
            double z = map.GetZMin() + z_xyi*map.GetStep();
            for(MapTask* t : plot_tasks) t->Z_Loop_Start(z);

            for (int xi = 0; xi <= map.GetXCells(); xi++)
            {
                double x = map.GetXMin()+xi*map.GetStep();

                for (int yi = 0; yi <= map.GetYCells(); yi++)
                {
                    double y = map.GetYMin()+yi*map.GetStep();
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
            double x = map.GetXMin()+xi*map.GetStep();
            for (int zi = 0; zi <= map.GetZCells(); zi++)
            {
                double z = map.GetZMin()+zi*map.GetStep();
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
            pads.DrawPadsDistortion((i+1)*5000/16,c_pads,map);
            //X17::DrawTrapezoid();
        }
    }

    return 0;
}

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

class Hist_YX_DX : public MapTask
{
    std::vector<TH2F*> v_yx_dx;

public:
    Hist_YX_DX(X17::Field<X17::MapPoint>* map) : MapTask(map) { }

    void Z_Loop_Start(const double& z) override
    {
        std::string h_dx_name = "h_yx_dx_" + std::to_string(z);
        std::string h_dx_title = "XY plane dx (z = " + std::to_string(z) + ");x [cm];y [cm]";
        v_yx_dx.push_back(new TH2F(h_dx_name.c_str(),h_dx_title.c_str(),map->GetYCells(),map->GetYMin(),map->GetYMax(),map->GetXCells(),map->GetXMin(),map->GetXMax()));
    }

    void XYZ_Loop(const double& x, const double& y, const double& z, const X17::MapPoint& current) override
    {
        if(current.x != 0) v_yx_dx.back()->Fill(y,x,current.x-x);
    }

    void PostLoop() override
    {
        for(TH2F* h : v_yx_dx) h->Write();
    }
};

class Hist_YX_DY : public MapTask
{
    std::vector<TH2F*> v_yx_dy;

public:
    Hist_YX_DY(X17::Field<X17::MapPoint>* map) : MapTask(map) { }

    void Z_Loop_Start(const double& z) override
    {
        std::string h_dy_name = "h_yx_dy_" + std::to_string(z);
        std::string h_dy_title = "XY plane dy (z = " + std::to_string(z) + ");x [cm];y [cm]";
        v_yx_dy.push_back(new TH2F(h_dy_name.c_str(),h_dy_title.c_str(),map->GetYCells(),map->GetYMin(),map->GetYMax(),map->GetXCells(),map->GetXMin(),map->GetXMax()));
    }

    void XYZ_Loop(const double& x, const double& y, const double& z, const X17::MapPoint& current) override
    {
        if(current.y != 0) v_yx_dy.back()->Fill(y,x,current.y-y);
    }

    void PostLoop() override
    {
        for(TH2F* h : v_yx_dy) h->Write();
    }
};


class Hist_YX_T1 : public MapTask
{
    double xmin,xmax,ymin,ymax;
    std::vector<TH2F*> v_yx_t1;

public:
    Hist_YX_T1(X17::Field<X17::MapPoint>* map, double xmin, double xmax, double ymin, double ymax) : MapTask(map),xmin(xmin),xmax(xmax),ymin(ymin),ymax(ymax) { }

    void Z_Loop_Start(const double& z) override
    {
        std::string h_t1_name = "h_yx_t1_" + std::to_string(z);
        std::string h_t1_title = "XY plane t1 (z = " + std::to_string(z) + ");x [cm];y [cm]";
        v_yx_t1.push_back(new TH2F(h_t1_name.c_str(),h_t1_title.c_str(),map->GetYCells(),map->GetYMin(),map->GetYMax(),map->GetXCells(),map->GetXMin(),map->GetXMax()));
    }

    void XYZ_Loop(const double& x, const double& y, const double& z, const X17::MapPoint& current) override
    {
        if(current.t != 0) v_yx_t1.back()->Fill(y,x,current.t);
    }

    void PostLoop() override
    {
        TCanvas* c_yx_t1 = new TCanvas("c_yx_t1","Map of drift times");
        c_yx_t1->Divide(4,4);
        for (int i = 0; i < v_yx_t1.size()-1; i++)
        {
            c_yx_t1->cd(i+1);
            v_yx_t1[i]->GetXaxis()->SetRangeUser(ymin,ymax);
            v_yx_t1[i]->GetYaxis()->SetRangeUser(xmin,xmax);
            v_yx_t1[i]->Draw("colz");
            double min_percent = 0.8+0.13*(v_yx_t1.size()-1-i)/(v_yx_t1.size()-1); //how high is minimum compared to maximum
            v_yx_t1[i]->SetMinimum(min_percent*v_yx_t1[i]->GetMaximum());
            v_yx_t1[i]->SetMaximum(min_percent*v_yx_t1[i]->GetMaximum()+500);
            X17::DrawTrapezoid();
        }
        c_yx_t1->Write();
    }
};

class Graph_YX : public MapTask
{
    double xmin,xmax,ymin,ymax;
    std::vector<TGraphErrors*> v_g_yx;
    std::vector<std::vector<TArrow*>> vv_g_yx_arrows;
    std::vector<TArrow*> arrows;

public:
    Graph_YX(X17::Field<X17::MapPoint>* map, double xmin, double xmax, double ymin, double ymax) : MapTask(map),xmin(xmin),xmax(xmax),ymin(ymin),ymax(ymax) { }

    void Z_Loop_Start(const double& z) override
    {
        std::string g_yx_name = "g_yx_" + std::to_string(z);
        std::string g_yx_title = "Map of electron readout positions, z = " + std::to_string(z);

        v_g_yx.push_back(new TGraphErrors());

        v_g_yx.back()->SetName(g_yx_name.c_str());
        v_g_yx.back()->SetTitle(g_yx_title.c_str());
        v_g_yx.back()->SetMarkerColor(2);
        v_g_yx.back()->SetMarkerStyle(6);
        v_g_yx.back()->GetXaxis()->SetTitle("y [cm]");
        v_g_yx.back()->GetYaxis()->SetTitle("x [cm]");
    }

    void XYZ_Loop(const double& x, const double& y, const double& z, const X17::MapPoint& current) override
    {
        if((current.x != 0) && (current.y != 0))
        {
            v_g_yx.back()->AddPoint(current.y,current.x);
            v_g_yx.back()->SetPointError(v_g_yx.back()->GetN()-1,current.ydev,current.xdev);
            TArrow* arrow = new TArrow(y,x,current.y,current.x);
            arrows.push_back(arrow);
        }
    }

    void Z_Loop_End() override
    {
        vv_g_yx_arrows.push_back(arrows);
    }

    void PostLoop() override
    {
        // for(TGraph* g : v_g_yx) g->Write();
        TCanvas* c_g_yx = new TCanvas("c_g_yx","Map of electron readout positions");
        c_g_yx->Divide(4,4);
        for (int i = 0; i < v_g_yx.size()-1; i++)
        {
            c_g_yx->cd(i+1);
            v_g_yx[i]->GetXaxis()->SetRangeUser(ymin,ymax);
            v_g_yx[i]->GetYaxis()->SetRangeUser(xmin,xmax);
            v_g_yx[i]->Draw("AP");
            X17::DrawTrapezoid();

            for(TArrow* arr : vv_g_yx_arrows[i]) 
            {
                if(arr->GetX2() > ymin && arr->GetX2() < ymax && arr->GetY2() > xmin && arr->GetY2() < xmax)
                {
                    arr->SetArrowSize(0.002);
                    arr->Draw();
                }
            }
        }
        c_g_yx->Write();
    }
};

class Graph_ZT : public MapTask
{
    TGraph* g_zt;

public:
    Graph_ZT(X17::Field<X17::MapPoint>* map) : MapTask(map) { }

    void Z_Loop_Start(const double& z) override
    {
        g_zt = new TGraph();

        g_zt->SetName("g_zt");
        g_zt->SetTitle("Original height vs drift time");
        g_zt->SetMarkerColor(2);
        g_zt->SetMarkerStyle(6);
        g_zt->GetXaxis()->SetTitle("z [cm]");
        g_zt->GetYaxis()->SetTitle("t [ns]");
    }

    void XYZ_Loop(const double& x, const double& y, const double& z, const X17::MapPoint& current) override
    {
        g_zt->AddPoint(z,current.t);
    }

    void PostLoop() override
    {
        TCanvas* c_g_zt = new TCanvas("c_g_zt","Original height vs drift time");
        g_zt->Draw("AP");
        c_g_zt->Write();
    }
};

class Graph_XZ : public MapTask
{
    TGraphErrors* g_xz;
    std::vector<TArrow*> v_g_xz_arrows;

public:
    Graph_XZ(X17::Field<X17::MapPoint>* map) : MapTask(map) { }

    void PreLoop() override
    {
        g_xz = new TGraphErrors();
        g_xz->SetName("g_xz");
        g_xz->SetTitle("Map of electron displacement (x = 0)");
        g_xz->SetMarkerColor(2);
        g_xz->SetMarkerStyle(6);
        g_xz->GetXaxis()->SetTitle("x [cm]");
        g_xz->GetYaxis()->SetTitle("z [cm]");
    }

    void ZX_Loop(const double& z, const double& x, const X17::MapPoint& current) override
    {
        g_xz->AddPoint(current.x,z);
        g_xz->SetPointError(g_xz->GetN()-1,current.xdev,0);
        TArrow* arrow = new TArrow(x,z,current.x,z);
        v_g_xz_arrows.push_back(arrow);
    }

    void PostLoop() override
    {
        TCanvas* c_g_xz = new TCanvas("c_g_xz","Map of ionization electron displacement");
        g_xz->Draw("AP");
        TLine* lleft  = new TLine(6.51,8,6.51,-8);
        TLine* lright = new TLine(14.61,8,14.61,-8);
        lleft->Draw();lright->Draw();
        for(TArrow* arr : v_g_xz_arrows) {arr->SetArrowSize(0.004); arr->Draw();}
        c_g_xz->Write();
    }
};

class Graph_XT : public MapTask
{
    TGraphErrors* g_xt;

public:
    Graph_XT(X17::Field<X17::MapPoint>* map) : MapTask(map) { }

    void PreLoop() override
    {
        g_xt = new TGraphErrors();
        g_xt->SetName("g_xt");
        g_xt->SetTitle("Map of drift times (y = 0)");
        g_xt->SetMarkerColor(2);
        g_xt->SetMarkerStyle(6);
        g_xt->GetXaxis()->SetTitle("x [cm]");
        g_xt->GetYaxis()->SetTitle("t [ns]");
    }

    void ZX_Loop(const double& z, const double& x, const X17::MapPoint& current) override
    {        
        g_xt->AddPoint(x,current.t);
        g_xt->SetPointError(g_xt->GetN()-1,current.xdev,current.tdev);
    }

    void PostLoop() override
    {
        TCanvas* c_g_xt = new TCanvas("c_g_xt","Map of ionization electron drift times");
        g_xt->Draw("AP");
        TLine* lleft  = new TLine(6.51,0,6.51,5000);
        TLine* lright = new TLine(14.61,0,14.61,5000);
        lleft->Draw();lright->Draw();
        c_g_xt->Write();
    }
};

class Hist_XZ_T1 : public MapTask
{
    TH2F* xz_t1;

public:
    Hist_XZ_T1(X17::Field<X17::MapPoint>* map) : MapTask(map) { }

    void PreLoop() override
    {
        xz_t1 = new TH2F("h_xz_t1","XZ plane t1;x [cm];z [cm]",map->GetXCells(),map->GetXMin(),map->GetXMax(),map->GetZCells(),map->GetZMin(),map->GetZMax());
    }

    void ZX_Loop(const double& z, const double& x, const X17::MapPoint& current) override
    {        
        xz_t1->Fill(x,z,current.t);
    }

    void PostLoop() override
    {
        xz_t1->Write();
    }
};