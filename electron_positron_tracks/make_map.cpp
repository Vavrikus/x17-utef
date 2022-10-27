// my dependencies
#include "../VectorField.h"

// C++ dependencies
#include <iostream>
#include <string>
#include <vector>

// ROOT dependencies
#include "TArrow.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TStyle.h"

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

//calculates standard deviation assuming gaussian distribution
double stdev(vector<double> values, double average)
{
    double sqdev_sum = 0;
    double dmin = values[0];
    double dmax = values[0];

    for (double d : values)
    {
        if(d<dmin) dmin = d;
        if(d>dmax) dmax = d;
        sqdev_sum += pow((d-average),2);
    }

    // cout << "min: " << dmin << " average: " << average << " max: " << dmax << " N: " << values.size() << " stdev: " << sqrt(sqdev_sum/values.size()) << "\n";

    return sqrt(sqdev_sum/values.size());
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
    map.SetDefault({0,0,0,0,0,0,0,0});
    map.InitField(-30,30,-8,8,0,15,1);

    // variables for checking with the previous position
    double x0_prev = 0; double y0_prev = 0; double z0_prev = 0;
    int same_prev = 1;
    
    double x1_avg = 0; double y1_avg = 0; double z1_avg = 0; double t1_avg = 0; // variables for averages
    vector<double> x1_vec,y1_vec,z1_vec,t1_vec;                                 // variables for data from each point
    
    // calculate averages and standart deviations
    for (int i = 0; i < map_data_in->GetEntries(); i++)
    {
        map_data_in->GetEntry(i);
        
        double percent_complete = 100*i*1.0/map_data_in->GetEntries();
        if (floor(percent_complete)-floor(percent_complete*(i-1.0)/i) == 1) cout << floor(percent_complete) << "\%\n";

        if ((x0 == x0_prev) && (y0 == y0_prev) && (z0 == z0_prev) && (i != map_data_in->GetEntries()-1))
        {
            x1_avg += x1; y1_avg += y1; z1_avg += z1; t1_avg += t1;
            x1_vec.push_back(x1);y1_vec.push_back(y1);z1_vec.push_back(z1);t1_vec.push_back(t1);
            same_prev++;
        }

        else
        {
            if(i == map_data_in->GetEntries()-1) 
            {
                same_prev++;
                x1_vec.push_back(x1);y1_vec.push_back(y1);z1_vec.push_back(z1);t1_vec.push_back(t1);
                x1_avg = x1; y1_avg = y1; z1_avg = z1; t1_avg = t1;
            }

            x1_avg /= same_prev; y1_avg /= same_prev; z1_avg /= same_prev; t1_avg /= same_prev;
            if(i != 0) 
            {
                map.SetPoint(x0_prev,y0_prev,z0_prev,{x1_avg,y1_avg,z1_avg,t1_avg,stdev(x1_vec,x1_avg),stdev(y1_vec,y1_avg),stdev(z1_vec,z1_avg),stdev(t1_vec,t1_avg)});
                x1_vec.clear();y1_vec.clear();z1_vec.clear();t1_vec.clear();
            }

            x1_vec.push_back(x1);y1_vec.push_back(y1);z1_vec.push_back(z1);t1_vec.push_back(t1);
            x1_avg = x1; y1_avg = y1; z1_avg = z1; t1_avg = t1;
            x0_prev = x0; y0_prev = y0; z0_prev = z0;
            same_prev = 1;
        }
    }

    // output
    TFile* outfile = new TFile("map.root","RECREATE");
    outfile->WriteObject(&map,"map");
    // outfile->Close();

    // plotting
    bool MakePlots = true;

    if(MakePlots)
    {
        vector<TH2F*> v_xz_dx;
        vector<TH2F*> v_xz_dz;
        vector<TGraph*> v_g_xz;
        vector<vector<TArrow*>> vv_g_xz_arrows;

        TH2F* yz_t1 = new TH2F("h_yz_t1","YZ plane t1;z [cm];y [cm]",map.zimax,map.zmin,map.zmax,map.yimax,map.ymin,map.ymax);

        gStyle->SetOptStat(0);

        double y_xz = 7;
        double x_yz = 0;

        // int y_xzi = round(y_xz-map.ymin)/map.step;
        int x_yzi = round(x_yz-map.xmin)/map.step;
        for (int y_xzi = 0; y_xzi <= map.yimax; y_xzi++)
        {
            double y = map.ymin + y_xzi*map.step;
            string h_dx_name = "h_xz_dx_" + to_string(y);
            string h_dz_name = "h_xz_dz_" + to_string(y);
            string g_xz_name = "g_xz_" + to_string(y);
            string h_dx_title = "XZ plane dx (y = " + to_string(y) + ");x [cm];z [cm]";
            string h_dz_title = "XZ plane dx (y = " + to_string(y) + ");x [cm];z [cm]";
            string g_xz_title = "Map of electron readout positions, y = " + to_string(y);

            v_xz_dx.push_back(new TH2F(h_dx_name.c_str(),h_dx_title.c_str(),map.ximax,map.xmin,map.xmax,map.zimax,map.zmin,map.zmax));
            v_xz_dz.push_back(new TH2F(h_dz_name.c_str(),h_dz_title.c_str(),map.ximax,map.xmin,map.xmax,map.zimax,map.zmin,map.zmax));
            v_g_xz.push_back(new TGraph());

            v_g_xz[y_xzi]->SetName(g_xz_name.c_str());
            v_g_xz[y_xzi]->SetTitle(g_xz_title.c_str());
            v_g_xz[y_xzi]->SetMarkerColor(2);
            v_g_xz[y_xzi]->SetMarkerStyle(6);
            v_g_xz[y_xzi]->GetXaxis()->SetTitle("x [cm]");
            v_g_xz[y_xzi]->GetYaxis()->SetTitle("z [cm]");

            vector<TArrow*> arrows;

            for (int zi = 0; zi <= map.zimax; zi++)
            {
                double z = map.zmin+zi*map.step;

                for (int xi = 0; xi <= map.ximax; xi++)
                {
                    double x = map.xmin+xi*map.step;
                    SensorData& s_curr = map.field[xi][y_xzi][zi];
                    if(s_curr.x1 != 0) v_xz_dx[y_xzi]->Fill(x,z,s_curr.x1-x);
                    if(s_curr.z1 != 0) v_xz_dz[y_xzi]->Fill(x,z,s_curr.z1-z);
                    if((s_curr.x1 != 0) && (s_curr.z1 != 0))
                    {
                        v_g_xz[y_xzi]->AddPoint(s_curr.x1,s_curr.z1);
                        TArrow* arrow = new TArrow(x,z,s_curr.x1,s_curr.z1);
                        arrows.push_back(arrow);
                    }
                }
            }

            vv_g_xz_arrows.push_back(arrows);
        }

        for (int zi = 0; zi <= map.zimax; zi++)
        {
            double z = map.zmin+zi*map.step;
            for (int yi = 0; yi <= map.yimax; yi++)
            {
                double y = map.ymin+yi*map.step;
                yz_t1->Fill(z,y,map.field[x_yzi][yi][zi].t1);
            }
        }
        yz_t1->Write();
        for(TH2F* h : v_xz_dx) h->Write();
        for(TH2F* h : v_xz_dz) h->Write();
        // for(TGraph* g : v_g_xz) g->Write();

        TCanvas* c_g_xz = new TCanvas("c_g_xz","Map of electron readout positions");
        double xmin = -10; double xmax = 10; double zmin = 5; double zmax = 16;
        c_g_xz->Divide(4,4);
        for (int i = 0; i < v_g_xz.size()-1; i++)
        {
            c_g_xz->cd(i+1);
            v_g_xz[i]->GetXaxis()->SetRangeUser(xmin,xmax);
            v_g_xz[i]->GetYaxis()->SetRangeUser(zmin,zmax);
            v_g_xz[i]->Draw("AP");
            TLine* l1 = new TLine(-7.45,14.61,7.45,14.61);
            TLine* l2 = new TLine(-2.25,6.51,2.25,6.51);
            TLine* l3 = new TLine(-7.45,14.61,-2.25,6.51);
            TLine* l4 = new TLine(7.45,14.61,2.25,6.51);
            l1->Draw();
            l2->Draw();
            l3->Draw();
            l4->Draw();
            for(TArrow* arr : vv_g_xz_arrows[i]) 
            {
                if(arr->GetX2() > xmin && arr->GetX2() < xmax && arr->GetY2() > zmin && arr->GetY2() < zmax)
                {
                    arr->SetArrowSize(0.002);
                    arr->Draw();
                }
            }
        }
        c_g_xz->Write();
        
    }

    return 0;
}

int main(int argc, char const *argv[])
{
    return make_map();
}