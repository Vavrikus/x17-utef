// my dependencies
#include "../VectorField.h"
#include "../X17Utilities.h"

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
#include "TGraph2D.h"
#include "TGraphErrors.h"
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
    // ADJUSTED FOR OLD COORDINATES ONLY!!!!!!!
    cout << "MAKE SURE TO CHANGE THE COORDINATE ASSIGNMENT IF USING NEW SIMULATIONS!!!!!!!\n";
    double x0, y0, z0, t0, e0;
    double x1, y1, z1, t1, e1;
    map_data_in->SetBranchAddress("z0",&x0);
    map_data_in->SetBranchAddress("x0",&y0);
    map_data_in->SetBranchAddress("y0",&z0);
    map_data_in->SetBranchAddress("t0",&t0);
    map_data_in->SetBranchAddress("e0",&e0);
    map_data_in->SetBranchAddress("z1",&x1);
    map_data_in->SetBranchAddress("x1",&y1);
    map_data_in->SetBranchAddress("y1",&z1);
    map_data_in->SetBranchAddress("t1",&t1);
    map_data_in->SetBranchAddress("e1",&e1);

    // prepare field for output
    Field<SensorData> map;
    map.SetDefault({0,0,0,0,0,0,0,0});
    map.InitField(0,15,-30,30,-8,8,1);

    // variables for checking with the previous position
    double x0_prev = 0; double y0_prev = 0; double z0_prev = 0;
    int same_prev = 1;
    
    double x1_avg = 0; double y1_avg = 0; double z1_avg = 0; double t1_avg = 0; // variables for averages
    vector<double> x1_vec,y1_vec,z1_vec,t1_vec;                                 // variables for data from each point
    
    // calculate averages and standard deviations
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

                // TCanvas* c = new TCanvas();
                // TGraph2D* point_spread = new TGraph2D();
                // string ps_title = "Point spread for x = " + to_string(x0_prev) + ", y = " + to_string(y0_prev) + ", z = " + to_string(z0_prev) + ";x [cm];y [cm];time [ns]";
                // point_spread->SetTitle(ps_title.c_str());
                // for (int i = 0; i < same_prev; i++) point_spread->AddPoint(x1_vec[i],y1_vec[i],t1_vec[i]);
                // point_spread->Draw("P");
                // point_spread->SetMarkerSize(2);
                // point_spread->SetMarkerStyle(2);

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
        vector<TH2F*> v_yx_dx;
        vector<TH2F*> v_yx_dy;
        vector<TGraphErrors*> v_g_yx;
        vector<vector<TArrow*>> vv_g_yx_arrows;
        vector<TH2F*> v_yx_t1;

        TGraphErrors* g_xz = new TGraphErrors();
        vector<TArrow*> v_g_xz_arrows;
        TGraphErrors* g_xt = new TGraphErrors();
        TGraph* g_zt = new TGraph();

        g_xz->SetName("g_xz");
        g_xz->SetTitle("Map of electron displacement (x = 0)");
        g_xz->SetMarkerColor(2);
        g_xz->SetMarkerStyle(6);
        g_xz->GetXaxis()->SetTitle("x [cm]");
        g_xz->GetYaxis()->SetTitle("z [cm]");

        g_xt->SetName("g_xt");
        g_xt->SetTitle("Map of drift times (y = 0)");
        g_xt->SetMarkerColor(2);
        g_xt->SetMarkerStyle(6);
        g_xt->GetXaxis()->SetTitle("x [cm]");
        g_xt->GetYaxis()->SetTitle("t [ns]");

        g_zt->SetName("g_zt");
        g_zt->SetTitle("Original height vs drift time");
        g_zt->SetMarkerColor(2);
        g_zt->SetMarkerStyle(6);
        g_zt->GetXaxis()->SetTitle("z [cm]");
        g_zt->GetYaxis()->SetTitle("t [ns]");

        TH2F* xz_t1 = new TH2F("h_xz_t1","XZ plane t1;x [cm];z [cm]",map.ximax,map.xmin,map.xmax,map.zimax,map.zmin,map.zmax);

        gStyle->SetOptStat(0);

        double z_xy = 7;
        double y_xz = 0;

        // int z_xyi = round(z_xy-map.zmin)/map.step;
        int y_xzi = round(y_xz-map.ymin)/map.step;
        for (int z_xyi = 0; z_xyi <= map.zimax; z_xyi++)
        {
            double z = map.zmin + z_xyi*map.step;
            string h_dx_name = "h_yx_dx_" + to_string(z);
            string h_dy_name = "h_yx_dy_" + to_string(z);
            string h_t1_name = "h_yx_t1_" + to_string(z);
            string g_yx_name = "g_yx_" + to_string(z);
            string h_dx_title = "XY plane dx (z = " + to_string(z) + ");x [cm];y [cm]";
            string h_dy_title = "XY plane dy (z = " + to_string(z) + ");x [cm];y [cm]";
            string h_t1_title = "XY plane t1 (z = " + to_string(z) + ");x [cm];y[cm]";
            string g_yx_title = "Map of electron readout positions, z = " + to_string(z);

            v_yx_dx.push_back(new TH2F(h_dx_name.c_str(),h_dx_title.c_str(),map.yimax,map.ymin,map.ymax,map.ximax,map.xmin,map.xmax));
            v_yx_dy.push_back(new TH2F(h_dy_name.c_str(),h_dy_title.c_str(),map.yimax,map.ymin,map.ymax,map.ximax,map.xmin,map.xmax));
            v_yx_t1.push_back(new TH2F(h_t1_name.c_str(),h_t1_title.c_str(),map.yimax,map.ymin,map.ymax,map.ximax,map.xmin,map.xmax));
            v_g_yx.push_back(new TGraphErrors());

            v_g_yx[z_xyi]->SetName(g_yx_name.c_str());
            v_g_yx[z_xyi]->SetTitle(g_yx_title.c_str());
            v_g_yx[z_xyi]->SetMarkerColor(2);
            v_g_yx[z_xyi]->SetMarkerStyle(6);
            v_g_yx[z_xyi]->GetXaxis()->SetTitle("y [cm]");
            v_g_yx[z_xyi]->GetYaxis()->SetTitle("x [cm]");

            vector<TArrow*> arrows;

            for (int xi = 0; xi <= map.ximax; xi++)
            {
                double x = map.xmin+xi*map.step;

                for (int yi = 0; yi <= map.yimax; yi++)
                {
                    double y = map.ymin+yi*map.step;
                    SensorData& s_curr = map.field[xi][yi][z_xyi];
                    if(s_curr.x1 != 0) v_yx_dx[z_xyi]->Fill(y,x,s_curr.x1-x);
                    if(s_curr.y1 != 0) v_yx_dy[z_xyi]->Fill(y,x,s_curr.y1-y);
                    if(s_curr.t1 != 0) v_yx_t1[z_xyi]->Fill(y,x,s_curr.t1);
                    if((s_curr.x1 != 0) && (s_curr.y1 != 0))
                    {
                        v_g_yx[z_xyi]->AddPoint(s_curr.y1,s_curr.x1);
                        v_g_yx[z_xyi]->SetPointError(v_g_yx[z_xyi]->GetN()-1,s_curr.y1dev,s_curr.x1dev);
                        TArrow* arrow = new TArrow(y,x,s_curr.y1,s_curr.x1);
                        arrows.push_back(arrow);
                    }
                    g_zt->AddPoint(z,s_curr.t1);
                }
            }

            vv_g_yx_arrows.push_back(arrows);
        }

        for (int xi = 0; xi <= map.ximax; xi++)
        {
            double x = map.xmin+xi*map.step;
            for (int zi = 0; zi <= map.zimax; zi++)
            {
                double z = map.zmin+zi*map.step;

                SensorData& s_curr = map.field[xi][y_xzi][zi];

                xz_t1->Fill(x,z,s_curr.t1);
                g_xz->AddPoint(s_curr.x1,z);
                g_xz->SetPointError(g_xz->GetN()-1,s_curr.x1dev,0);
                g_xt->AddPoint(x,s_curr.t1);
                g_xt->SetPointError(g_xt->GetN()-1,s_curr.x1dev,s_curr.t1dev);

                TArrow* arrow = new TArrow(x,z,s_curr.x1,z);
                v_g_xz_arrows.push_back(arrow);
            }
        }
        xz_t1->Write();
        for(TH2F* h : v_yx_dx) h->Write();
        for(TH2F* h : v_yx_dy) h->Write();
        // for(TGraph* g : v_g_xy) g->Write();

        TCanvas* c_g_yx = new TCanvas("c_g_yx","Map of electron readout positions");
        double ymin = -10; double ymax = 10; double xmin = 5; double xmax = 16;
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

        TCanvas* c_g_xz = new TCanvas("c_g_xz","Map of ionization electron displacement");
        g_xz->Draw("AP");
        TLine* lleft  = new TLine(6.51,8,6.51,-8);
        TLine* lright = new TLine(14.61,8,14.61,-8);
        lleft->Draw();lright->Draw();
        for(TArrow* arr : v_g_xz_arrows) {arr->SetArrowSize(0.004); arr->Draw();}
        c_g_xz->Write();

        TCanvas* c_g_xt = new TCanvas("c_g_xt","Map of ionization electron drift times");
        g_xt->Draw("AP");
        lleft  = new TLine(6.51,0,6.51,5000);
        lright = new TLine(14.61,0,14.61,5000);
        lleft->Draw();lright->Draw();
        c_g_xt->Write();

        TCanvas* c_g_zt = new TCanvas("c_g_zt","Original height vs drift time");
        g_zt->Draw("AP");
        c_g_zt->Write();

        outfile->Close();

        TFile* mapfile = new TFile("map.root");
        Field<SensorData>* map = (Field<SensorData>*)mapfile->Get("map");
        TCanvas* c_pads = new TCanvas("c_pads", "Pads distortion for different times.");
        c_pads->Divide(4,4);
        for (int i = 0; i < 16; i++)
        {
            c_pads->cd(i+1);
            TGraph* gr = new TGraph(); gr->AddPoint(-X17::yhigh,X17::xmin);gr->AddPoint(X17::yhigh,X17::xmax);gr->Draw("AP");
            X17::DrawPadsDistortion((i+1)*5000/16,c_pads,map);
            //X17::DrawTrapezoid();
        }
        
        
    }

    return 0;
}

int main(int argc, char const *argv[])
{
    return make_map();
}