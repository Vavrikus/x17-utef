#pragma once

#include "TGraph.h"
#include "TGraph2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLegend.h"
#include "TPolyLine3D.h"
#include "TStyle.h"
#include "TTree.h"

#include "CircleFit2D.h"
#include "CircleFit3D.h"
#include "RK4.h"
#include "TrackLoop.h"
#include "../X17Utilities.h"

// File included at the end of TrackLoop.h to avoid circular dependencies and extra translation units

/// @brief Class for plotting drift time vs distance to readout with linear fit
class DriftTimeTask : public PlotTask
{
    std::string canvasName  = "c_drift";
    std::string canvasTitle = "Drift time";

    void PreLoop() {}
    void Loop() {}

    void PostLoop()
    {
        TTree* electrons = loop->GetSingleTrack();
        electrons->Draw("t1:8-z0","z1>7.0");
        TGraph* tz = new TGraph(electrons->GetSelectedRows(), electrons->GetV2(), electrons->GetV1());
        tz->SetTitle("Drift time as function of distance;distance to readout [cm]; time [ns]");
        tz->SetMarkerStyle(2);
        tz->SetMarkerSize(0.4);
        tz->Draw("ap");
        tz->Fit("pol1","","",8.05,11);

        // // getting fit parameters
        // TF1* tz_fit = (TF1*) tz->GetListOfFunctions()->FindObject("pol1");
        // double a0 = tz_fit->GetParameter(0);
        // double a1 = tz_fit->GetParameter(1);
        // double b0 = -a0/a1; //inverse polynomial param
        // double b1 = 1.0/a1; //inverse polynomial param
    }
};

/// @brief Class for plotting XZ of original and reconstructed points of the track
class XZPlotTask : public PlotTask
{
    std::string canvasName  = "c_track_xz";
    std::string canvasTitle = "Electron track reconstruction";

    TGraph *xz, *xz_reco;

    void PreLoop() {xz = new TGraph(); xz_reco = new TGraph();}

    void Loop() 
    {
        IonElectron e = loop->GetCurrentElectron();
        SensorData reco = loop->GetCurrentReco();

        xz->AddPoint(e.x0,8-e.z0);
        xz_reco->AddPoint(reco.x1,8-reco.z1);

        // SensorData reco = RecoPoint(map,x1,z1,t1,0.001);
        // xz_reco->AddPoint(z1,b0+b1*t1);
    }

    void PostLoop()
    {
        xz_reco->SetTitle("Electron track reconstruction;x [cm]; distance to readout [cm]");
        xz_reco->SetMarkerStyle(2);
        xz_reco->SetMarkerSize(0.4);
        xz_reco->Draw("ap");

        xz->SetMarkerColor(2);
        xz->SetMarkerStyle(7);
        xz->SetMarkerSize(1.2);
        xz->Draw("p same");

        TLegend* leg_xz = new TLegend(0.129,0.786,0.360,0.887);
        leg_xz->AddEntry(xz,"original");
        leg_xz->AddEntry(xz_reco,"reconstructed");
        leg_xz->Draw("same");

        // fitting both tracks with circles
        double min = 0;  // 7
        double max = 15; // 10.5
        TF1* circle_fit  = FitCircle2(xz,min,max);
        TF1* circle_fit2 = FitCircle2(xz_reco,min,max);

        TGraph* magnetic_x  = new TGraph();
        TGraph* magnetic_x2 = new TGraph();

        VectorField* magfield = loop->GetMagField();
        double step = 0.1;
        RecoEnergy(circle_fit,magfield,magnetic_x,min,max,step);
        RecoEnergy(circle_fit2,magfield,magnetic_x2,min,max,step);
    }
};

/// @brief Class for plotting XY of original and reconstructed points of the track
class XYPlotTask : public PlotTask
{
    std::string canvasName  = "c_track_xy";
    std::string canvasTitle = "Electron track reconstruction";

    TGraph *xy, *xy_reco;

    void PreLoop() {xy = new TGraph(); xy_reco = new TGraph();}

    void Loop() 
    {
        IonElectron e = loop->GetCurrentElectron();
        SensorData reco = loop->GetCurrentReco();

        xy->AddPoint(e.x0,e.y0);
        xy_reco->AddPoint(reco.x1,reco.y1);
    }

    void PostLoop()
    {    
        xy_reco->SetTitle("Electron track reconstruction;x [cm]; y [cm]");
        xy_reco->SetMarkerStyle(2);
        xy_reco->SetMarkerSize(0.4);
        xy_reco->Draw("ap");

        xy->SetMarkerColor(2);
        xy->SetMarkerStyle(7);
        xy->SetMarkerSize(1.2);
        xy->Draw("p same");

        TLegend* leg_xy = new TLegend(0.129,0.786,0.360,0.887);
        leg_xy->AddEntry(xy,"original");
        leg_xy->AddEntry(xy_reco,"reconstructed");
        leg_xy->Draw("same");
    }
};

/// @brief Class for plotting residues of reconstructed interaction points
class GraphResTask : public PlotTask
{
    std::string canvasName  = "c_fit_res";
    std::string canvasTitle = "Reconstruction residuals";

    TGraph *gx_res, *gy_res, *gz_res, *gr_res;

    void PreLoop() 
    {
        gx_res = new TGraph();
        gy_res = new TGraph();
        gz_res = new TGraph();
        gr_res = new TGraph();
    }

    void Loop() 
    {
        IonElectron e = loop->GetCurrentElectron();
        SensorData reco = loop->GetCurrentReco();

        gx_res->AddPoint(e.x0,reco.x1-e.x0);
        gy_res->AddPoint(e.x0,reco.y1-e.y0);
        gz_res->AddPoint(e.x0,reco.z1-e.z0);
        gr_res->AddPoint(e.x0,sqrt(pow((reco.x1-e.x0),2)+pow((reco.y1-e.y0),2)+pow((reco.z1-e.z0),2)));
    }

    void PostLoop()
    {
        canvas->Divide(2,2);
        
        canvas->cd(1);
        gx_res->SetTitle("X residuals;x [cm];Δx [cm]");
        gx_res->SetMarkerStyle(2);
        gx_res->SetMarkerSize(0.4);
        gx_res->Draw("ap");
        
        canvas->cd(2);
        gy_res->SetTitle("Y residuals;x [cm];Δy [cm]");
        gy_res->SetMarkerStyle(2);
        gy_res->SetMarkerSize(0.4);
        gy_res->Draw("ap");
        
        canvas->cd(3);
        gz_res->SetTitle("Z residuals;x [cm];Δz [cm]");
        gz_res->SetMarkerStyle(2);
        gz_res->SetMarkerSize(0.4);
        gz_res->Draw("ap");
        
        canvas->cd(4);
        gr_res->SetTitle("Residuals;x [cm];distance [cm]");
        gr_res->SetMarkerStyle(2);
        gr_res->SetMarkerSize(0.4);
        gr_res->Draw("ap");       
    }
};

/// @brief Class for plotting residues of reconstructed interaction points
class HistResTask : public PlotTask
{
    std::string canvasName  = "c_fit_res2";
    std::string canvasTitle = "Reconstruction residuals";

    TH1F *hx_res, *hy_res, *hz_res, *hr_res;

    void PreLoop() 
    {
        hx_res = new TH1F("hx_res","X residuals",25,-0.5,0.5);
        hy_res = new TH1F("hy_res","Y residuals",25,-0.5,0.5);
        hz_res = new TH1F("hz_res","Z residuals",25,-0.5,0.5);
        hr_res = new TH1F("hr_res","Residuals",25,0,1);
    }

    void Loop() 
    {
        IonElectron e = loop->GetCurrentElectron();
        SensorData reco = loop->GetCurrentReco();

        hx_res->Fill(reco.x1-e.x0);
        hy_res->Fill(reco.y1-e.y0);
        hz_res->Fill(reco.z1-e.z0);
        hr_res->Fill(sqrt(pow((reco.x1-e.x0),2)+pow((reco.y1-e.y0),2)+pow((reco.z1-e.z0),2)));        
    }

    void PostLoop()
    {
        canvas->Divide(2,2);
        canvas->cd(1); hx_res->Draw();
        canvas->cd(2); hy_res->Draw();
        canvas->cd(3); hz_res->Draw();
        canvas->cd(4); hr_res->Draw();       
    }
};

/// @brief Class for padded reconstruction
class RecoPadsTask : public PlotTask
{
    friend class CircleAndRKFitTask;

    std::string canvasName  = "c_track_xyz";
    std::string canvasTitle = "Electron track reconstruction with pads and time bins";

    static constexpr int timebins = 50;
    int padhits[X17::channels][timebins];
    TGraph2D *g_xyz, *g_xyz_reco;

    double height = -2.5;

    void PreLoop() 
    {
        for (int i = 0; i < X17::channels; i++) for (int j = 0; j < timebins; j++) {padhits[i][j] = 0;}

        g_xyz      = new TGraph2D();
        g_xyz_reco = new TGraph2D();
    }

    void Loop() 
    {
        IonElectron e = loop->GetCurrentElectron();

        int channel = X17::GetPad(e.x1,e.y1);
        int timebin = e.t1/100.0;

        if(timebin > timebins - 1) cerr << "ERROR: Invalid timebin: " << timebin << endl;
        if(channel == -1) cerr << "ERROR: No pad hit found. Coordinates x,y: " << e.x1 << ", " << e.y1 << endl;

        padhits[channel-1][timebin]++;
        g_xyz->AddPoint(e.x0,e.y0,e.z0);
    }

    void PostLoop()
    {
        // Reconstruction with pads
        for (int i = 0; i < X17::channels; i++)
        {
            for (int j = 0; j < timebins; j++)
            {
                if(padhits[i][j] != 0)
                {
                    double time = 100 * j + 50;
                    double xpad,ypad;
                    X17::GetPadCenter(i+1,xpad,ypad);

                    SensorData reco = loop->GetMap()->Invert(xpad,ypad,time);
                    g_xyz_reco->AddPoint(reco.x1,reco.y1,reco.z1);
                }
            }        
        }
        
        // histogram for scalling axes
        TH3F* scale = new TH3F("scale","Electron track reconstruction;x [cm]; y [cm];z [cm]",1,X17::xmin,X17::xmax,1,-X17::yhigh,X17::yhigh,1,height,0);
        scale->Draw("");
        gStyle->SetOptStat(0);
        scale->GetXaxis()->SetTitleOffset(1.5);
        scale->GetZaxis()->SetTitleOffset(1.5);
        
        g_xyz_reco->SetTitle("Electron track reconstruction;x [cm]; y [cm];z [cm]");
        g_xyz_reco->SetMarkerStyle(2);
        g_xyz_reco->SetMarkerSize(0.4);
        g_xyz_reco->Draw("p same");

        g_xyz->SetMarkerColor(2);
        g_xyz->SetMarkerStyle(7);
        g_xyz->SetMarkerSize(1.2);
        g_xyz->Draw("p same");

        X17::DrawPads3D(height);
    }

public:
    RecoPadsTask(double pad_height = -2.5) : height(pad_height) {}
};

/// @brief Class for plotting circle and Runge-Kutta track fits
class CircleAndRKFitTask : public PlotTask
{
    CircleFit3D& cfit3d = CircleFit3D::NewCircleFit({X17::xmin,0,0},{1,0,0});
    RecoPadsTask* reco_task = nullptr;

    void PreLoop() {}

    void Loop()
    {
        IonElectron e = loop->GetCurrentElectron();
        cfit3d.AddPoint(e.x0,e.y0,e.z0,1);
    }

    void PostLoop()
    {
        // cfit3d.Prefit();
        cfit3d.SetFitter(4);
        cfit3d.FitCircle3D();
        cfit3d.PrintFitParams();

        VectorField* magfield = loop->GetMagField();

        TGraph2D* g_cfit3d = cfit3d.GetGraph(0.1,0);
        TGraph* g_en = cfit3d.GetEnergyGraph(magfield);

        g_cfit3d->SetLineColor(kGreen);
        g_cfit3d->SetLineWidth(2);
        g_cfit3d->Draw("LINE same");

        cout << "Reconstructed energy: " << cfit3d.GetEnergy(magfield,true) << " (middle field), " << cfit3d.GetEnergy(magfield,false) << " (average field).\n";

        vector<DataPoint> reco_data;
        CircleFit3D& cfit3d_reco = CircleFit3D::NewCircleFit({0,0,0},{1,0,0});

        for (int i = 0; i < X17::channels; i++)
        {
            for (int j = 0; j < reco_task->timebins; j++)
            {
                if(reco_task->padhits[i][j] != 0)
                {
                    double time = 100 * j + 50;
                    double xpad,ypad;
                    X17::GetPadCenter(i+1,xpad,ypad);

                    SensorData reco = loop->GetMap()->Invert(xpad,ypad,time);
                    cfit3d_reco.AddPoint(reco.x1,reco.y1,reco.z1,reco_task->padhits[i][j]);
                    reco_data.emplace_back(reco.x1,reco.y1,reco.z1,reco_task->padhits[i][j]);
                }
            }        
        }

        auto reco_markers = GetDataMarkers(reco_data);
        for (auto m : reco_markers) m->Draw("same");


        // cfit3d_reco.Prefit();
        cfit3d_reco.SetFitter(4);
        cfit3d_reco.FitCircle3D();
        cfit3d_reco.PrintFitParams();

        TGraph2D* g_cfit3d_reco = cfit3d_reco.GetGraph(0.1,0);
        TGraph* g_en_reco = cfit3d_reco.GetEnergyGraph(magfield);
        g_cfit3d_reco->SetLineColor(kBlue);
        g_cfit3d_reco->SetLineWidth(2);
        g_cfit3d_reco->Draw("LINE same");

        cout << "Reconstructed energy: " << cfit3d_reco.GetEnergy(magfield,true) << " (middle field), " << cfit3d_reco.GetEnergy(magfield,false) << " (average field).\n";

        RKFit* rkfit = new RKFit(magfield,true,1E-13,{0,0,0},{1,0,0},cfit3d_reco.GetData());
        rkfit->SetEnergy(cfit3d_reco.GetEnergy(magfield,true));
        rkfit->SetFitter();
        rkfit->FitRK();
        rkfit->PrintFitParams();
        TGraph2D* g_rkfit = rkfit->GetFitGraph();
        g_rkfit->SetLineColor(kMagenta);
        g_rkfit->SetLineWidth(2);
        g_rkfit->Draw("LINE same");

        TLegend* leg_xyz = new TLegend(0.741,0.742,0.956,0.931);
        leg_xyz->AddEntry(reco_task->g_xyz,"original");
        leg_xyz->AddEntry(g_cfit3d,"fit");
        leg_xyz->AddEntry(reco_task->g_xyz_reco,"reconstructed");
        leg_xyz->AddEntry(g_cfit3d_reco,"fit");
        leg_xyz->AddEntry(g_rkfit,"Runge-Kutta fit");
        leg_xyz->Draw("same");

        TCanvas* c_fit3d = new TCanvas("c_fit3d","Fit 3D");
            g_cfit3d_reco->SetTitle("Fit 3D;x [cm]; y [cm];z [cm]");
            g_cfit3d_reco->Draw("LINE");
            g_cfit3d->Draw("LINE same");
            g_rkfit->Draw("LINE same");
            g_cfit3d_reco->GetYaxis()->SetNdivisions(508);
            g_cfit3d_reco->GetXaxis()->SetTitleOffset(1.5);
            g_cfit3d_reco->GetZaxis()->SetTitleOffset(1.5);

        TCanvas* c_energy = new TCanvas("c_energy","Energy along track");
            g_en->SetMarkerColor(kRed);
            g_en->SetMarkerSize(1);
            g_en->SetMarkerStyle(2);
            g_en->Draw("AP");
            g_en_reco->SetMarkerColor(kBlue);
            g_en_reco->SetMarkerSize(1);
            g_en_reco->SetMarkerStyle(2);
            g_en_reco->Draw("P same");

            TLegend* leg_en = new TLegend(0.129,0.786,0.360,0.887);
            leg_en->AddEntry(g_en,"fit original");
            leg_en->AddEntry(g_en_reco,"fit reconstructed");
            leg_en->Draw("same");
    }

public:
    CircleAndRKFitTask(RecoPadsTask* t) : reco_task(t) {makeNewCanvas = false;}
};

/// @brief Class for plotting Runge-Kutta simulated energy vs circle fit reconstruction
class CircleFitEnergyTask : public PlotTask
{
    std::string canvasName  = "c_cfit_energy";
    std::string canvasTitle = "Original vs reconstructed energy";

    TH2F* h_energy;
    std::vector<TH1F*> h_energy_diff;
    int bins;

    std::vector<TGraph2D*> tracks,fits;
    double Ediff_max = 2e+6;

    void PreLoop()
    {
        const char* h_energy_title = "Original RK vs circle fit energy;RK energy [MeV];circle fit energy [MeV]";
        h_energy = new TH2F("h_energy",h_energy_title,bins,4,12,bins,0,15);

        for (int i = 4; i < 12; i++)
        {
            std::string name  = "h_energy_diff_" + to_string(i);
            std::string title = "Difference between measured and fitted energy;Energy difference [MeV];Count";
            h_energy_diff.push_back(new TH1F(name.c_str(), title.c_str(),25,-2+i+0.5,2+i+0.5));
        }

        CircleFit3D::GetCircleFit().SetFitter(4,false);
    }

    void Loop()
    {
        const TrackRK* track = loop->GetCurrentTrack();

        CircleFit3D& cfit = CircleFit3D::NewCircleFit(track->origin,track->orientation);
        tracks.push_back(new TGraph2D());
        cfit.SetAlpha(track->electron);

        for (auto p : track->points) {cfit.AddPoint(p);tracks.back()->AddPoint(p.x,p.y,p.z);}//g_2d->AddPoint(p.x,p.y,p.z);}
        
        // cfit.Prefit();
        cfit.FitCircle3D();
        // cfit.PrintFitParams();

        // g_fit = cfit.GetGraph();

        h_energy->Fill(track->kin_energy/1e+6,cfit.GetEnergy(loop->GetMagField())/1e+6);
        h_energy_diff[floor(track->kin_energy/1e+6)-4]->Fill((track->kin_energy+cfit.GetEnergy(loop->GetMagField())-track->kin_energy)/1e+6);
        if (abs(cfit.GetEnergy(loop->GetMagField())-track->kin_energy) > Ediff_max)
        {
            cout << "\n Electon track? --> " << track->electron << "\n";
            cfit.PrintFitParams();
            fits.push_back(cfit.GetGraph());
        }
        if (abs(cfit.GetEnergy(loop->GetMagField())-track->kin_energy) < Ediff_max) tracks.pop_back();         
    }

    void PostLoop()
    {
        h_energy->Draw("colz");
        TF1 *f1 = new TF1("f1","x",4,12);
        f1->Draw("same");
        // TCanvas* c = new TCanvas();
        // g_2d->Draw("LINE");
        // g_2d->SetLineColor(kBlue);
        // g_2d->SetTitle("RK track (blue) and its fit (red);x [cm];y [cm];z [cm]");
        // g_fit->Draw("LINE same");
        // g_fit->SetLineColor(kRed);
        TCanvas* c0 = new TCanvas("c_energy_diff","Energy difference");
        TH1F* h1 = c0->DrawFrame(-2,0,14,50);
        for(TH1F* h : h_energy_diff) h->Draw("same plc hist");
        gStyle->SetOptStat(1111);

        TCanvas* c = new TCanvas("c_cfit_failed","Tracks with failed circle fit");
        // histogram for scalling axes
        TH3F* scale = new TH3F("scale","Tracks with failed circle fit;x [cm];y [cm];z [cm]",1,X17::xmin,X17::xmax,1,-X17::yhigh,X17::yhigh,1,X17::zmin,X17::zmax);
        scale->Draw("");
        gStyle->SetOptStat(0);
        scale->GetXaxis()->SetTitleOffset(1.5);
        scale->GetZaxis()->SetTitleOffset(1.5);

        for (TGraph2D* g : tracks) 
        {
            g->SetLineColor(kRed);
            g->Draw("LINE same");
        }

        for (TGraph2D* g : fits) 
        {
            g->SetLineColor(kBlue);
            g->Draw("LINE same");
        }

        cout << "n_failed_fits: " << tracks.size() << "\n";
    }

public:
    CircleFitEnergyTask(const int& bins = 100) : bins(bins) {}
};


/// @brief Class for plotting Runge-Kutta simulated tracks with lower than selected energy
class PlotSelectionTask : public PlotTask
{
    std::string canvasName  = "c_cfit_failed";
    std::string canvasTitle = "Tracks with failed circle fit";

    std::vector<TGraph2D*> tracks;
    double E_max;

    void PreLoop()
    {
        CircleFit3D::GetCircleFit().SetFitter();
    }

    void Loop()
    {
        const TrackRK* track = loop->GetCurrentTrack();

        CircleFit3D& cfit = CircleFit3D::NewCircleFit(track->origin,track->orientation);
        tracks.push_back(new TGraph2D());

        for (auto p : track->points) {cfit.AddPoint(p);tracks.back()->AddPoint(p.x,p.y,p.z);}
        
        cfit.FitCircle3D();

        if (cfit.GetEnergy(loop->GetMagField()) > E_max) tracks.pop_back(); 
    }

    void PostLoop()
    {
        // histogram for scalling axes
        TH3F* scale = new TH3F("scale","Tracks with failed circle fit;x [cm];y [cm];z [cm]",1,X17::xmin,X17::xmax,1,-X17::yhigh,X17::yhigh,1,X17::zmin,X17::zmax);
        scale->Draw("");
        gStyle->SetOptStat(0);
        scale->GetXaxis()->SetTitleOffset(1.5);
        scale->GetZaxis()->SetTitleOffset(1.5);

        for (TGraph2D* g : tracks) 
        {
            g->SetLineColor(kRed);
            g->Draw("LINE same");
        }

        cout << "n_failed_fits: " << tracks.size() << "\n";
    }

public:
    PlotSelectionTask(double E_max = 3.5e+6) : E_max(E_max) {}
};