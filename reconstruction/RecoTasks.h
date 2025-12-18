#pragma once

// C++ dependencies
#include <string>

// ROOT dependencies
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TPolyLine3D.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

// X17 dependencies
#include "CircleFit2D.h"
#include "CircleFit3D.h"
#include "Field.h"
#include "PadLayout.h"
#include "RK4.h"
#include "TrackLoop.h"
#include "Reconstruction.h"
#include "X17Utilities.h"

// Tasks, that were moved to separate files
#include "RecoPadsTask.h"
#include "RKFitCircleTask.h"
#include "MapRecoCompareTask.h"
#include "MapRecoTask.h"
#include "MicroCircleAndRKFitTask.h"
#include "MicroFitAndSaveTask.h"

using namespace X17;

/// @brief Class for plotting drift time vs distance to readout with linear fit.
class DriftTimeTask : public RecoTask
{
    void PostElectronLoop() override
    {
        std::string c_name = "c_drift" + std::to_string(m_loop->curr_track_index);
        TCanvas* c = new TCanvas(c_name.c_str(),"Drift time");

        TTree* electrons = m_loop->curr_micro_tree;
        electrons->Draw("t1:8-z0","z1>7.0");
        TGraph* tz = new TGraph(electrons->GetSelectedRows(), electrons->GetV2(), electrons->GetV1());
        tz->SetTitle("Drift time as function of distance;distance to readout [cm]; time [ns]");
        tz->SetMarkerStyle(2);
        tz->SetMarkerSize(0.4);
        tz->Draw("ap");
        tz->Fit("pol1","","",8.05,11);

        c->Write();
    }
};

/// @brief Class for plotting XZ of original and reconstructed points of the track.
class XZPlotTask : public RecoTask
{
    TGraph *xz, *xz_reco;

    void PreElectronLoop() override
    {
        xz = new TGraph();
        xz_reco = new TGraph();
    }

    void ElectronLoop() override
    {
        MicroPoint micro = m_loop->curr_micro;
        RecoPoint reco = m_loop->curr_reco;

        xz->AddPoint(micro.x0(),8-micro.z0());
        xz_reco->AddPoint(reco.x(),8-reco.z());

        if (abs(reco.z()-micro.z0()) > 0.5) 
        {
            std::cout << "High z deviation: " << reco.z()-micro.z0() << "\n";
            std::cout << "Simulated coordinates:     x0 = " << micro.x0() << ", y0 = " << micro.y0() << ", z0 = " << micro.z0() << ", t0 = " << micro.t0() << ", e0 = " << micro.e0 << "\n";
            std::cout << "                           x1 = " << micro.x1() << ", y1 = " << micro.y1() << ", z1 = " << micro.z1() << ", t1 = " << micro.t1() << ", e1 = " << micro.e1 << "\n";
            std::cout << "Reconstructed coordinates: x  = " << reco.x()   << ", y  = " << reco.y()   << ", z  = " << reco.z()   << ", count = " << reco.count << "\n\n";
        }
    }

    void PostElectronLoop() override
    {
        std::string c_name = "c_track_xz" + std::to_string(m_loop->curr_track_index);
        TCanvas* c = new TCanvas(c_name.c_str(),"Electron track reconstruction");

        xz_reco->SetTitle("Electron track reconstruction;x [cm]; distance to readout [cm]");
        xz_reco->SetMarkerStyle(2);
        xz_reco->SetMarkerSize(0.4);
        xz_reco->Draw("ap");

        xz->SetMarkerColor(2);
        xz->SetMarkerStyle(7);
        xz->SetMarkerSize(1.2);
        xz->Draw("p same");

        TLegend* leg_xz = new TLegend(0.129,0.786,0.360,0.887);
        leg_xz->AddEntry(xz,"ionization vertices","p");
        leg_xz->AddEntry(xz_reco,"reconstructed","p");
        leg_xz->Draw("same");

        // Fitting both tracks with circles.
        // double min = 0;
        // double max = 15;
        // TF1* circle_fit  = FitCircle2(xz,min,max);
        // TF1* circle_fit2 = FitCircle2(xz_reco,min,max);

        // TGraph* magnetic_x  = new TGraph();
        // TGraph* magnetic_x2 = new TGraph();

        // Field<Vector>* magfield = m_loop->magfield;
        // double step = 0.1;
        // RecoEnergy(circle_fit,*magfield,magnetic_x,min,max,step);
        // RecoEnergy(circle_fit2,*magfield,magnetic_x2,min,max,step);

        c->Write();
    }
};

/// @brief Class for plotting XY of original and reconstructed points of the track.
class XYPlotTask : public RecoTask
{
    TGraph *xy, *xy_reco;

    void PreElectronLoop() override
    {
        xy = new TGraph();
        xy_reco = new TGraph();
    }

    void ElectronLoop() override
    {
        MicroPoint micro = m_loop->curr_micro;
        RecoPoint reco = m_loop->curr_reco;

        xy->AddPoint(micro.x0(),micro.y0());
        xy_reco->AddPoint(reco.x(),reco.y());
    }

    void PostElectronLoop() override
    {
        std::string c_name = "c_track_xy" + std::to_string(m_loop->curr_track_index);
        TCanvas* c = new TCanvas(c_name.c_str(),"Electron track reconstruction");

        xy_reco->SetTitle("Electron track reconstruction;x [cm]; y [cm]");
        xy_reco->SetMarkerStyle(2);
        xy_reco->SetMarkerSize(0.4);
        xy_reco->Draw("ap");

        xy->SetMarkerColor(2);
        xy->SetMarkerStyle(7);
        xy->SetMarkerSize(1.2);
        xy->Draw("p same");

        TLegend* leg_xy = new TLegend(0.129,0.786,0.360,0.887);
        leg_xy->AddEntry(xy,"ionization vertices","p");
        leg_xy->AddEntry(xy_reco,"reconstructed","p");
        leg_xy->Draw("same");

        c->Write();
    }
};

/// @brief Class for plotting residues of reconstructed interaction points.
class GraphResTask : public RecoTask
{
    TGraph *gx_res, *gy_res, *gz_res, *gr_res;

    void PreElectronLoop() override
    {
        gx_res = new TGraph();
        gy_res = new TGraph();
        gz_res = new TGraph();
        gr_res = new TGraph();
    }

    void ElectronLoop() override
    {
        MicroPoint micro = m_loop->curr_micro;
        RecoPoint reco = m_loop->curr_reco;

        gx_res->AddPoint(micro.x0(), reco.x() - micro.x0());
        gy_res->AddPoint(micro.x0(), reco.y() - micro.y0());
        gz_res->AddPoint(micro.x0(), reco.z() - micro.z0());
        gr_res->AddPoint(micro.x0(), std::sqrt(pow((reco.x() - micro.x0()),2) + pow((reco.y() - micro.y0()),2) + pow((reco.z()-micro.z0()),2)));
    }

    void PostElectronLoop() override
    {
        std::string c_name = "c_fit_res" + std::to_string(m_loop->curr_track_index);
        TCanvas* c = new TCanvas(c_name.c_str(),"Reconstruction residuals");

        c->Divide(2,2);
        
        c->cd(1);
        gx_res->SetTitle("X residuals;x [cm];#Deltax [cm]");
        gx_res->SetMarkerStyle(2);
        gx_res->SetMarkerSize(0.4);
        gx_res->Draw("ap");
        
        c->cd(2);
        gy_res->SetTitle("Y residuals;x [cm];#Deltay [cm]");
        gy_res->SetMarkerStyle(2);
        gy_res->SetMarkerSize(0.4);
        gy_res->Draw("ap");
        
        c->cd(3);
        gz_res->SetTitle("Z residuals;x [cm];#Deltaz [cm]");
        gz_res->SetMarkerStyle(2);
        gz_res->SetMarkerSize(0.4);
        gz_res->Draw("ap");
        
        c->cd(4);
        gr_res->SetTitle("Residuals;x [cm];distance [cm]");
        gr_res->SetMarkerStyle(2);
        gr_res->SetMarkerSize(0.4);
        gr_res->Draw("ap");

        c->Write();
    }
};

/// @brief Class for plotting residues of reconstructed interaction points.
class HistResTask : public RecoTask
{
    TH1F *hx_res, *hy_res, *hz_res, *hr_res;

    void PreElectronLoop() override
    {
        hx_res = new TH1F("hx_res","X residuals;x deviation [cm];# of electrons",25,-0.2,0.2);
        hy_res = new TH1F("hy_res","Y residuals;y deviation [cm];# of electrons",25,-0.2,0.2);
        hz_res = new TH1F("hz_res","Z residuals;z deviation [cm];# of electrons",25,-0.2,0.2);
        hr_res = new TH1F("hr_res","Residuals;Deviation [cm];# of electrons",25,0,0.25);
    }

    void ElectronLoop() override
    {
        MicroPoint micro = m_loop->curr_micro;
        RecoPoint reco = m_loop->curr_reco;

        hx_res->Fill(reco.x() - micro.x0());
        hy_res->Fill(reco.y() - micro.y0());
        hz_res->Fill(reco.z() - micro.z0());
        hr_res->Fill(std::sqrt(pow((reco.x() - micro.x0()),2) + pow((reco.y() - micro.y0()),2) + pow((reco.z() - micro.z0()),2)));
    }

    void PostElectronLoop() override
    {
        std::string c_name = "c_fit_res2" + std::to_string(m_loop->curr_track_index);
        TCanvas* c = new TCanvas(c_name.c_str(),"Reconstruction residuals");
        c->Divide(2,2);
        c->cd(1); hx_res->Draw();
        c->cd(2); hy_res->Draw();
        c->cd(3); hz_res->Draw();
        c->cd(4); hr_res->Draw();

        c->Write();   
    }
};

/// @brief Class for plotting Runge-Kutta simulated tracks with lower than selected energy.
class PlotSelectionTask : public RecoTask
{
    std::vector<TGraph2D*> tracks;
    double E_max;
    // CircleFit3D* cfit = nullptr;
    RKFitCircleTask* cfit_energy = nullptr;

    void PreTrackLoop() override
    {
        // cfit = new CircleFit3D();
        // cfit->SetFitter();
    }

    void PreElectronLoop() override
    {
        const TrackRK* track = m_loop->curr_rk;

        // cfit = new CircleFit3D(track->origin,track->orientation);
        tracks.push_back(new TGraph2D());
    }

    void ElectronLoop() override
    {
        const TrackRK* track = m_loop->curr_rk;
        RKPoint p = m_loop->curr_rkpoint;

        // cfit->AddPoint(p);
        tracks.back()->AddPoint(p.x(),p.y(),p.z());
    }

    void PostElectronLoop() override
    {   
        // cfit->FitCircle3D();

        if (cfit_energy->cfit->GetEnergy(*(m_loop->magfield)) > E_max) tracks.pop_back(); 
    }

    void PostTrackLoop() override
    {
        using namespace X17::constants;
        TCanvas* c = new TCanvas("c_cfit_failed","Tracks with failed circle fit");

        // A histogram for scalling of the axes.
        TH3F* scale = new TH3F("scale","Tracks with failed circle fit;x [cm];y [cm];z [cm]",1,xmin,xmax,1,-yhigh,yhigh,1,zmin,zmax);
        scale->Draw("");
        gStyle->SetOptStat(0);
        scale->GetXaxis()->SetTitleOffset(1.5);
        scale->GetZaxis()->SetTitleOffset(1.5);

        for (TGraph2D* g : tracks) 
        {
            g->SetLineColor(kRed);
            g->Draw("LINE same");
        }

        c->Write();

        std::cout << "PlotSelectionTask: failed fits: " << tracks.size() << "\n";
    }

public:
    PlotSelectionTask(RKFitCircleTask* cfit, double E_max = 3.5e+6) : E_max(E_max), cfit_energy(cfit) {}
};

/// @brief Task for plotting several Runge-Kutta tracks with a color palette.
class PlotForwardTask : public RecoTask
{
    std::vector<TLine*> track_segments;
    double x_prev = X17::constants::xmin;
    double z_prev = 0;
    int color;

    void PreElectronLoop() override
    {
        double norm_energy = (m_loop->curr_rk->kin_energy - 3e+6)/10e+6;

        int color_index = norm_energy * (gStyle->GetNumberOfColors() - 1);
        color = gStyle->GetColorPalette(color_index);
    }

    void ElectronLoop() override
    {
        const TrackRK* track = m_loop->curr_rk;
        RKPoint p = m_loop->curr_rkpoint;
        TLine* line = new TLine(x_prev,z_prev,p.x(),p.z());
        line->SetLineColor(color);
        x_prev = p.x();
        z_prev = p.z();
        if (p.x() < 15) track_segments.push_back(line);
    }

    void PostElectronLoop() override
    {
        x_prev = X17::constants::xmin;
        z_prev = 0;
    }
    
    void PostTrackLoop() override
    {
        using namespace X17::constants;

        TCanvas* c = new TCanvas("c_forward","Forward tracks with different energies",1.5*500/0.72,1.5*500/7*8.5/0.8);
        c->SetTopMargin(0.07);
        c->SetBottomMargin(0.13);
        c->SetRightMargin(0.2);
        TH2F* scale = new TH2F("scale",";x [cm];z [cm];E [MeV]",10,xmin,15,10,-6,1);
        scale->SetStats(0);
        scale->SetBinContent(1,3);
        scale->SetMinimum(3);
        scale->SetMaximum(13);
        scale->GetXaxis()->SetLabelSize(0.04);
        scale->GetXaxis()->SetTitleSize(0.05);
        scale->GetYaxis()->SetLabelSize(0.04);
        scale->GetYaxis()->SetTitleSize(0.05);
        scale->GetZaxis()->SetLabelSize(0.04);
        scale->GetZaxis()->SetTitleSize(0.05);
        scale->GetXaxis()->SetTitleOffset(0.9);
        scale->GetYaxis()->SetTitleOffset(0.7); 
        scale->GetZaxis()->SetTitleOffset(1.1);
        scale->Draw("colz");

        for (TLine* l : track_segments)
        {
            l->SetLineWidth(2);
            l->Draw("same");
        }
        c->Write();
    }
};