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
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"
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

        if(abs(reco.z()-micro.z0()) > 0.5) 
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
        gr_res->AddPoint(micro.x0(), sqrt(pow((reco.x() - micro.x0()),2) + pow((reco.y() - micro.y0()),2) + pow((reco.z()-micro.z0()),2)));
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
        hr_res->Fill(sqrt(pow((reco.x() - micro.x0()),2) + pow((reco.y() - micro.y0()),2) + pow((reco.z() - micro.z0()),2)));
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

/// @brief Class for padded reconstruction.
class RecoPadsTask : public RecoTask
{
    friend class CircleAndRKFitTask;
    friend class FitAndSaveTask;

    static constexpr int timebins = 200;
    int padhits[constants::channels][timebins];
    TGraph2D *g_xyz, *g_xyz_reco;
    TCanvas* c_reco;

    double height; // Height at which the pads will be drawn.

    void PreElectronLoop() override
    {
        for (int i = 0; i < constants::channels; i++) for (int j = 0; j < timebins; j++) {padhits[i][j] = 0;}

        if (m_loop->make_track_plots)
        {
            g_xyz      = new TGraph2D();
            g_xyz_reco = new TGraph2D();
        }
    }

    void ElectronLoop() override
    {
        MicroPoint micro = m_loop->curr_micro;

        int channel = DefaultLayout::GetDefaultLayout().GetPad(micro.x1(),micro.y1());
        int timebin = micro.t1() / 100.0;

        if (timebin > timebins - 1) std::cerr << "ERROR: Invalid timebin: " << timebin << std::endl;
        if (channel == -1) std::cerr << "ERROR: No pad hit found. Coordinates x,y: " << micro.x1() << ", " << micro.y1() << std::endl;

        padhits[channel-1][timebin]++;
        
        if (m_loop->make_track_plots) g_xyz->AddPoint(micro.x0(),micro.y0(),micro.z0());
    }

    void PostElectronLoop() override
    {
        using namespace X17::constants;

        // Reconstruction with pads.
        for (int i = 0; i < channels; i++)
        {
            for (int j = 0; j < timebins; j++)
            {
                if(padhits[i][j] != 0)
                {
                    double time = 100 * j + 50;
                    double xpad,ypad;
                    DefaultLayout::GetDefaultLayout().GetPadCenter(i+1,xpad,ypad);

                    RecoPoint reco = Reconstruct(*(m_loop->map),EndPoint(xpad,ypad,zmax,time));
                    if (m_loop->make_track_plots) g_xyz_reco->AddPoint(reco.x(),reco.y(),reco.z());
                }
            }        
        }
        
        if (m_loop->make_track_plots)
        {
            std::string c_reco_name = "c_track_xyz" + std::to_string(m_loop->curr_track_index);
            c_reco = new TCanvas(c_reco_name.c_str(),"Electron track reconstruction with pads and time bins");

            // A histogram for scalling of the axes.
            TH3F* scale = new TH3F("scale","Electron track reconstruction;x [cm]; y [cm];z [cm]",1,xmin,xmax,1,-2,2,1,height,-height);
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

            DefaultLayout::GetDefaultLayout().DrawPads3D(height);
        }
    }

public:
    RecoPadsTask(double pad_height = -2.5) : height(pad_height) {}
};

/// @brief Class for plotting circle and Runge-Kutta track fits.
class CircleAndRKFitTask : public RecoTask
{
    CircleFit3D* cfit3d = nullptr;
    RecoPadsTask* reco_task = nullptr;

    TH3F* h_all;
    TH2F* h_theta_phi;
    TH2F* h_theta_energy;
    TH2F* h_phi_energy;
    TH2F* h_deltaenergy_energy;
    TH1F* h_delta_energy;

    void PreTrackLoop() override
    {
        // Binning.
        int angle_bins  = 21;
        int energy_bins = 11;

        // Ranges for simulation.
        double theta_max = (180/TMath::Pi())*atan((constants::win_height/2)/constants::xmin); // The maximal simulated theta [rad].
        double theta_min = -theta_max;                                                        // The minimal simulated theta [rad].
        double phi_max = (180/TMath::Pi())*atan((constants::win_width/2)/constants::xmin);    // The maximal simulated phi [rad].
        double phi_min = -phi_max;                                                            // The minimal simulated phi [rad].
        double E_max = 13;                                                                    // The maximal simulated energy [MeV].
        double E_min = 3;                                                                     // THe minimal simulated energy [MeV].

        // Adjusting boundaries for the extra bin.
        phi_max   = phi_min   + (phi_max  - phi_min)   * angle_bins  / (angle_bins  - 1);
        theta_max = theta_min + (theta_max- theta_min) * angle_bins  / (angle_bins  - 1);
        E_max     = E_min     + (E_max    - E_min)     * energy_bins / (energy_bins - 1);

        // Histogram titles.
        const char* h_all_title          = "Energy resolution (#DeltaE/E);Phi [deg];Theta [deg];Simulated energy [MeV];#DeltaE/E [\%]";
        const char* h_theta_phi_title    = "Energy resolution (#DeltaE/E);Phi [deg];Theta [deg];#DeltaE/E [\%]";
        const char* h_theta_energy_title = "Energy resolution (#DeltaE/E);Theta [deg];Simulated energy [MeV];#DeltaE/E [\%]";
        const char* h_phi_energy_title   = "Energy resolution (#DeltaE/E);Phi [deg];Simulated energy [MeV];#DeltaE/E [\%]";

        // Histogram initializations.
        h_all          = new TH3F("h_all",h_all_title,angle_bins,phi_min,phi_max,angle_bins,theta_min,theta_max,energy_bins,E_min,E_max);
        h_theta_phi    = new TH2F("h_theta_phi",h_theta_phi_title,angle_bins,phi_min,phi_max,angle_bins,theta_min,theta_max);
        h_theta_energy = new TH2F("h_theta_energy",h_theta_energy_title,angle_bins,theta_min,theta_max,energy_bins,E_min,E_max);
        h_phi_energy   = new TH2F("h_phi_energy",h_phi_energy_title,angle_bins,phi_min,phi_max,energy_bins,E_min,E_max);
        h_deltaenergy_energy = new TH2F("h_deltaenergy_energy","",energy_bins,E_min,E_max,100,-50,50);
    }

    void PreElectronLoop() override
    {
        if(m_loop->curr_loop == TrackLoop::SINGLE) cfit3d = new CircleFit3D({constants::xmin,0,0},{1,0,0});
        else if(m_loop->curr_loop == TrackLoop::MULTI) cfit3d = new CircleFit3D(m_loop->curr_microtrack->origin,m_loop->curr_microtrack->orientation);
    }

    void ElectronLoop() override
    {
        MicroPoint micro = m_loop->curr_micro;
        cfit3d->AddPoint(micro.x0(),micro.y0(),micro.z0(),1);
    }

    void PostElectronLoop() override
    {
        using namespace X17::constants;

        if(m_loop->curr_loop == TrackLoop::MULTI) std::cout << "Simulated energy: " << m_loop->curr_microtrack->kin_energy << " eV.";

        cfit3d->SetFitter(4,false);
        if(m_loop->curr_loop == TrackLoop::MULTI) cfit3d->SetParameters(0,10,1.5,m_loop->curr_microtrack->electron);
        cfit3d->FitCircle3D();
        cfit3d->PrintFitParams();

        Field<Vector>* magfield = m_loop->magfield;

        std::cout << "Reconstructed energy (circle fit no pads): " << cfit3d->GetEnergy(*magfield,true) << " (middle field), " << cfit3d->GetEnergy(*magfield,false) << " (average field).\n";

        std::vector<RecoPoint> reco_data;
        CircleFit3D* cfit3d_reco;
        if(m_loop->curr_loop == TrackLoop::SINGLE) cfit3d_reco = new CircleFit3D({0,0,0},{1,0,0});
        else if(m_loop->curr_loop == TrackLoop::MULTI) cfit3d_reco = new CircleFit3D(m_loop->curr_microtrack->origin,m_loop->curr_microtrack->orientation);

        for (int i = 0; i < channels; i++)
        {
            for (int j = 0; j < reco_task->timebins; j++)
            {
                if(reco_task->padhits[i][j] != 0)
                {
                    double time = 100 * j + 50;
                    double xpad,ypad;
                    DefaultLayout::GetDefaultLayout().GetPadCenter(i+1,xpad,ypad);

                    RecoPoint reco = Reconstruct(*(m_loop->map),EndPoint(xpad,ypad,zmax,time));
                    cfit3d_reco->AddPoint(reco.x(),reco.y(),reco.z(),reco_task->padhits[i][j]);
                    reco_data.emplace_back(reco.x(),reco.y(),reco.z(),reco_task->padhits[i][j]);
                }
            }        
        }

        auto reco_markers = GetDataMarkers(reco_data);
        for (auto m : reco_markers) m->Draw("same");

        cfit3d_reco->SetFitter(4,false);
        if(m_loop->curr_loop == TrackLoop::MULTI) cfit3d_reco->SetParameters(0,10,1.5,m_loop->curr_microtrack->electron);
        cfit3d_reco->FitCircle3D();
        cfit3d_reco->PrintFitParams();

        std::cout << "Reconstructed energy (circle fit with pads): " << cfit3d_reco->GetEnergy(*magfield,true) << " (middle field), " << cfit3d_reco->GetEnergy(*magfield,false) << " (average field).\n";

        RKFit* rkfit;
        if(m_loop->curr_loop == TrackLoop::SINGLE) rkfit = new RKFit(magfield,true,1E-13,{0,0,0},{1,0,0},cfit3d_reco->GetData());
        else if(m_loop->curr_loop == TrackLoop::MULTI) rkfit = new RKFit(magfield,m_loop->curr_microtrack->electron,1E-13,m_loop->curr_microtrack->origin,m_loop->curr_microtrack->orientation,cfit3d_reco->GetData());
        rkfit->SetEnergy(cfit3d_reco->GetEnergy(*magfield,true));
        rkfit->SetFitter();
        rkfit->FitRK();
        rkfit->PrintFitParams();

        if (m_loop->make_track_plots)
        {
            TGraph2D* g_cfit3d = cfit3d->GetGraph(0.1,0);
            TGraph* g_en = cfit3d->GetEnergyGraph(*magfield);
            g_cfit3d->SetLineColor(kGreen);
            g_cfit3d->SetLineWidth(2);
            // g_cfit3d->Draw("LINE same");


            TGraph2D* g_cfit3d_reco = cfit3d_reco->GetGraph(0.1,0);
            TGraph* g_en_reco = cfit3d_reco->GetEnergyGraph(*magfield);
            g_cfit3d_reco->SetLineColor(kBlue);
            g_cfit3d_reco->SetLineWidth(2);
            g_cfit3d_reco->Draw("LINE same");


            TGraph2D* g_rkfit = rkfit->GetFitGraph();
            g_rkfit->SetLineColor(kMagenta);
            g_rkfit->SetLineWidth(2);
            g_rkfit->Draw("LINE same");

            TLegend* leg_xyz = new TLegend(0.741,0.742,0.956,0.931);
            leg_xyz->AddEntry(reco_task->g_xyz,"original trajectory","p");
            // leg_xyz->AddEntry(g_cfit3d,"original trajectory fit by circle");
            leg_xyz->AddEntry(reco_task->g_xyz_reco,"reconstructed trajectory (with pads)","p");
            leg_xyz->AddEntry(g_cfit3d_reco,"reconstructed trajectory fit by circle");
            leg_xyz->AddEntry(g_rkfit,"reconstructed trajectory fit by Runge-Kutta");
            leg_xyz->Draw("same");

            reco_task->c_reco->Write();


            std::string c_fit3d_name = "c_fit3d" + std::to_string(m_loop->curr_track_index);
            TCanvas* c_fit3d = new TCanvas(c_fit3d_name.c_str(),"Fit 3D");
                g_cfit3d_reco->SetTitle("Fit 3D;x [cm];y [cm];z [cm]");
                g_cfit3d_reco->Draw("LINE");
                g_cfit3d->Draw("LINE same");
                g_rkfit->Draw("LINE same");
                g_cfit3d_reco->GetYaxis()->SetNdivisions(508);
                g_cfit3d_reco->GetXaxis()->SetTitleOffset(1.5);
                g_cfit3d_reco->GetZaxis()->SetTitleOffset(1.5);
            
            c_fit3d->Write();

            std::string c_energy_name = "c_energy" + std::to_string(m_loop->curr_track_index);
            TCanvas* c_energy = new TCanvas(c_energy_name.c_str(),"Energy along track");
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

            c_energy->Write();
        }

        double track_E      = m_loop->curr_microtrack->kin_energy;
        double track_theta  = (180/TMath::Pi())*asin(m_loop->curr_microtrack->orientation.z);
        double track_phi    = (180/TMath::Pi())*acos(m_loop->curr_microtrack->orientation.x/cos(asin(m_loop->curr_microtrack->orientation.z)))*sign(m_loop->curr_microtrack->orientation.y);
        double E_resolution = 100*(rkfit->GetEnergy()-track_E)/track_E;

        if (m_loop->curr_loop == TrackLoop::MULTI)
        {
            h_all->Fill(track_phi,track_theta,track_E/1e+6,E_resolution);
            h_deltaenergy_energy->Fill(track_E/1e+6,E_resolution);
            // h_theta_phi->Fill(track_phi,track_theta,track_E,E_resolution);
            // h_theta_energy->Fill(track_phi,track_theta,track_E,E_resolution);
            // h_phi_energy->Fill(track_phi,track_theta,track_E,E_resolution);
        }
    }

    void PostTrackLoop() override
    {
        h_all->Write();
        h_deltaenergy_energy->Write();
    }

public:
    CircleAndRKFitTask(RecoPadsTask* t) : reco_task(t) { }
};










/// @brief Class for plotting Runge-Kutta simulated energy vs circle fit reconstruction.
class CircleFitEnergyTask : public RecoTask
{
    TH2F* h_energy;
    std::vector<TH1F*> h_energy_diff;
    int bins;
    CircleFit3D* cfit = nullptr;

    std::vector<TGraph2D*> tracks,fits;
    double Ediff_max = 2e+6;

    void PreTrackLoop() override
    {
        const char* h_energy_title = "Original RK vs circle fit energy;RK energy [MeV];circle fit energy [MeV]";
        h_energy = new TH2F("h_energy",h_energy_title,bins,4,12,bins,0,15);

        for (int i = 4; i < 12; i++)
        {
            std::string name  = "h_energy_diff_" + std::to_string(i);
            std::string title = "Difference between measured and fitted energy;Energy difference [MeV];Count";
            h_energy_diff.push_back(new TH1F(name.c_str(), title.c_str(),25,-2+i+0.5,2+i+0.5));
        }

        cfit = new CircleFit3D();
        cfit->SetFitter(4,false);
    }

    void PreElectronLoop() override
    {
        const TrackRK* track = m_loop->curr_rk;

        CircleFit3D* cfit = new CircleFit3D(track->origin,track->orientation);
        tracks.push_back(new TGraph2D());
        cfit->SetAlpha(track->electron);
    }

    void ElectronLoop() override
    {
        const TrackRK* track = m_loop->curr_rk;
        RKPoint p = m_loop->curr_rkpoint;

        cfit->AddPoint(p);tracks.back()->AddPoint(p.x(),p.y(),p.z());//g_2d->AddPoint(p.x,p.y,p.z);
    }

    void PostElectronLoop() override
    {
        const TrackRK* track = m_loop->curr_rk;
        cfit->FitCircle3D();
        // cfit->PrintFitParams();

        // g_fit = cfit->GetGraph();

        h_energy->Fill(track->kin_energy / 1e+6,cfit->GetEnergy(*(m_loop->magfield)) / 1e+6);
        h_energy_diff[floor(track->kin_energy / 1e+6) - 4]->Fill((track->kin_energy+cfit->GetEnergy(*(m_loop->magfield))-track->kin_energy) / 1e+6);
        if (abs(cfit->GetEnergy(*(m_loop->magfield))-track->kin_energy) > Ediff_max)
        {
            std::cout << "\n Electron track? --> " << track->electron << "\n";
            cfit->PrintFitParams();
            fits.push_back(cfit->GetGraph());
        }
        if (abs(cfit->GetEnergy(*(m_loop->magfield))-track->kin_energy) < Ediff_max) tracks.pop_back();         
    }

    void PostLoop()
    {
        using namespace X17::constants;
        TCanvas* c = new TCanvas("c_cfit_energy","Original vs reconstructed energy");

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

        TCanvas* c1 = new TCanvas("c_cfit_failed","Tracks with failed circle fit");
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

        for (TGraph2D* g : fits) 
        {
            g->SetLineColor(kBlue);
            g->Draw("LINE same");
        }

        std::cout << "n_failed_fits: " << tracks.size() << "\n";
    }

public:
    CircleFitEnergyTask(const int& bins = 100) : bins(bins) {}
};


/// @brief Class for plotting Runge-Kutta simulated tracks with lower than selected energy.
class PlotSelectionTask : public RecoTask
{
    std::vector<TGraph2D*> tracks;
    double E_max;
    CircleFit3D* cfit = nullptr;

    void PreTrackLoop() override
    {
        cfit = new CircleFit3D();
        cfit->SetFitter();
    }

    void PreElectronLoop() override
    {
        const TrackRK* track = m_loop->curr_rk;

        cfit = new CircleFit3D(track->origin,track->orientation);
        tracks.push_back(new TGraph2D());
    }

    void ElectronLoop() override
    {
        const TrackRK* track = m_loop->curr_rk;
        RKPoint p = m_loop->curr_rkpoint;

        cfit->AddPoint(p);
        tracks.back()->AddPoint(p.x(),p.y(),p.z());
    }

    void PostElectronLoop() override
    {   
        cfit->FitCircle3D();

        if (cfit->GetEnergy(*(m_loop->magfield)) > E_max) tracks.pop_back(); 
    }

    void PostLoop()
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

        std::cout << "n_failed_fits: " << tracks.size() << "\n";
    }

public:
    PlotSelectionTask(double E_max = 3.5e+6) : E_max(E_max) {}
};










class FitAndSaveTask : public RecoTask
{
    CircleFit3D* cfit3d = nullptr;
    RecoPadsTask* reco_task = nullptr;

    TTree* tracks_info;
    TrackInfo curr_info;

    void PreTrackLoop() override
    {
        tracks_info = new TTree("tracks_info","Simulated and reconstructed track information.");
        tracks_info->Branch("track_info",&curr_info);
    }

    void PreElectronLoop() override
    {
        if(m_loop->curr_loop == TrackLoop::SINGLE) cfit3d = new CircleFit3D({constants::xmin,0,0},{1,0,0});
        else if(m_loop->curr_loop == TrackLoop::MULTI) cfit3d = new CircleFit3D(m_loop->curr_microtrack->origin,m_loop->curr_microtrack->orientation);
    }

    void ElectronLoop() override
    {
        MicroPoint micro = m_loop->curr_micro;
        cfit3d->AddPoint(micro.x0(),micro.y0(),micro.z0(),1);
    }

    void PostElectronLoop() override
    {
        using namespace X17::constants;

        std::cout << "Simulated energy: " << m_loop->curr_microtrack->kin_energy << " eV.";

        cfit3d->SetFitter(4,false);
        cfit3d->SetParameters(0,10,1.5,m_loop->curr_microtrack->electron);
        cfit3d->FitCircle3D();
        cfit3d->PrintFitParams();

        Field<Vector>* magfield = m_loop->magfield;
        double cfit3d_E_mid = cfit3d->GetEnergy(*magfield,true);
        double cfit3d_E_avg = cfit3d->GetEnergy(*magfield,false);

        std::cout << "Reconstructed energy (circle fit no pads): " << cfit3d_E_mid << " (middle field), " << cfit3d_E_avg << " (average field).\n";

        std::vector<RecoPoint> reco_data;
        CircleFit3D* cfit3d_reco;
        if(m_loop->curr_loop == TrackLoop::SINGLE) cfit3d_reco = new CircleFit3D({0,0,0},{1,0,0});
        else if(m_loop->curr_loop == TrackLoop::MULTI) cfit3d_reco = new CircleFit3D(m_loop->curr_microtrack->origin,m_loop->curr_microtrack->orientation);

        for (int i = 0; i < channels; i++)
        {
            for (int j = 0; j < reco_task->timebins; j++)
            {
                if(reco_task->padhits[i][j] != 0)
                {
                    double time = 100 * j + 50;
                    double xpad,ypad;
                    DefaultLayout::GetDefaultLayout().GetPadCenter(i+1,xpad,ypad);

                    RecoPoint reco = Reconstruct(*(m_loop->map),EndPoint(xpad,ypad,zmax,time));
                    cfit3d_reco->AddPoint(reco.x(),reco.y(),reco.z(),reco_task->padhits[i][j]);
                    reco_data.emplace_back(reco.x(),reco.y(),reco.z(),reco_task->padhits[i][j]);
                }
            }        
        }

        auto reco_markers = GetDataMarkers(reco_data);
        for (auto m : reco_markers) m->Draw("same");

        cfit3d_reco->SetFitter(4,false);
        if(m_loop->curr_loop == TrackLoop::MULTI) cfit3d_reco->SetParameters(0,10,1.5,m_loop->curr_microtrack->electron);
        cfit3d_reco->FitCircle3D();
        cfit3d_reco->PrintFitParams();

        double cfit3d_reco_E_mid = cfit3d_reco->GetEnergy(*magfield,true);
        double cfit3d_reco_E_avg = cfit3d_reco->GetEnergy(*magfield,false);

        std::cout << "Reconstructed energy (circle fit with pads): " << cfit3d_reco_E_mid << " (middle field), " << cfit3d_reco_E_avg << " (average field).\n";

        RKFit* rkfit;
        if(m_loop->curr_loop == TrackLoop::SINGLE) rkfit = new RKFit(magfield,true,1E-13,{0,0,0},{1,0,0},cfit3d_reco->GetData());
        else if(m_loop->curr_loop == TrackLoop::MULTI) rkfit = new RKFit(magfield,m_loop->curr_microtrack->electron,1E-13,m_loop->curr_microtrack->origin,m_loop->curr_microtrack->orientation,cfit3d_reco->GetData());
        rkfit->SetEnergy(cfit3d_reco_E_mid);
        rkfit->SetFitter();
        rkfit->FitRK();
        rkfit->PrintFitParams();

        double track_E      = m_loop->curr_microtrack->kin_energy;
        double track_theta  = (180/TMath::Pi())*asin(m_loop->curr_microtrack->orientation.z);
        double track_phi    = (180/TMath::Pi())*acos(m_loop->curr_microtrack->orientation.x/cos(asin(m_loop->curr_microtrack->orientation.z)))*sign(m_loop->curr_microtrack->orientation.y);
        double E_resolution = 100*(rkfit->GetEnergy()-track_E)/track_E;

        curr_info.electron   = m_loop->curr_microtrack->electron;
        curr_info.theta      = track_theta;
        curr_info.phi        = track_phi;
        curr_info.kin_energy = track_E/1e+6;

        curr_info.cfit_nopads_energy_mid = cfit3d_E_mid/1e+6;
        curr_info.cfit_nopads_energy_avg = cfit3d_E_avg/1e+6;
        curr_info.cfit_pads_energy_mid   = cfit3d_reco_E_mid/1e+6;
        curr_info.cfit_pads_energy_avg   = cfit3d_reco_E_avg/1e+6;
        curr_info.rkfit_energy           = rkfit->GetEnergy()/1e+6;
        curr_info.rkfit_energy_err       = rkfit->GetEnergyError()/1e+6;

        tracks_info->Fill();
    }

    void PostTrackLoop() override
    {
        tracks_info->Write();
    }

public:
    FitAndSaveTask(RecoPadsTask* t) : reco_task(t) { }
};