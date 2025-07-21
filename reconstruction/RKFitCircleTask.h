#pragma once

// ROOT dependencies
#include "TF1.h"
#include "TLatex.h"
#include "TStyle.h"

// X17 dependencies
#include "CircleFit3D.h"
#include "Track.h"
#include "TrackLoop.h"
#include "X17Utilities.h"


/// @brief Class for plotting Runge-Kutta simulated energy vs circle fit reconstruction.
class RKFitCircleTask : public X17::RecoTask
{
    friend class PlotSelectionTask;

    TH2F* h_energy, *h_energy2;
    TH2F* h_diff_energy, *h_diff_energy2;
    TH2F* h_diff_theta, *h_diff_theta2;
    TH2F* h_diff_phi, *h_diff_phi2;

    std::vector<double> v_diff_mid, v_diff_avg;

    std::vector<TH1F*> h_energy_diff;
    int bins;
    X17::CircleFit3D* cfit = nullptr;

    std::vector<TGraph2D*> tracks,fits;
    double Ediff_max = 2e+6;

    void PreTrackLoop() override
    {
        v_diff_mid.reserve(100000);
        v_diff_avg.reserve(100000); 

        const char* h_energy_title = ";RK energy [MeV];Circle fit energy [MeV]";
        h_energy = new TH2F("h_energy",h_energy_title,bins,3,13,bins,0,15);
        h_energy2 = new TH2F("h_energy2",h_energy_title,bins,3,13,bins,0,15);

        const char* h_energy_title2 = ";RK energy [MeV];#frac{E_{rec}-E}{E} [\%]";
        h_diff_energy = new TH2F("h_diff_energy",h_energy_title2,bins,3,13,bins,-25,25);
        h_diff_energy2 = new TH2F("h_diff_energy2",h_energy_title2,bins,3,13,bins,-30,10);

        const char* h_energy_title3 = ";Theta [deg];#frac{E_{rec}-E}{E} [\%]";
        h_diff_theta = new TH2F("h_diff_theta",h_energy_title3,bins,-17.1,17.1,bins,-25,25);
        h_diff_theta2 = new TH2F("h_diff_theta2",h_energy_title3,bins,-17.1,17.1,bins,-40,10);

        const char* h_energy_title4 = ";Phi [deg];#frac{E_{rec}-E}{E} [\%]";
        h_diff_phi = new TH2F("h_diff_phi",h_energy_title4,bins,-16.3,16.3,bins,-25,25);
        h_diff_phi2 = new TH2F("h_diff_phi2",h_energy_title4,bins,-16.3,16.3,bins,-30,10);

        // Simulation from 3 to 13 MeV
        // for (int i = 3; i < 13; i++)
        // {
        //     std::string name  = "h_energy_diff_" + std::to_string(i);
        //     std::string title = "Difference between measured and fitted energy;Energy difference [MeV];Count";
        //     h_energy_diff.push_back(new TH1F(name.c_str(), title.c_str(),25,-2+i+0.5,2+i+0.5));
        // }

        cfit = new X17::CircleFit3D();
        cfit->SetFitter(4,false);
    }

    void PreElectronLoop() override
    {
        const X17::TrackRK* track = m_loop->curr_rk;

        cfit->SetOrigOrient(track->origin,track->orientation);
        tracks.push_back(new TGraph2D());
        cfit->SetFitter(4,false);
        cfit->SetParameters(0,20,0.5,track->electron);
    }

    void ElectronLoop() override
    {
        const X17::TrackRK* track = m_loop->curr_rk;
        X17::RKPoint p = m_loop->curr_rkpoint;

        cfit->AddPoint(p);
        tracks.back()->AddPoint(p.x(),p.y(),p.z());//g_2d->AddPoint(p.x,p.y,p.z);
    }

    void PostElectronLoop() override
    {
        const X17::TrackRK* track = m_loop->curr_rk;
        std::cout << "Prefit sumsq: " << cfit->GetSumSq() << std::endl;
        cfit->FitCircle3D();
        // cfit->PrintFitParams();

        // g_fit = cfit->GetGraph();

        double fit_energy_mid = cfit->GetEnergy(*(m_loop->magfield),true);
        double fit_energy_avg = cfit->GetEnergy(*(m_loop->magfield),false);

        double rel_diff_mid = 100*(fit_energy_mid - track->kin_energy) / track->kin_energy;
        double rel_diff_avg = 100*(fit_energy_avg - track->kin_energy) / track->kin_energy;
        
        double theta = 180/M_PI*std::asin(track->orientation.z);
        double phi = 180/M_PI*std::acos(track->orientation.x/std::cos(std::asin(track->orientation.z)))*sign(track->orientation.y);
        if(track->electron){
        h_energy->Fill(track->kin_energy / 1e+6, fit_energy_mid / 1e+6);
        h_energy2->Fill(track->kin_energy / 1e+6, fit_energy_avg / 1e+6);
        v_diff_mid.push_back(rel_diff_mid);
        v_diff_avg.push_back(rel_diff_avg);

        h_diff_energy->Fill(track->kin_energy / 1e+6, rel_diff_mid);
        h_diff_energy2->Fill(track->kin_energy / 1e+6, rel_diff_avg);

        h_diff_theta->Fill(theta, rel_diff_mid);
        h_diff_theta2->Fill(theta, rel_diff_avg);
        h_diff_phi->Fill(phi, rel_diff_mid);
        h_diff_phi2->Fill(phi, rel_diff_avg);}

        // int index = floor(track->kin_energy / 1e+6) - 3; // Zeroth bin is 3+ MeV
        // h_energy_diff[index]->Fill((fit_energy_mid-track->kin_energy) / 1e+6);

        if (abs(fit_energy_mid - track->kin_energy) > Ediff_max)
        {
            std::cout << "FAILED TRACK:\n";
            std::cout << "   electron: " << track->electron << " Ek: " << track->kin_energy;
            std::cout << " origin: (" << track->origin.x << "," << track->origin.y << "," << track->origin.z << ")\n";
            std::cout << "   orientation: (" << track->orientation.x << "," << track->orientation.y << "," << track->orientation.z << ")\n";
            std::cout << "   theta: " << theta << " phi: " << phi << "\n";
            std::cout << "   fit energy: " << fit_energy_mid << "\n";

            cfit->PrintFitParams();
            fits.push_back(cfit->GetGraph());
        }
        else tracks.pop_back();
        cfit->ResetData();
    }

    void PostTrackLoop() override
    {
        using namespace X17::constants;
        double d_width = 700;
        double d_height = 500;
        double extra_margin_w = 0.03;
        double extra_margin_h = 0.075;
        double margin_w = 0.1 + extra_margin_w;
        double margin_h = 0.1 + extra_margin_h;
        double margin_top = 0.05;

        double extra_margin_rl = 0.035;

        double width = d_width/((0.9-extra_margin_w-extra_margin_rl)/0.9);
        double height = d_height/((0.9-extra_margin_h)/0.9);

        double x_offset = 1.18;
        double y_offset = 1;

        gStyle->SetNumberContours(255);
        for (auto* h : {h_energy,h_energy2,h_diff_energy,h_diff_energy2,h_diff_theta,h_diff_theta2,h_diff_phi,h_diff_phi2})
        {
            h->SetStats(0);
            ApplyThesisStyle(h);
            h->GetXaxis()->SetTitleOffset(y_offset);
            h->GetYaxis()->SetTitleOffset(x_offset-0.05);
            h->SetContour(256);
        }

        TCanvas* c = new TCanvas("c_cfit_energy_mid","Original vs reconstructed energy",width,height);
        c->SetLeftMargin(margin_h);
        c->SetBottomMargin(margin_w);
        c->SetTopMargin(margin_top);
        c->SetRightMargin(0.1 + extra_margin_rl);
        
        h_energy->Draw("colz");
        TF1 *f1 = new TF1("f1","x",3,13);
        f1->Draw("same");
        c->Write();

        c = new TCanvas("c_cfit_energy_avg","Original vs reconstructed energy",width,height);
        c->SetLeftMargin(margin_h);
        c->SetBottomMargin(margin_w);
        c->SetTopMargin(margin_top);
        h_energy2->Draw("colz");
        TF1 *f2 = new TF1("f2","x",3,13);
        f2->Draw("same");
        c->Write();

        c = new TCanvas("c_cfit_diff_energy_mid","Original vs reconstructed energy",width,height);
        c->SetLeftMargin(margin_h);
        c->SetBottomMargin(margin_w);
        c->SetTopMargin(margin_top);
        c->SetRightMargin(0.1 + extra_margin_rl);

        h_diff_energy->Draw("colz");
        TF1 *f3 = new TF1("f3","0",3,13);
        f3->Draw("same");
        c->Write();

        c = new TCanvas("c_cfit_diff_energy_avg","Original vs reconstructed energy",width,height);
        c->SetLeftMargin(margin_h);
        c->SetBottomMargin(margin_w);
        c->SetTopMargin(margin_top);
        c->SetRightMargin(0.1 + extra_margin_rl);

        h_diff_energy2->Draw("colz");
        TF1 *f4 = new TF1("f4","0",3,13);
        f4->Draw("same");
        c->Write();

        c = new TCanvas("c_cfit_diff_theta_mid","Original vs reconstructed energy",width,height);
        c->SetLeftMargin(margin_h);
        c->SetBottomMargin(margin_w);
        c->SetTopMargin(margin_top);
        c->SetRightMargin(0.1 + extra_margin_rl);

        h_diff_theta->Draw("colz");
        TF1 *f5 = new TF1("f5","0",3,13);
        f5->Draw("same");
        c->Write();

        c = new TCanvas("c_cfit_diff_theta_avg","Original vs reconstructed energy",width,height);
        c->SetLeftMargin(margin_h);
        c->SetBottomMargin(margin_w);
        c->SetTopMargin(margin_top);
        c->SetRightMargin(0.1 + extra_margin_rl);

        h_diff_theta2->Draw("colz");
        TF1 *f6 = new TF1("f6","0",3,13);
        f6->Draw("same");
        c->Write();

        c = new TCanvas("c_cfit_diff_phi_mid","Original vs reconstructed energy",width,height);
        c->SetLeftMargin(margin_h);
        c->SetBottomMargin(margin_w);
        c->SetTopMargin(margin_top);
        c->SetRightMargin(0.1 + extra_margin_rl);

        h_diff_phi->Draw("colz");
        TF1 *f7 = new TF1("f7","0",3,13);
        f7->Draw("same");
        c->Write();

        c = new TCanvas("c_cfit_diff_phi_avg","Original vs reconstructed energy",width,height);
        c->SetLeftMargin(margin_h);
        c->SetBottomMargin(margin_w);
        c->SetTopMargin(margin_top);
        c->SetRightMargin(0.1 + extra_margin_rl);

        h_diff_phi2->Draw("colz");
        TF1 *f8 = new TF1("f8","0",3,13);
        f8->Draw("same");
        c->Write();

        std::cout << "Final histograms...\n";
        
        double bins_mid = GetBinsFreedmanDiaconis(v_diff_mid,-25,25);
        double bins_avg = GetBinsFreedmanDiaconis(v_diff_avg,-20,0);
        
        TH1F* h_diff_mid = new TH1F("h_diff_mid",";#frac{E_{rec}-E}{E} [\%];Number of tracks",bins_mid,-25,25);
        TH1F* h_diff_avg = new TH1F("h_diff_avg",";#frac{E_{rec}-E}{E} [\%];Number of tracks",bins_avg,-20,0);
        h_diff_mid->SetStats(0);
        h_diff_avg->SetStats(0);
        ApplyThesisStyle(h_diff_mid);
        ApplyThesisStyle(h_diff_avg);
        h_diff_mid->GetXaxis()->SetTitleOffset(x_offset);
        h_diff_mid->GetYaxis()->SetTitleOffset(y_offset+0.2);
        h_diff_avg->GetXaxis()->SetTitleOffset(x_offset);
        h_diff_avg->GetYaxis()->SetTitleOffset(y_offset+0.2);

        double width2 = d_width/((0.9-extra_margin_w-extra_margin_rl)/0.9);
        double height2 = d_height/((0.9-extra_margin_h)/0.9);
        
        for (double diff : v_diff_mid) h_diff_mid->Fill(diff);
        for (double diff : v_diff_avg) h_diff_avg->Fill(diff);
        
        c = new TCanvas("c_cfit_diff_mid","Original vs reconstructed energy",width2,height2);
        c->SetLeftMargin(margin_w+extra_margin_rl);
        c->SetBottomMargin(margin_h);
        c->SetTopMargin(margin_top);
        c->SetRightMargin(margin_top);
        h_diff_mid->Draw("hist");
        double fwhm = GetFWHM(h_diff_mid,true);
        double mean = h_diff_mid->GetMean();
        double stdev = h_diff_mid->GetStdDev();
        
        // Add text labels with information
        TLatex* label = new TLatex();
        label->SetNDC();
        label->SetTextSize(0.055);
        
        label->DrawLatex(0.67, 0.86, Form("Mean: %.2f %%", mean));
        label->DrawLatex(0.67, 0.785, Form("FWHM: %.2f %%", fwhm));
        label->DrawLatex(0.67, 0.71, Form("RMS: %.2f %%", stdev));
        c->Write();
        
        c = new TCanvas("c_cfit_diff_avg","Original vs reconstructed energy",width2,height2);
        c->SetLeftMargin(margin_w+extra_margin_rl);
        c->SetBottomMargin(margin_h);
        c->SetTopMargin(margin_top);
        c->SetRightMargin(margin_top);
        h_diff_avg->Draw("hist");
        fwhm = GetFWHM(h_diff_avg,true);
        mean = h_diff_avg->GetMean();
        stdev = h_diff_avg->GetStdDev();

        // Add text labels with information
        label = new TLatex();
        label->SetNDC();
        label->SetTextSize(0.055);

        label->DrawLatex(0.2, 0.86, Form("Mean: %.2f %%", mean));
        label->DrawLatex(0.2, 0.785, Form("FWHM: %.2f %%", fwhm));
        label->DrawLatex(0.2, 0.71, Form("RMS: %.2f %%", stdev));
        c->Write();

        // TCanvas* c = new TCanvas();
        // g_2d->Draw("LINE");
        // g_2d->SetLineColor(kBlue);
        // g_2d->SetTitle("RK track (blue) and its fit (red);x [cm];y [cm];z [cm]");
        // g_fit->Draw("LINE same");
        // g_fit->SetLineColor(kRed);

        // TCanvas* c0 = new TCanvas("c_energy_diff","Energy difference");
        // TH1F* h1 = c0->DrawFrame(-2,0,14,50);
        // for(TH1F* h : h_energy_diff) h->Draw("same plc hist");
        // gStyle->SetOptStat(1111);
        // c0->Write();

        // TCanvas* c1 = new TCanvas("c_cfit_failed","Tracks with failed circle fit");
        // // A histogram for scalling of the axes.
        // TH3F* scale = new TH3F("scale","Tracks with failed circle fit;x [cm];y [cm];z [cm]",1,xmin,xmax,1,-yhigh,yhigh,1,zmin,zmax);
        // scale->Draw("");
        // gStyle->SetOptStat(0);
        // scale->GetXaxis()->SetTitleOffset(1.5);
        // scale->GetZaxis()->SetTitleOffset(1.5);

        // for (TGraph2D* g : tracks) 
        // {
        //     g->SetLineColor(kRed);
        //     g->Draw("LINE same");
        // }

        // for (TGraph2D* g : fits) 
        // {
        //     g->SetLineColor(kBlue);
        //     g->Draw("LINE same");
        // }

        std::cout << "n_failed_fits: " << tracks.size() << "\n";
        // c1->Write();
    }

public:
    RKFitCircleTask(const int& bins = 100) : bins(bins) {}
};
