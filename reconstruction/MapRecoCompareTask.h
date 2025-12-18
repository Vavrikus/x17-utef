#pragma once

// ROOT dependencies
#include "TLatex.h"

// X17 dependencies
#include "Reconstruction.h"
#include "TrackLoop.h"
#include "Utilities.h"
#include "X17Utilities.h"

class MapRecoCompareTask : public X17::RecoTask
{
    std::vector<double> v_res_x;
    std::vector<double> v_res_y;
    std::vector<double> v_res_z;
    std::vector<double> v_res_total;
    std::string canvas_name;
    bool electrons, positrons;
    int phi_bin, theta_bin, energy_bin;
    double max_res = -1;

    void PreElectronLoop() override
    {
        using namespace X17::constants;
        // phi_bin    = floor(angle_bins  * (m_loop->curr_microtrack->varphi() - phi_min) / (phi_max - phi_min));
        // theta_bin  = floor(angle_bins  * (m_loop->curr_microtrack->theta() - theta_min) / (theta_max - theta_min));
        // energy_bin = floor(energy_bins * (m_loop->curr_microtrack->kin_energy - E_min) / (E_max - E_min));
    }

    void ElectronLoop() override
    {
        X17::MicroPoint micro = m_loop->curr_micro;
        X17::RecoPoint reco_new = X17::Reconstruct(m_loop->map,micro);
        X17::RecoPoint reco_old = X17::ReconstructOld(m_loop->map,micro.end,1E-9,false);

        bool curr_electron = m_loop->curr_microtrack->electron;
        if (curr_electron  && !electrons) return;
        if (!curr_electron && !positrons) return;

        X17::Vector oldnew_res = reco_old.AsVector() - reco_new.AsVector();

        v_res_x.push_back(oldnew_res.x);
        v_res_y.push_back(oldnew_res.y);
        v_res_z.push_back(oldnew_res.z);
        v_res_total.push_back(oldnew_res.Magnitude());

        if (oldnew_res.Magnitude() > max_res) max_res = oldnew_res.Magnitude();
    }

    void PostTrackLoop() override
    {
        std::cout << "Number of electrons: " << v_res_x.size() << "\n";
        TCanvas* c = new TCanvas(canvas_name.c_str(),"Reconstruction residuals",1400,1000);
        c->Divide(2,2);
        int xbins = GetBinsFreedmanDiaconis(v_res_x,-max_res,max_res);
        int ybins = GetBinsFreedmanDiaconis(v_res_y,-max_res,max_res);
        int zbins = GetBinsFreedmanDiaconis(v_res_z,-max_res,max_res);
        int rbins = GetBinsFreedmanDiaconis(v_res_total,0,max_res);

        auto hist_style = [](TH1F* h){
            ApplyThesisStyle(h);
            h->SetStats(0);
            h->GetYaxis()->SetTitleOffset(1.35);
            h->GetXaxis()->SetTitleOffset(0.98);
        };

        double left   = 0.15;
        double right  = 0.05;
        double top    = 0.08;
        double bottom = 0.12;

        std::vector<double> vecs[4] = {v_res_x,v_res_y,v_res_z,v_res_total};
        TH1F* hists[4];
        std::string titles[4] = {"X residuals;#Deltax [cm];# of electrons",
                                 "Y residuals;#Deltay [cm];# of electrons",
                                 "Z residuals;#Deltaz [cm];# of electrons",
                                 "Residuals;Deviation [cm];# of electrons"};
        int bins[4] = {xbins,ybins,zbins,rbins};

        for (int i = 0; i < 4; i++)
        {
            TPad* pad = (TPad*)c->cd(i+1);
            hists[i] = new TH1F("",titles[i].c_str(), bins[i],(i == 3 ? 0 : -0.4), 0.4);
            for (int j = 0; j < vecs[i].size(); j++)
                hists[i]->Fill(vecs[i][j]);
            hist_style(hists[i]);
            pad->SetMargin(left,right,bottom,top);
            hists[i]->Draw("hist");
        }

        c->Write();
        delete c;
        for (TH1F* h : hists)
            delete h;

        std::string c_literals[3] = {"_x","_y","_z"};
        std::string c_titles[3] = {"X residuals","Y residuals","Z residuals"};
        std::string h_titles[3] = {";#Deltax [cm];# of electrons",";#Deltay [cm];# of electrons",";#Deltaz [cm];# of electrons"};

        for (int i = 0; i < 3; i++)
        {
            TCanvas* c = new TCanvas((canvas_name + c_literals[i]).c_str(),c_titles[i].c_str());
            TPad* pad = (TPad*)c->cd();
            TH1F* h = new TH1F("",h_titles[i].c_str(),bins[i],-0.4,0.4);
            for (int j = 0; j < vecs[i].size(); j++)
                h->Fill(vecs[i][j]);
            hist_style(h);
            pad->SetMargin(left,right,bottom,top);
            h->Draw("hist");
            h->Fit("gaus");

            TF1* fit = h->GetFunction("gaus");
            
            // TF1* fit = GetSkewGaus(-0.4,0.4);
            // fit->SetParLimits(2,0,1);
            // fit->SetParameters(10000,0,0.05,0);
            fit->SetNpx(1000);
            // h->Fit(fit,"L");
            fit->Draw("same");

            double mean, sigma, skew, fwhm;
            // GetSkewGausStats(fit,mean,sigma,skew,fwhm);
            mean  = fit->GetParameter(1);
            sigma = fit->GetParameter(2);
            skew = 0;
            fwhm = 2*sigma*std::sqrt(2*std::log(2));

            TLatex* label = new TLatex();
            label->SetNDC();
            label->SetTextSize(0.055);
            label->DrawLatex(0.68, 0.85, Form("Mean: %.2f mm", 10*mean));
            label->DrawLatex(0.68, 0.78, Form("Sigma: %.2f mm", 10*sigma));
            // label->DrawLatex(0.68, 0.71, Form("Skew: %.2f", skew));

            c->Write();
            delete c;
        }
    }

public:
    MapRecoCompareTask(std::string name, bool electrons, bool positrons) : canvas_name(name), electrons(electrons), positrons(positrons) {}
};