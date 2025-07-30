#pragma once

// X17 dependencies
#include "Reconstruction.h"
#include "TrackLoop.h"
#include "Utilities.h"
#include "X17Utilities.h"

class MapRecoTask : public X17::RecoTask
{
    std::vector<double> v_res_x;
    std::vector<double> v_res_y;
    std::vector<double> v_res_z;
    std::vector<double> v_res_total;
    std::string canvas_name;
    bool electrons, positrons;
    int phi_bin, theta_bin, energy_bin;

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
        X17::RecoPoint reco = X17::Reconstruct(m_loop->map,micro);

        bool curr_electron = m_loop->curr_microtrack->electron;
        if (curr_electron && !electrons) return;
        if (!curr_electron && !positrons) return;

        v_res_x.push_back(reco.x() - micro.x0());
        v_res_y.push_back(reco.y() - micro.y0());
        v_res_z.push_back(reco.z() - micro.z0());
        v_res_total.push_back(std::sqrt(pow((reco.x() - micro.x0()),2) + pow((reco.y() - micro.y0()),2) + pow((reco.z() - micro.z0()),2)));
    }

    void PostTrackLoop() override
    {
        std::cout << "Number of electrons: " << v_res_x.size() << "\n";
        TCanvas* c = new TCanvas(canvas_name.c_str(),"Reconstruction residuals",1400,1000);
        c->Divide(2,2);
        int xbins = GetBinsFreedmanDiaconis(v_res_x,-0.4,0.4);
        int ybins = GetBinsFreedmanDiaconis(v_res_y,-0.4,0.4);
        int zbins = GetBinsFreedmanDiaconis(v_res_z,-0.4,0.4);
        int rbins = GetBinsFreedmanDiaconis(v_res_total,0,0.4);

        auto hist_style = [](TH1F* h){
            ApplyThesisStyle(h);
            h->SetStats(0);
            h->GetYaxis()->SetTitleOffset(1.35);
            h->GetXaxis()->SetTitleOffset(0.98);
        };

        double left = 0.15;
        double right = 0.05;
        double top = 0.08;
        double bottom = 0.12;

        TPad* pad = (TPad*)c->cd(1);
        TH1F* hx_res = new TH1F("","X residuals;#Deltax [cm];# of electrons",xbins,-0.4,0.4);
        for (int i = 0; i < v_res_x.size(); i++)
            hx_res->Fill(v_res_x[i]);
        hist_style(hx_res);
        pad->SetMargin(left,right,bottom,top);
        hx_res->Draw("hist");

        pad = (TPad*)c->cd(2);
        TH1F* hy_res = new TH1F("","Y residuals;#Deltay [cm];# of electrons",ybins,-0.4,0.4);
        for (int i = 0; i < v_res_y.size(); i++)
            hy_res->Fill(v_res_y[i]);
        hist_style(hy_res);
        pad->SetMargin(left,right,bottom,top);
        hy_res->Draw("hist");

        pad = (TPad*)c->cd(3);
        TH1F* hz_res = new TH1F("","Z residuals;#Deltaz [cm];# of electrons",zbins,-0.4,0.4);
        for (int i = 0; i < v_res_z.size(); i++)
            hz_res->Fill(v_res_z[i]);
        hist_style(hz_res);
        pad->SetMargin(left,right,bottom,top);
        hz_res->Draw("hist");

        pad = (TPad*)c->cd(4);
        TH1F* hr_res = new TH1F("","Residuals;Deviation [cm];# of electrons",rbins,0,0.4);
        for (int i = 0; i < v_res_total.size(); i++)
            hr_res->Fill(v_res_total[i]);
        hist_style(hr_res);
        pad->SetMargin(left,right,bottom,top);
        hr_res->Draw("hist");
        c->Write();
        delete c;

        TCanvas* c_x = new TCanvas((canvas_name + "_x").c_str(),"X residuals");
        pad = (TPad*)c_x->cd();
        TH1F* hx_res_fit = new TH1F("",";#Deltax [cm];# of electrons",xbins,-0.4,0.4);
        for (int i = 0; i < v_res_x.size(); i++)
            hx_res_fit->Fill(v_res_x[i]);
        hist_style(hx_res_fit);
        pad->SetMargin(left,right,bottom,top);
        hx_res_fit->Draw("hist");
        hx_res_fit->Fit("gaus");

        TF1* fit = hx_res_fit->GetFunction("gaus");
        fit->SetNpx(1000);
        fit->Draw("same");
        double mean = fit->GetParameter(1);
        double sigma = fit->GetParameter(2);

        TLatex* label = new TLatex();
        label->SetNDC();
        label->SetTextSize(0.055);
        label->DrawLatex(0.68, 0.85, Form("Mean: %.2f mm", 10*mean));
        label->DrawLatex(0.68, 0.78, Form("Sigma: %.2f mm", 10*sigma));

        c_x->Write();
        delete c_x;

        TCanvas* c_y = new TCanvas((canvas_name + "_y").c_str(),"Y residuals");
        pad = (TPad*)c_y->cd();
        TH1F* hy_res_fit = new TH1F("",";#Deltay [cm];# of electrons",ybins,-0.4,0.4);
        for (int i = 0; i < v_res_y.size(); i++)
            hy_res_fit->Fill(v_res_y[i]);
        hist_style(hy_res_fit);
        pad->SetMargin(left,right,bottom,top);
        hy_res_fit->Draw("hist");
        hy_res_fit->Fit("gaus");

        fit = hy_res_fit->GetFunction("gaus");
        fit->SetNpx(1000);
        fit->Draw("same");
        mean = fit->GetParameter(1);
        sigma = fit->GetParameter(2);

        label = new TLatex();
        label->SetNDC();
        label->SetTextSize(0.055);
        label->DrawLatex(0.68, 0.85, Form("Mean: %.2f mm", 10*mean));
        label->DrawLatex(0.68, 0.78, Form("Sigma: %.2f mm", 10*sigma));

        c_y->Write();
        delete c_y;

        TCanvas* c_z = new TCanvas((canvas_name + "_z").c_str(),"Z residuals");
        pad = (TPad*)c_z->cd();
        TH1F* hz_res_fit = new TH1F("",";#Deltaz [cm];# of electrons",zbins,-0.4,0.4);
        for (int i = 0; i < v_res_z.size(); i++)
            hz_res_fit->Fill(v_res_z[i]);
        hist_style(hz_res_fit);
        pad->SetMargin(left,right,bottom,top);
        hz_res_fit->Draw("hist");
        hz_res_fit->Fit("gaus");

        fit = hz_res_fit->GetFunction("gaus");
        fit->SetNpx(1000);
        fit->Draw("same");
        mean = fit->GetParameter(1);
        sigma = fit->GetParameter(2);

        label = new TLatex();
        label->SetNDC();
        label->SetTextSize(0.055);
        label->DrawLatex(0.68, 0.85, Form("Mean: %.2f mm", 10*mean));
        label->DrawLatex(0.68, 0.78, Form("Sigma: %.2f mm", 10*sigma));

        c_z->Write();
        delete c_z;

        delete hx_res;
        delete hy_res;
        delete hz_res;
        delete hr_res;
    }

public:
    MapRecoTask(std::string name, bool electrons, bool positrons) : canvas_name(name), electrons(electrons), positrons(positrons) {}
};