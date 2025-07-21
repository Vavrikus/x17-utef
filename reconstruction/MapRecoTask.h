#pragma once

// X17 dependencies
#include "Reconstruction.h"
#include "TrackLoop.h"
#include "Utilities.h"

class MapRecoTask : public X17::RecoTask
{
    std::vector<double> v_res_x;
    std::vector<double> v_res_y;
    std::vector<double> v_res_z;
    std::vector<double> v_res_total;
    std::string canvas_name;
    bool electrons, positrons;

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
        TCanvas* c = new TCanvas(canvas_name.c_str(),"Reconstruction residuals",1400,1000);
        c->Divide(2,2);
        int xbins = GetBinsFreedmanDiaconis(v_res_x,-0.4,0.4);
        int ybins = GetBinsFreedmanDiaconis(v_res_y,-0.4,0.4);
        int zbins = GetBinsFreedmanDiaconis(v_res_z,-0.4,0.4);
        int rbins = GetBinsFreedmanDiaconis(v_res_total,0,0.4);

        c->cd(1);
        TH1F* hx_res = new TH1F("","X residuals;#Deltax [cm];# of electrons",xbins,-0.4,0.4);
        for (int i = 0; i < v_res_x.size(); i++)
            hx_res->Fill(v_res_x[i]);
        hx_res->Draw("hist");

        c->cd(2);
        TH1F* hy_res = new TH1F("","Y residuals;#Deltay [cm];# of electrons",ybins,-0.4,0.4);
        for (int i = 0; i < v_res_y.size(); i++)
            hy_res->Fill(v_res_y[i]);
        hy_res->Draw("hist");

        c->cd(3);
        TH1F* hz_res = new TH1F("","Z residuals;#Deltaz [cm];# of electrons",zbins,-0.4,0.4);
        for (int i = 0; i < v_res_z.size(); i++)
            hz_res->Fill(v_res_z[i]);
        hz_res->Draw("hist");

        c->cd(4);
        TH1F* hr_res = new TH1F("","Residuals;Deviation [cm];# of electrons",rbins,0,0.4);
        for (int i = 0; i < v_res_total.size(); i++)
            hr_res->Fill(v_res_total[i]);
        hr_res->Draw("hist");
        c->Write();

        delete c;
        delete hx_res;
        delete hy_res;
        delete hz_res;
        delete hr_res;
    }

public:
    MapRecoTask(std::string name, bool electrons, bool positrons) : canvas_name(name), electrons(electrons), positrons(positrons) {}
};