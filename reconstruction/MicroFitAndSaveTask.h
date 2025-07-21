#pragma once

// ROOT dependencies
#include "TMath.h"

// X17 dependencies
#include "CircleFit3D.h"
#include "RK4.h"
#include "TrackLoop.h"

// Task dependencies
#include "RecoPadsTask.h"

class MicroFitAndSaveTask : public X17::RecoTask
{
    X17::CircleFit3D* cfit3d = nullptr;
    X17::CircleFit3D* cfit3d_reco = nullptr;
    RecoPadsTask* reco_task = nullptr;

    TTree* tracks_info;
    X17::TrackInfo curr_info;

    void PreTrackLoop() override
    {
        tracks_info = new TTree("tracks_info","Simulated and reconstructed track information.");
        tracks_info->Branch("track_info",&curr_info);

        // cfit3d = new CircleFit3D();
        cfit3d_reco = new X17::CircleFit3D();
    }

    void PreElectronLoop() override
    {
        // if (m_loop->curr_loop == X17::TrackLoop::SINGLE) cfit3d = new CircleFit3D({constants::xmin,0,0},{1,0,0});
        // else if (m_loop->curr_loop == X17::TrackLoop::MULTI) cfit3d->SetOrigOrient(m_loop->curr_microtrack->origin,m_loop->curr_microtrack->orientation);
    }

    void ElectronLoop() override
    {
        X17::MicroPoint micro = m_loop->curr_micro;
        // cfit3d->AddPoint(micro.x0(),micro.y0(),micro.z0(),1);
    }

    void PostElectronLoop() override
    {
        using namespace X17::constants;

        std::cout << "Simulated energy: " << m_loop->curr_microtrack->kin_energy << " eV.";

        // cfit3d->SetFitter(4,false);
        // cfit3d->SetParameters(0,10,1.5,m_loop->curr_microtrack->electron);
        // cfit3d->FitCircle3D();
        // cfit3d->PrintFitParams();

        X17::Field<X17::Vector>*& magfield = m_loop->magfield;
        // double cfit3d_E_mid = cfit3d->GetEnergy(*magfield,true);
        // double cfit3d_E_avg = cfit3d->GetEnergy(*magfield,false);

        // std::cout << "Reconstructed energy (circle fit no pads): " << cfit3d_E_mid << " (middle field), " << cfit3d_E_avg << " (average field).\n";

        std::vector<X17::RecoPoint> reco_data;
        
        if (m_loop->curr_loop == X17::TrackLoop::SINGLE) cfit3d_reco = new X17::CircleFit3D({0,0,0},{1,0,0});
        else if (m_loop->curr_loop == X17::TrackLoop::MULTI) cfit3d_reco->SetOrigOrient(m_loop->curr_microtrack->origin,m_loop->curr_microtrack->orientation);

        for (int i = 0; i < channels; i++)
        {
            for (int j = 0; j < reco_task->timebins; j++)
            {
                if (reco_task->padhits[i][j] != 0)
                {
                    double time = 100 * j + 50;
                    double xpad,ypad;
                    X17::DefaultLayout::GetDefaultLayout().GetPadCenter(i+1,xpad,ypad);

                    X17::RecoPoint reco = Reconstruct(m_loop->map,X17::EndPoint(xpad,ypad,zmax,time));
                    cfit3d_reco->AddPoint(reco.x(),reco.y(),reco.z(),reco_task->padhits[i][j]);
                    reco_data.emplace_back(reco.x(),reco.y(),reco.z(),reco_task->padhits[i][j]);
                }
            }        
        }

        // auto reco_markers = GetDataMarkers(reco_data);
        // for (auto m : reco_markers) m->Draw("same");

        cfit3d_reco->SetFitter(4,false);
        if (m_loop->curr_loop == X17::TrackLoop::MULTI) cfit3d_reco->SetParameters(0,10,1.5,m_loop->curr_microtrack->electron);
        cfit3d_reco->FitCircle3D();
        // cfit3d_reco->PrintFitParams();

        double cfit3d_reco_E_mid = cfit3d_reco->GetEnergy(*magfield,true);
        double cfit3d_reco_E_avg = cfit3d_reco->GetEnergy(*magfield,false);

        std::cout << "Reconstructed energy (circle fit with pads): " << cfit3d_reco_E_mid << " (middle field), " << cfit3d_reco_E_avg << " (average field).\n";

        double energy_est = std::isnan(cfit3d_reco_E_avg) ? cfit3d_reco_E_mid : cfit3d_reco_E_avg;
        energy_est = std::clamp(energy_est,2e+6,15e+6);

        X17::RKFit* rkfit;
        if (m_loop->curr_loop == X17::TrackLoop::SINGLE) rkfit = new X17::RKFit(magfield,true,1E-13*16.66,{0,0,0},{1,0,0},cfit3d_reco->GetData());
        else if (m_loop->curr_loop == X17::TrackLoop::MULTI) rkfit = new X17::RKFit(magfield,m_loop->curr_microtrack->electron,1E-13*16.66,m_loop->curr_microtrack->origin,m_loop->curr_microtrack->orientation,cfit3d_reco->GetData());
        rkfit->SetEnergy(energy_est);
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

        curr_info.cfit_nopads_energy_mid = -1; //cfit3d_E_mid/1e+6;
        curr_info.cfit_nopads_energy_avg = -1; //cfit3d_E_avg/1e+6;
        curr_info.cfit_pads_energy_mid   = cfit3d_reco_E_mid/1e+6;
        curr_info.cfit_pads_energy_avg   = cfit3d_reco_E_avg/1e+6;
        curr_info.rkfit_energy           = rkfit->GetEnergy()/1e+6;
        curr_info.rkfit_energy_err       = rkfit->GetEnergyError()/1e+6;

        tracks_info->Fill();
        // cfit3d->ResetData();
        cfit3d_reco->ResetData();
    }

    void PostTrackLoop() override
    {
        tracks_info->Write();
    }

public:
    MicroFitAndSaveTask(RecoPadsTask* t) : reco_task(t) { }
};