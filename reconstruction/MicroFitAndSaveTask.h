#pragma once

// X17 dependencies
#include "CircleFit3D.h"
#include "RK4.h"
#include "TrackLoop.h"

// Task dependencies
#include "RecoPadsTask.h"

class MicroFitAndSaveTask : public X17::RecoTask
{
  X17::CircleFit3D* m_cfit3d      = nullptr;
  X17::CircleFit3D* m_cfit3d_reco = nullptr;
  RecoPadsTask* m_reco_task       = nullptr;

  TTree* m_tracks_info;
  X17::TrackInfo m_curr_info;

public:
  void PreTrackLoop() override
  {
    m_tracks_info = new TTree("tracks_info", "Simulated and reconstructed track information.");
    m_tracks_info->Branch("track_info", &m_curr_info);

    // cfit3d = new CircleFit3D();
    m_cfit3d_reco = new X17::CircleFit3D();
  }

  void PreElectronLoop() override
  {
    // if (loop->curr_loop == X17::TrackLoop::SINGLE) cfit3d = new CircleFit3D({constants::xmin,0,0},{1,0,0});
    // else if (loop->curr_loop == X17::TrackLoop::MULTI)
    // cfit3d->SetOrigOrient(loop->curr_microtrack->origin,loop->curr_microtrack->orientation);
  }

  void ElectronLoop() override
  {
    X17::MicroPoint micro = GetLoop_()->curr_micro;
    // cfit3d->AddPoint(micro.x0(),micro.y0(),micro.z0(),1);
  }

  void PostElectronLoop() override
  {
    using namespace X17::constants;

    X17::TrackLoop* loop = GetLoop_();

    std::cout << "Simulated energy: " << loop->curr_microtrack->kin_energy << " eV.";

    // cfit3d->SetFitter(4,false);
    // cfit3d->SetParameters(0,10,1.5,loop->curr_microtrack->electron);
    // cfit3d->FitCircle3D();
    // cfit3d->PrintFitParams();

    X17::Field<X17::Vector>*& magfield = loop->magfield;
    // double cfit3d_E_mid = cfit3d->GetEnergy(*magfield,true);
    // double cfit3d_E_avg = cfit3d->GetEnergy(*magfield,false);

    // std::cout << "Reconstructed energy (circle fit no pads): " << cfit3d_E_mid << " (middle field), " << cfit3d_E_avg
    // << " (average field).\n";

    std::vector<X17::RecoPoint> reco_data;

    if (loop->curr_loop == X17::TrackLoop::SINGLE)
      m_cfit3d_reco = new X17::CircleFit3D({ 0, 0, 0 }, { 1, 0, 0 });
    else if (loop->curr_loop == X17::TrackLoop::MULTI)
      m_cfit3d_reco->SetOrigOrient(loop->curr_microtrack->origin, loop->curr_microtrack->orientation);

    for (int i = 0; i < channels; i++)
    {
      for (int j = 0; j < RecoPadsTask::s_timebins; j++)
      {
        if (m_reco_task->m_padhits[i][j] != 0)
        {
          double time = 100 * j + 50;
          double xpad, ypad;
          X17::DefaultLayout::GetDefaultLayout().GetPadCenter(i + 1, xpad, ypad);

          X17::RecoPoint reco = Reconstruct(loop->map, X17::EndPoint(xpad, ypad, zmax, time));
          m_cfit3d_reco->AddPoint(reco.x(), reco.y(), reco.z(), m_reco_task->m_padhits[i][j]);
          reco_data.emplace_back(reco.x(), reco.y(), reco.z(), m_reco_task->m_padhits[i][j]);
        }
      }
    }

    // auto reco_markers = GetDataMarkers(reco_data);
    // for (auto m : reco_markers) m->Draw("same");

    m_cfit3d_reco->SetFitter(4, false);
    if (loop->curr_loop == X17::TrackLoop::MULTI)
      m_cfit3d_reco->SetParameters(0, 10, 1.5, loop->curr_microtrack->electron);
    m_cfit3d_reco->FitCircle3D();
    // cfit3d_reco->PrintFitParams();

    double cfit3d_reco_E_mid = m_cfit3d_reco->GetEnergy(*magfield, true);
    double cfit3d_reco_E_avg = m_cfit3d_reco->GetEnergy(*magfield, false);

    std::cout << "Reconstructed energy (circle fit with pads): " << cfit3d_reco_E_mid << " (middle field), "
              << cfit3d_reco_E_avg << " (average field).\n";

    double energy_est = std::isnan(cfit3d_reco_E_avg) ? cfit3d_reco_E_mid : cfit3d_reco_E_avg;
    energy_est        = std::clamp(energy_est, 2e+6, 15e+6);

    X17::RKFit* rkfit;
    if (loop->curr_loop == X17::TrackLoop::SINGLE)
      rkfit = new X17::RKFit(magfield, true, 1E-13 * 16.66, { 0, 0, 0 }, { 1, 0, 0 }, m_cfit3d_reco->GetData());
    else if (loop->curr_loop == X17::TrackLoop::MULTI)
      rkfit = new X17::RKFit(magfield, loop->curr_microtrack->electron, 1E-13 * 16.66, loop->curr_microtrack->origin,
                             loop->curr_microtrack->orientation, m_cfit3d_reco->GetData());
    rkfit->SetEnergy(energy_est);
    rkfit->SetFitter();
    rkfit->FitRK();
    rkfit->PrintFitParams();

    double track_E      = loop->curr_microtrack->kin_energy;
    double track_theta  = loop->curr_microtrack->theta();
    double track_phi    = loop->curr_microtrack->varphi();
    double E_resolution = 100 * (rkfit->GetEnergy() - track_E) / track_E;

    m_curr_info.electron   = loop->curr_microtrack->electron;
    m_curr_info.theta      = track_theta;
    m_curr_info.phi        = track_phi;
    m_curr_info.kin_energy = track_E / 1e+6;

    m_curr_info.cfit_nopads_energy_mid = -1; // cfit3d_E_mid/1e+6;
    m_curr_info.cfit_nopads_energy_avg = -1; // cfit3d_E_avg/1e+6;
    m_curr_info.cfit_pads_energy_mid   = cfit3d_reco_E_mid / 1e+6;
    m_curr_info.cfit_pads_energy_avg   = cfit3d_reco_E_avg / 1e+6;
    m_curr_info.rkfit_energy           = rkfit->GetEnergy() / 1e+6;
    m_curr_info.rkfit_energy_err       = rkfit->GetEnergyError() / 1e+6;

    m_tracks_info->Fill();
    // cfit3d->ResetData();
    m_cfit3d_reco->ResetData();
  }

  void PostTrackLoop() override { m_tracks_info->Write(); }

public:
  MicroFitAndSaveTask(RecoPadsTask* t)
    : m_reco_task(t)
  {
  }
};