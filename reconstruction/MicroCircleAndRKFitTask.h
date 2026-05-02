#pragma once

// ROOT dependencies
#include "TLegend.h"
#include "TMath.h"

// X17 dependencies
#include "CircleFit3D.h"
#include "RK4.h"
#include "TrackLoop.h"

// Task dependencies
#include "RecoPadsTask.h"

/// @brief Class for plotting circle and Runge-Kutta track fits.
class MicroCircleAndRKFitTask : public X17::RecoTask
{
  X17::CircleFit3D* m_cfit3d = nullptr; // Circle fit of the original track.
  RecoPadsTask* m_reco_task  = nullptr;

  TH3F* m_h_all;
  TH2F* m_h_theta_phi;
  TH2F* m_h_theta_energy;
  TH2F* m_h_phi_energy;
  TH2F* m_h_deltaenergy_energy;
  TH1F* m_h_delta_energy;

public:
  void PreTrackLoop() override
  {
    using namespace X17::constants;

    X17::TrackLoop* loop = GetLoop_();

    // Binning.
    int angle_bins  = 21;
    int energy_bins = 11;

    // Ranges for simulation.
    double theta_max = (180 / TMath::Pi()) * atan((win_height / 2) / xmin); // The maximal simulated theta [rad].
    double theta_min = -theta_max;                                          // The minimal simulated theta [rad].
    double phi_max   = (180 / TMath::Pi()) * atan((win_width / 2) / xmin);  // The maximal simulated phi [rad].
    double phi_min   = -phi_max;                                            // The minimal simulated phi [rad].
    double E_max     = 13;                                                  // The maximal simulated energy [MeV].
    double E_min     = 3;                                                   // THe minimal simulated energy [MeV].

    // Adjusting boundaries for the extra bin.
    phi_max   = phi_min + (phi_max - phi_min) * angle_bins / (angle_bins - 1);
    theta_max = theta_min + (theta_max - theta_min) * angle_bins / (angle_bins - 1);
    E_max     = E_min + (E_max - E_min) * energy_bins / (energy_bins - 1);

    // Histogram titles.
    const char* h_all_title
      = "Energy resolution (#DeltaE/E);Phi [deg];Theta [deg];Simulated energy [MeV];#DeltaE/E [\%]";
    const char* h_theta_phi_title = "Energy resolution (#DeltaE/E);Phi [deg];Theta [deg];#DeltaE/E [\%]";
    const char* h_theta_energy_title
      = "Energy resolution (#DeltaE/E);Theta [deg];Simulated energy [MeV];#DeltaE/E [\%]";
    const char* h_phi_energy_title = "Energy resolution (#DeltaE/E);Phi [deg];Simulated energy [MeV];#DeltaE/E [\%]";

    // Histogram initializations.
    m_h_all = new TH3F("h_all", h_all_title, angle_bins, phi_min, phi_max, angle_bins, theta_min, theta_max,
                       energy_bins, E_min, E_max);
    m_h_theta_phi
      = new TH2F("h_theta_phi", h_theta_phi_title, angle_bins, phi_min, phi_max, angle_bins, theta_min, theta_max);
    m_h_theta_energy
      = new TH2F("h_theta_energy", h_theta_energy_title, angle_bins, theta_min, theta_max, energy_bins, E_min, E_max);
    m_h_phi_energy
      = new TH2F("h_phi_energy", h_phi_energy_title, angle_bins, phi_min, phi_max, energy_bins, E_min, E_max);
    m_h_deltaenergy_energy = new TH2F("h_deltaenergy_energy", "", energy_bins, E_min, E_max, 100, -50, 50);
  }

  void PreElectronLoop() override
  {
    X17::TrackLoop* loop = GetLoop_();

    if (loop->curr_loop == X17::TrackLoop::SINGLE)
      m_cfit3d = new X17::CircleFit3D({ X17::constants::xmin, 0, 0 }, { 1, 0, 0 });
    else if (loop->curr_loop == X17::TrackLoop::MULTI)
      m_cfit3d = new X17::CircleFit3D(loop->curr_microtrack->origin, loop->curr_microtrack->orientation);
  }

  void ElectronLoop() override
  {
    X17::TrackLoop* loop = GetLoop_();

    X17::MicroPoint micro = loop->curr_micro;
    m_cfit3d->AddPoint(micro.x0(), micro.y0(), micro.z0(), 1);
  }

  void PostElectronLoop() override
  {
    using namespace X17::constants;

    X17::TrackLoop* loop = GetLoop_();

    if (loop->curr_loop == X17::TrackLoop::MULTI)
      std::cout << "Simulated energy: " << loop->curr_microtrack->kin_energy << " eV.";

    // Circle fit of the original track.
    m_cfit3d->SetFitter(4, false);
    if (loop->curr_loop == X17::TrackLoop::MULTI)
      m_cfit3d->SetParameters(0, 20, 0.6, loop->curr_microtrack->electron);
    m_cfit3d->FitCircle3D();
    m_cfit3d->PrintFitParams();

    X17::Field<X17::Vector>* magfield = loop->magfield;
    std::cout << "Reconstructed energy (circle fit no pads): " << m_cfit3d->GetEnergy(*magfield, true)
              << " (middle field), ";
    std::cout << m_cfit3d->GetEnergy(*magfield, false) << " (average field).\n";

    // Circle fit of the reconstructed track.
    X17::CircleFit3D* cfit3d_reco;
    if (loop->curr_loop == X17::TrackLoop::SINGLE)
      cfit3d_reco = new X17::CircleFit3D({ 0, 0, 0 }, { 1, 0, 0 });
    else if (loop->curr_loop == X17::TrackLoop::MULTI)
      cfit3d_reco = new X17::CircleFit3D(loop->curr_microtrack->origin, loop->curr_microtrack->orientation);

    for (X17::RecoPoint p : m_reco_task->m_reco_data)
      cfit3d_reco->AddPoint(p);

    cfit3d_reco->SetFitter(4, false);
    if (loop->curr_loop == X17::TrackLoop::MULTI)
      cfit3d_reco->SetParameters(0, 20, 0.6, loop->curr_microtrack->electron);
    cfit3d_reco->FitCircle3D();
    cfit3d_reco->PrintFitParams();

    std::cout << "Reconstructed energy (circle fit with pads): " << cfit3d_reco->GetEnergy(*magfield, true)
              << " (middle field), ";
    std::cout << cfit3d_reco->GetEnergy(*magfield, false) << " (average field).\n";

    // RK fit of the reconstructed track.
    X17::RKFit* rkfit;
    if (loop->curr_loop == X17::TrackLoop::SINGLE)
      rkfit = new X17::RKFit(magfield, true, 1E-13 * 16.66, { 0, 0, 0 }, { 1, 0, 0 }, cfit3d_reco->GetData());
    else if (loop->curr_loop == X17::TrackLoop::MULTI)
      rkfit = new X17::RKFit(magfield, loop->curr_microtrack->electron, 1E-13 * 16.66, loop->curr_microtrack->origin,
                             loop->curr_microtrack->orientation, cfit3d_reco->GetData());
    rkfit->SetEnergy(cfit3d_reco->GetEnergy(*magfield, true));
    rkfit->SetFitter();
    rkfit->FitRK();
    rkfit->PrintFitParams();

    if (loop->make_track_plots)
    {
      TGraph2D* g_cfit3d    = m_cfit3d->GetGraph(0.1, 0);
      TGraph* g_en          = m_cfit3d->GetEnergyGraph(*magfield);
      TPolyLine3D* l_cfit3d = GetLine3D(g_cfit3d);
      l_cfit3d->SetLineColor(kGreen);
      l_cfit3d->SetLineWidth(2);
      // l_cfit3d->Draw("same");

      TGraph2D* g_cfit3d_reco    = cfit3d_reco->GetGraph(0.1, 0);
      TGraph* g_en_reco          = cfit3d_reco->GetEnergyGraph(*magfield);
      TPolyLine3D* l_cfit3d_reco = GetLine3D(g_cfit3d_reco);
      l_cfit3d_reco->SetLineColor(kBlue);
      l_cfit3d_reco->SetLineWidth(2);
      l_cfit3d_reco->Draw("same");

      TGraph2D* g_rkfit    = rkfit->GetFitGraph();
      TPolyLine3D* l_rkfit = GetLine3D(g_rkfit);
      l_rkfit->SetLineColor(kMagenta);
      l_rkfit->SetLineWidth(2);
      l_rkfit->Draw("same");

      TLegend* leg_xyz = new TLegend(0.741, 0.742, 0.956, 0.931);
      leg_xyz->AddEntry(m_reco_task->m_g_xyz, "original trajectory", "p");
      // leg_xyz->AddEntry(g_cfit3d,"original trajectory fit by circle");
      leg_xyz->AddEntry(m_reco_task->m_g_xyz_reco, "reconstructed trajectory (with pads)", "p");
      leg_xyz->AddEntry(g_cfit3d_reco, "reconstructed trajectory fit by circle");
      leg_xyz->AddEntry(g_rkfit, "reconstructed trajectory fit by Runge-Kutta");
      // leg_xyz->Draw("same");

      m_reco_task->m_c_reco->Write(nullptr, TObject::kOverwrite);

      std::string c_fit3d_name = "c_fit3d" + std::to_string(loop->curr_track_index);
      TCanvas* c_fit3d         = new TCanvas(c_fit3d_name.c_str(), "Fit 3D");
      g_cfit3d_reco->SetTitle("Fit 3D;x [cm];y [cm];z [cm]");
      g_cfit3d_reco->Draw("LINE");
      g_cfit3d->Draw("LINE same");
      g_rkfit->Draw("LINE same");
      g_cfit3d_reco->GetYaxis()->SetNdivisions(508);
      g_cfit3d_reco->GetXaxis()->SetTitleOffset(1.5);
      g_cfit3d_reco->GetZaxis()->SetTitleOffset(1.5);

      c_fit3d->Write();

      std::string c_energy_name = "c_energy" + std::to_string(loop->curr_track_index);
      TCanvas* c_energy         = new TCanvas(c_energy_name.c_str(), "Energy along track");
      g_en->SetMarkerColor(kRed);
      g_en->SetMarkerSize(1);
      g_en->SetMarkerStyle(2);
      g_en->Draw("AP");
      g_en_reco->SetMarkerColor(kBlue);
      g_en_reco->SetMarkerSize(1);
      g_en_reco->SetMarkerStyle(2);
      g_en_reco->Draw("P same");

      TLegend* leg_en = new TLegend(0.129, 0.786, 0.360, 0.887);
      leg_en->AddEntry(g_en, "fit original");
      leg_en->AddEntry(g_en_reco, "fit reconstructed");
      leg_en->Draw("same");

      c_energy->Write();
    }

    double track_E      = loop->curr_microtrack->kin_energy;
    double track_theta  = (180 / TMath::Pi()) * asin(loop->curr_microtrack->orientation.z);
    double track_phi    = (180 / TMath::Pi())
                          * acos(loop->curr_microtrack->orientation.x / cos(asin(loop->curr_microtrack->orientation.z)))
                          * sign(loop->curr_microtrack->orientation.y);
    double E_resolution = 100 * (rkfit->GetEnergy() - track_E) / track_E;

    if (loop->curr_loop == X17::TrackLoop::MULTI)
    {
      m_h_all->Fill(track_phi, track_theta, track_E / 1e+6, E_resolution);
      m_h_deltaenergy_energy->Fill(track_E / 1e+6, E_resolution);
      // h_theta_phi->Fill(track_phi,track_theta,track_E,E_resolution);
      // h_theta_energy->Fill(track_phi,track_theta,track_E,E_resolution);
      // h_phi_energy->Fill(track_phi,track_theta,track_E,E_resolution);
    }
  }

  void PostTrackLoop() override
  {
    m_h_all->Write();
    m_h_deltaenergy_energy->Write();
  }

public:
  MicroCircleAndRKFitTask(RecoPadsTask* t)
    : m_reco_task(t)
  {
  }
};