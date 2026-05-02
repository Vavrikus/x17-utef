#pragma once

// C++ dependencies
#include <string>

// ROOT dependencies
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TTree.h"

// X17 dependencies
#include "CircleFit3D.h"
#include "TrackLoop.h"
#include "X17Utilities.h"

// Tasks, that were moved to separate files
#include "MapRecoCompareTask.h"      // IWYU pragma: keep
#include "MapRecoTask.h"             // IWYU pragma: keep
#include "MicroCircleAndRKFitTask.h" // IWYU pragma: keep
#include "MicroFitAndSaveTask.h"     // IWYU pragma: keep
#include "RecoPadsTask.h"            // IWYU pragma: keep
#include "RKFitCircleTask.h"         // IWYU pragma: keep

using namespace X17;

/// @brief Class for plotting drift time vs distance to readout with linear fit.
class DriftTimeTask : public RecoTask
{
public:
  void PostElectronLoop() override
  {
    X17::TrackLoop* loop = GetLoop_();

    std::string c_name = "c_drift" + std::to_string(loop->curr_track_index);
    TCanvas* c         = new TCanvas(c_name.c_str(), "Drift time");

    TTree* electrons = loop->curr_micro_tree;
    electrons->Draw("t1:8-z0", "z1>7.0");
    TGraph* tz = new TGraph(static_cast<int>(electrons->GetSelectedRows()), electrons->GetV2(), electrons->GetV1());
    tz->SetTitle("Drift time as function of distance;distance to readout [cm]; time [ns]");
    tz->SetMarkerStyle(2);
    tz->SetMarkerSize(0.4);
    tz->Draw("ap");
    tz->Fit("pol1", "", "", 8.05, 11);

    c->Write();
  }
};

/// @brief Class for plotting XZ of original and reconstructed points of the track.
class XZPlotTask : public RecoTask
{
  TGraph *m_xz, *m_xz_reco;

public:
  void PreElectronLoop() override
  {
    m_xz      = new TGraph();
    m_xz_reco = new TGraph();
  }

  void ElectronLoop() override
  {
    X17::TrackLoop* loop = GetLoop_();

    MicroPoint micro = loop->curr_micro;
    RecoPoint reco   = loop->curr_reco;

    m_xz->AddPoint(micro.x0(), 8 - micro.z0());
    m_xz_reco->AddPoint(reco.x(), 8 - reco.z());

    if (std::abs(reco.z() - micro.z0()) > 0.5)
    {
      std::cout << "High z deviation: " << reco.z() - micro.z0() << "\n";
      std::cout << "Simulated coordinates:     x0 = " << micro.x0() << ", y0 = " << micro.y0()
                << ", z0 = " << micro.z0() << ", t0 = " << micro.t0() << ", e0 = " << micro.e0 << "\n";
      std::cout << "                           x1 = " << micro.x1() << ", y1 = " << micro.y1()
                << ", z1 = " << micro.z1() << ", t1 = " << micro.t1() << ", e1 = " << micro.e1 << "\n";
      std::cout << "Reconstructed coordinates: x  = " << reco.x() << ", y  = " << reco.y() << ", z  = " << reco.z()
                << ", count = " << reco.count << "\n\n";
    }
  }

  void PostElectronLoop() override
  {
    X17::TrackLoop* loop = GetLoop_();

    std::string c_name = "c_track_xz" + std::to_string(loop->curr_track_index);
    TCanvas* c         = new TCanvas(c_name.c_str(), "Electron track reconstruction");

    m_xz_reco->SetTitle("Electron track reconstruction;x [cm]; distance to readout [cm]");
    m_xz_reco->SetMarkerStyle(2);
    m_xz_reco->SetMarkerSize(0.4);
    m_xz_reco->Draw("ap");

    m_xz->SetMarkerColor(2);
    m_xz->SetMarkerStyle(7);
    m_xz->SetMarkerSize(1.2);
    m_xz->Draw("p same");

    TLegend* leg_xz = new TLegend(0.129, 0.786, 0.360, 0.887);
    leg_xz->AddEntry(m_xz, "ionization vertices", "p");
    leg_xz->AddEntry(m_xz_reco, "reconstructed", "p");
    leg_xz->Draw("same");

    // Fitting both tracks with circles.
    // double min = 0;
    // double max = 15;
    // TF1* circle_fit  = FitCircle2(xz,min,max);
    // TF1* circle_fit2 = FitCircle2(xz_reco,min,max);

    // TGraph* magnetic_x  = new TGraph();
    // TGraph* magnetic_x2 = new TGraph();

    // Field<Vector>* magfield = loop->magfield;
    // double step = 0.1;
    // RecoEnergy(circle_fit,*magfield,magnetic_x,min,max,step);
    // RecoEnergy(circle_fit2,*magfield,magnetic_x2,min,max,step);

    c->Write();
  }
};

/// @brief Class for plotting XY of original and reconstructed points of the track.
class XYPlotTask : public RecoTask
{
  TGraph *m_xy, *m_xy_reco;

public:
  void PreElectronLoop() override
  {
    m_xy      = new TGraph();
    m_xy_reco = new TGraph();
  }

  void ElectronLoop() override
  {
    X17::TrackLoop* loop = GetLoop_();

    MicroPoint micro = loop->curr_micro;
    RecoPoint reco   = loop->curr_reco;

    m_xy->AddPoint(micro.x0(), micro.y0());
    m_xy_reco->AddPoint(reco.x(), reco.y());
  }

  void PostElectronLoop() override
  {
    X17::TrackLoop* loop = GetLoop_();

    std::string c_name = "c_track_xy" + std::to_string(loop->curr_track_index);
    TCanvas* c         = new TCanvas(c_name.c_str(), "Electron track reconstruction");

    m_xy_reco->SetTitle("Electron track reconstruction;x [cm]; y [cm]");
    m_xy_reco->SetMarkerStyle(2);
    m_xy_reco->SetMarkerSize(0.4);
    m_xy_reco->Draw("ap");

    m_xy->SetMarkerColor(2);
    m_xy->SetMarkerStyle(7);
    m_xy->SetMarkerSize(1.2);
    m_xy->Draw("p same");

    TLegend* leg_xy = new TLegend(0.129, 0.786, 0.360, 0.887);
    leg_xy->AddEntry(m_xy, "ionization vertices", "p");
    leg_xy->AddEntry(m_xy_reco, "reconstructed", "p");
    leg_xy->Draw("same");

    c->Write();
  }
};

/// @brief Class for plotting residues of reconstructed interaction points.
class GraphResTask : public RecoTask
{
  TGraph *m_gx_res, *m_gy_res, *m_gz_res, *m_gr_res;

public:
  void PreElectronLoop() override
  {
    m_gx_res = new TGraph();
    m_gy_res = new TGraph();
    m_gz_res = new TGraph();
    m_gr_res = new TGraph();
  }

  void ElectronLoop() override
  {
    X17::TrackLoop* loop = GetLoop_();

    MicroPoint micro = loop->curr_micro;
    RecoPoint reco   = loop->curr_reco;

    m_gx_res->AddPoint(micro.x0(), reco.x() - micro.x0());
    m_gy_res->AddPoint(micro.x0(), reco.y() - micro.y0());
    m_gz_res->AddPoint(micro.x0(), reco.z() - micro.z0());
    m_gr_res->AddPoint(micro.x0(), std::sqrt(pow((reco.x() - micro.x0()), 2) + pow((reco.y() - micro.y0()), 2)
                                             + pow((reco.z() - micro.z0()), 2)));
  }

  void PostElectronLoop() override
  {
    X17::TrackLoop* loop = GetLoop_();

    std::string c_name = "c_fit_res" + std::to_string(loop->curr_track_index);
    TCanvas* c         = new TCanvas(c_name.c_str(), "Reconstruction residuals");

    c->Divide(2, 2);

    c->cd(1);
    m_gx_res->SetTitle("X residuals;x [cm];#Deltax [cm]");
    m_gx_res->SetMarkerStyle(2);
    m_gx_res->SetMarkerSize(0.4);
    m_gx_res->Draw("ap");

    c->cd(2);
    m_gy_res->SetTitle("Y residuals;x [cm];#Deltay [cm]");
    m_gy_res->SetMarkerStyle(2);
    m_gy_res->SetMarkerSize(0.4);
    m_gy_res->Draw("ap");

    c->cd(3);
    m_gz_res->SetTitle("Z residuals;x [cm];#Deltaz [cm]");
    m_gz_res->SetMarkerStyle(2);
    m_gz_res->SetMarkerSize(0.4);
    m_gz_res->Draw("ap");

    c->cd(4);
    m_gr_res->SetTitle("Residuals;x [cm];distance [cm]");
    m_gr_res->SetMarkerStyle(2);
    m_gr_res->SetMarkerSize(0.4);
    m_gr_res->Draw("ap");

    c->Write();
  }
};

/// @brief Class for plotting residues of reconstructed interaction points.
class HistResTask : public RecoTask
{
  TH1F *m_hx_res, *m_hy_res, *m_hz_res, *m_hr_res;

public:
  void PreElectronLoop() override
  {
    m_hx_res = new TH1F("hx_res", "X residuals;x deviation [cm];# of electrons", 25, -0.2, 0.2);
    m_hy_res = new TH1F("hy_res", "Y residuals;y deviation [cm];# of electrons", 25, -0.2, 0.2);
    m_hz_res = new TH1F("hz_res", "Z residuals;z deviation [cm];# of electrons", 25, -0.2, 0.2);
    m_hr_res = new TH1F("hr_res", "Residuals;Deviation [cm];# of electrons", 25, 0, 0.25);
  }

  void ElectronLoop() override
  {
    X17::TrackLoop* loop = GetLoop_();

    MicroPoint micro = loop->curr_micro;
    RecoPoint reco   = loop->curr_reco;

    m_hx_res->Fill(reco.x() - micro.x0());
    m_hy_res->Fill(reco.y() - micro.y0());
    m_hz_res->Fill(reco.z() - micro.z0());
    m_hr_res->Fill(
      std::sqrt(pow((reco.x() - micro.x0()), 2) + pow((reco.y() - micro.y0()), 2) + pow((reco.z() - micro.z0()), 2)));
  }

  void PostElectronLoop() override
  {
    X17::TrackLoop* loop = GetLoop_();

    std::string c_name = "c_fit_res2" + std::to_string(loop->curr_track_index);
    TCanvas* c         = new TCanvas(c_name.c_str(), "Reconstruction residuals");
    c->Divide(2, 2);
    c->cd(1);
    m_hx_res->Draw();
    c->cd(2);
    m_hy_res->Draw();
    c->cd(3);
    m_hz_res->Draw();
    c->cd(4);
    m_hr_res->Draw();

    c->Write();
  }
};

/// @brief Class for plotting Runge-Kutta simulated tracks with lower than selected energy.
class PlotSelectionTask : public RecoTask
{
  std::vector<TGraph2D*> m_tracks;
  double m_e_max;
  // CircleFit3D* cfit = nullptr;
  RKFitCircleTask* m_cfit_energy = nullptr;

public:
  void PreTrackLoop() override
  {
    // cfit = new CircleFit3D();
    // cfit->SetFitter();
  }

  void PreElectronLoop() override
  {
    const TrackRK* track = GetLoop_()->curr_rk;

    // cfit = new CircleFit3D(track->origin,track->orientation);
    m_tracks.push_back(new TGraph2D());
  }

  void ElectronLoop() override
  {
    X17::TrackLoop* loop = GetLoop_();

    const TrackRK* track = loop->curr_rk;
    RKPoint p            = loop->curr_rkpoint;

    // cfit->AddPoint(p);
    m_tracks.back()->AddPoint(p.x(), p.y(), p.z());
  }

  void PostElectronLoop() override
  {
    // cfit->FitCircle3D();

    if (m_cfit_energy->m_cfit->GetEnergy(*GetLoop_()->magfield) > m_e_max)
      m_tracks.pop_back();
  }

  void PostTrackLoop() override
  {
    using namespace X17::constants;
    TCanvas* c = new TCanvas("c_cfit_failed", "Tracks with failed circle fit");

    // A histogram for scalling of the axes.
    TH3F* scale = new TH3F("scale", "Tracks with failed circle fit;x [cm];y [cm];z [cm]", 1, xmin, xmax, 1, -yhigh,
                           yhigh, 1, zmin, zmax);
    scale->Draw("");
    gStyle->SetOptStat(0);
    scale->GetXaxis()->SetTitleOffset(1.5);
    scale->GetZaxis()->SetTitleOffset(1.5);

    for (TGraph2D* g : m_tracks)
    {
      g->SetLineColor(kRed);
      g->Draw("LINE same");
    }

    c->Write();

    std::cout << "PlotSelectionTask: failed fits: " << m_tracks.size() << "\n";
  }

public:
  PlotSelectionTask(RKFitCircleTask* cfit, double E_max = 3.5e+6)
    : m_e_max(E_max), m_cfit_energy(cfit)
  {
  }
};

/// @brief Task for plotting several Runge-Kutta tracks with a color palette.
class PlotForwardTask : public RecoTask
{
  std::vector<TLine*> m_track_segments;
  double m_x_prev = X17::constants::xmin;
  double m_z_prev = 0;
  Color_t m_color;

public:
  void PreElectronLoop() override
  {
    X17::TrackLoop* loop = GetLoop_();

    double norm_energy = (loop->curr_rk->kin_energy - 3e+6) / 10e+6;

    int color_index = static_cast<int>(norm_energy * (gStyle->GetNumberOfColors() - 1));
    m_color         = static_cast<Color_t>(gStyle->GetColorPalette(color_index));
  }

  void ElectronLoop() override
  {
    X17::TrackLoop* loop = GetLoop_();

    const TrackRK* track = loop->curr_rk;
    RKPoint p            = loop->curr_rkpoint;
    TLine* line          = new TLine(m_x_prev, m_z_prev, p.x(), p.z());
    line->SetLineColor(m_color);
    m_x_prev = p.x();
    m_z_prev = p.z();
    if (p.x() < 15)
      m_track_segments.push_back(line);
  }

  void PostElectronLoop() override
  {
    m_x_prev = X17::constants::xmin;
    m_z_prev = 0;
  }

  void PostTrackLoop() override
  {
    using namespace X17::constants;

    int c_width  = iround(1.5 * 500 / 0.72);
    int c_height = iround(1.5 * 500 / 7. * 8.5 / 0.8);

    TCanvas* c = new TCanvas("c_forward", "Forward tracks with different energies", c_width, c_height);
    c->SetTopMargin(0.07);
    c->SetBottomMargin(0.13);
    c->SetRightMargin(0.2);
    TH2F* scale = new TH2F("scale", ";x [cm];z [cm];E [MeV]", 10, xmin, 15, 10, -6, 1);
    scale->SetStats(false);
    scale->SetBinContent(1, 3);
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

    for (TLine* l : m_track_segments)
    {
      l->SetLineWidth(2);
      l->Draw("same");
    }
    c->Write();
  }
};