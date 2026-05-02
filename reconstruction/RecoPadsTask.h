#pragma once

// ROOT dependencies
#include <algorithm>

#include "TH2Poly.h"
#include "TLine.h"
#include "TStyle.h"

// X17 dependencies
#include "PadLayout.h"
#include "Reconstruction.h"
#include "TrackLoop.h"
#include "X17Utilities.h"

/// @brief Class for padded reconstruction.
class RecoPadsTask : public X17::RecoTask
{
  friend class MicroCircleAndRKFitTask;
  friend class MicroFitAndSaveTask;

  static constexpr int s_timebins = 200;
  int m_padhits[X17::constants::channels][s_timebins];
  std::vector<X17::RecoPoint> m_reco_data;
  TGraph2D *m_g_xyz, *m_g_xyz_reco;
  TH2Poly* m_pads;
  TH2Poly* m_e_pads;
  TH2Poly* m_p_pads;
  TH3F* m_scale;
  TCanvas* m_c_reco;

  double m_height; // Height at which the pads will be drawn.

public:
  void PreTrackLoop() override
  {
    for (TH2Poly** poly : { &m_pads, &m_e_pads, &m_p_pads })
    {
      *poly = new TH2Poly();

      for (int i = 1; i <= X17::constants::channels; i++)
      {
        double x1, y1, x2, y2;
        X17::DefaultLayout::GetDefaultLayout().GetPadCorners(i, x1, y1, x2, y2, true);
        (*poly)->AddBin(y1, x1, y2, x2);
      }
    }
  }

  void PreElectronLoop() override
  {
    X17::TrackLoop* loop = GetLoop_();
    for (auto& m_padhit : m_padhits)
      for (int& j : m_padhit)
        j = 0;

    if (loop->make_track_plots)
    {
      m_g_xyz      = new TGraph2D();
      m_g_xyz_reco = new TGraph2D();
    }
  }

  void ElectronLoop() override
  {
    X17::TrackLoop* loop  = GetLoop_();
    X17::MicroPoint micro = loop->curr_micro;

    int channel = X17::DefaultLayout::GetDefaultLayout().GetPad(micro.x1(), micro.y1());
    int timebin = iround(micro.t1() / 100.0);

    if (timebin > s_timebins - 1)
      std::cerr << "ERROR: Invalid timebin: " << timebin << '\n';
    if (channel == -1)
      std::cerr << "ERROR: No pad hit found. Coordinates x,y: " << micro.x1() << ", " << micro.y1() << '\n';

    m_padhits[channel - 1][timebin]++;

    if (loop->make_track_plots)
      m_g_xyz->AddPoint(micro.x0(), micro.y0(), micro.z0());
  }

  void PostElectronLoop() override
  {
    using namespace X17::constants;

    X17::TrackLoop* loop = GetLoop_();

    // Reconstruction with pads.
    int max_charge   = 0;
    X17::Vector vmin = { 100, 100, 100 };
    X17::Vector vmax = { -100, -100, -100 };
    std::vector<TPolyLine3D*> marker_lines;

    for (auto& m_padhit : m_padhits)
      for (int j : m_padhit)
        max_charge = std::max(j, max_charge);

    // std::cout << "Max charge: " << max_charge << std::endl;

    for (int i = 0; i < channels; i++)
    {
      // ReportProgress(i,channels);
      for (int j = 0; j < s_timebins; j++)
      {
        int charge = m_padhits[i][j];
        if (charge == 0)
          continue;

        double time = 100 * j + 50;
        double xpad, ypad;
        X17::DefaultLayout::GetDefaultLayout().GetPadCenter(i + 1, xpad, ypad);

        X17::TrackLoop* loop = GetLoop_();
        X17::RecoPoint reco  = Reconstruct(loop->map, X17::EndPoint(xpad, ypad, zmax, time));
        m_reco_data.emplace_back(reco.x(), reco.y(), reco.z(), charge);
        if (loop->make_track_plots)
          m_g_xyz_reco->AddPoint(reco.x(), reco.y(), reco.z());

        // Drawing markers.
        bool old_reco      = false;
        const int pts      = old_reco ? 4 : 100;
        double norm_charge = static_cast<double>(charge) / max_charge;
        int color_index    = static_cast<int>(norm_charge * (gStyle->GetNumberOfColors() - 1));
        Color_t color      = static_cast<Color_t>(gStyle->GetColorPalette(color_index));

        double x1, y1, x2, y2;
        X17::DefaultLayout::GetDefaultLayout().GetPadCorners(i + 1, x1, y1, x2, y2);
        vmin.x = std::min(vmin.x, x1);
        vmin.y = std::min(vmin.y, y1);
        vmax.x = std::max(vmax.x, x2);
        vmax.y = std::max(vmax.y, y2);

        X17::Vector corners[4] = { { x1, y1, -8 }, { x1, y2, -8 }, { x2, y1, -8 }, { x2, y2, -8 } };

        // Corner lines
        for (auto v_corner : corners)
        {
          TPolyLine3D* corner_line = new TPolyLine3D(pts);
          corner_line->SetLineColor(color);

          for (int it = 0; it < pts; it++)
          {
            double t = time - 50 + 100.0 * it / (pts - 1.0);
            X17::EndPoint t_corner(v_corner, t);
            X17::RecoPoint reco2   = old_reco ? X17::ReconstructOld(loop->map, t_corner, 1e-6, false)
                                              : X17::Reconstruct(loop->map, t_corner);
            X17::Vector line_point = reco.AsVector() + norm_charge * (reco2.AsVector() - reco.AsVector());
            corner_line->SetPoint(it, line_point.x, line_point.y, line_point.z);

            vmin.z = std::min(line_point.z, vmin.z);
            vmax.z = std::max(line_point.z, vmax.z);
          }

          marker_lines.push_back(corner_line);
        }

        // Plane lines
        for (double t : { time - 50, time + 50 })
        {
          TPolyLine3D* lines[4]
            = { new TPolyLine3D(pts), new TPolyLine3D(pts), new TPolyLine3D(pts), new TPolyLine3D(pts) };

          for (auto* line : lines)
          {
            line->SetLineColor(color);
            marker_lines.push_back(line);
          }

          for (int k = 0; k < pts; k++)
          {
            std::vector<X17::RecoPoint> points;
            if (old_reco)
            {
              points = {
                X17::ReconstructOld(loop->map, { x1 + k * (x2 - x1) / (pts - 1), y1, -8, t }, 1e-6, false),
                X17::ReconstructOld(loop->map, { x2, y1 + k * (y2 - y1) / (pts - 1), -8, t }, 1e-6, false),
                X17::ReconstructOld(loop->map, { x1 + k * (x2 - x1) / (pts - 1), y2, -8, t }, 1e-6, false),
                X17::ReconstructOld(loop->map, { x1, y1 + k * (y2 - y1) / (pts - 1), -8, t }, 1e-6, false),
              };
            }

            else
            {
              points = {
                X17::Reconstruct(loop->map, { x1 + k * (x2 - x1) / (pts - 1), y1, -8, t }),
                X17::Reconstruct(loop->map, { x2, y1 + k * (y2 - y1) / (pts - 1), -8, t }),
                X17::Reconstruct(loop->map, { x1 + k * (x2 - x1) / (pts - 1), y2, -8, t }),
                X17::Reconstruct(loop->map, { x1, y1 + k * (y2 - y1) / (pts - 1), -8, t }),
              };
            }

            for (int l = 0; l < 4; l++)
            {
              X17::Vector line_point = reco.AsVector() + norm_charge * (points[l].AsVector() - reco.AsVector());
              lines[l]->SetPoint(k, line_point.x, line_point.y, line_point.z);
            }
          }
        }
      }
    }

    // Margin.
    vmin -= { 0.2, 0.2, 0.2 };
    vmax += { 0.2, 0.2, 0.2 };

    // Old code for drawing markers.
    // auto reco_markers = GetDataMarkers(reco_data);
    // for (auto m : reco_markers) m->Draw("same");

    if (loop->make_track_plots)
    {
      std::string c_reco_name = "c_track_xyz_" + std::to_string(loop->curr_track_index);
      m_c_reco = new TCanvas(c_reco_name.c_str(), "Electron track reconstruction with pads and time bins",
                             iround(700 * 1.33 * 1.5), iround(500 * 1.5));
      ApplyThesisStyle(m_c_reco);
      m_c_reco->SetRightMargin(0.15);

      // A histogram for scalling of the axes.
      m_scale = new TH3F("scale", ";x [cm];y [cm];z [cm]", 1, vmin.x, vmax.x, 1, vmin.y, vmax.y, 1, vmin.z, vmax.z);
      ApplyThesisStyle(m_scale);
      m_scale->SetStats(false);
      m_scale->SetBinContent(1, 0);
      m_scale->SetMinimum(0);
      m_scale->SetMaximum(max_charge);
      m_scale->GetXaxis()->SetTitleOffset(1.1);
      m_scale->GetYaxis()->SetTitleOffset(1.15);
      m_scale->GetZaxis()->SetTitleOffset(0.8);
      m_scale->GetYaxis()->SetNdivisions(505);
      m_scale->Draw("BOX2Z");
      // TPaletteAxis* palette = (TPaletteAxis*)scale->GetListOfFunctions()->FindObject("palette");
      // palette->SetX1NDC(0.9);
      // palette->SetX2NDC(0.92);

      m_g_xyz_reco->SetTitle("Electron track reconstruction;x [cm];y [cm];z [cm]");
      m_g_xyz_reco->SetMarkerStyle(2);
      m_g_xyz_reco->SetMarkerSize(0.4);
      m_g_xyz_reco->Draw("p same");

      m_g_xyz->SetMarkerColor(2);
      m_g_xyz->SetMarkerStyle(7);
      m_g_xyz->SetMarkerSize(1.2);
      m_g_xyz->Draw("p same");

      for (auto* line : marker_lines)
        line->Draw("same");

      m_c_reco->Write();

      // delete c_reco;
      // delete scale;

      // DefaultLayout::GetDefaultLayout().DrawPads3D(height);
    }

    if (!loop->make_track_plots)
      for (auto* line : marker_lines)
        delete line;

    for (int i = 1; i <= X17::constants::channels; i++)
      for (int j = 0; j < s_timebins; j++)
      {
        double curr_el_charge = m_pads->GetBinContent(i);
        m_pads->SetBinContent(i, curr_el_charge + m_padhits[i - 1][j]);
        if (loop->curr_microtrack->electron)
        {
          curr_el_charge = m_e_pads->GetBinContent(i);
          m_e_pads->SetBinContent(i, curr_el_charge + m_padhits[i - 1][j]);
        }
        else
        {
          curr_el_charge = m_p_pads->GetBinContent(i);
          m_p_pads->SetBinContent(i, curr_el_charge + m_padhits[i - 1][j]);
        }
      }
  }

  void PostTrackLoop() override
  {
    TH2Poly* polys[3]    = { m_pads, m_e_pads, m_p_pads };
    std::string names[3] = { "pads", "e_pads", "p_pads" };
    for (int i = 0; i < 3; i++)
    {
      TCanvas* c = new TCanvas(names[i].c_str(), names[i].c_str(), 700, 500);
      polys[i]->SetTitle(";y [cm];x [cm]");
      ApplyThesisStyle(c);
      polys[i]->Draw();
      // Phi angles
      int angle_bins = 21;
      double phi_max = atan((X17::constants::win_width / 2) / X17::constants::xmin); // The maximal simulated phi [rad].
      double phi_min = -phi_max;                                                     // The minimal simulated phi [rad].

      for (int i = 0; i < angle_bins; i++)
      {
        using namespace X17::constants;
        double phi  = phi_min + (phi_max - phi_min) * i / (angle_bins - 1);
        TLine* line = new TLine(xmin * tan(phi), xmin, xmax * tan(phi), xmax);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        line->Draw("same");
      }
      c->Write();
    }
  }

public:
  RecoPadsTask(double pad_height = -2.5)
    : m_height(pad_height)
  {
  }
};