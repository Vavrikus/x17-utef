// C++ dependencies
#include <memory>
#include <vector>

// ROOT dependencies
#include "Rtypes.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TH3.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TTree.h"

// X17 dependencies
#include "AppManager.h"
#include "Field.h"
#include "Matrix.h"
#include "PadLayout.h"
#include "Points.h"
#include "RK4.h"
#include "Track.h"
#include "Utilities.h"
#include "Vector.h"
#include "X17Utilities.h"

using namespace X17::constants;

namespace
{
  int rk_tracks()
  {
    X17::AppManager man("rk_tracks", 3, "Simulate tracks using Runge-Kutta method.");

    constexpr int n_tracks = 25;            // The number of tracks to be simulated by Runge-Kutta.
    constexpr double step  = 1E-13 * 16.66; // The step of Runge-Kutta [s] (original times 8 MeV gamma factor).

    // Loading the magnetic field.

    std::unique_ptr<X17::MagField> magfield = X17::AppManager::LoadMagField();

    // Some necessary variables for simulating and saving track parameters.
    TTree* simulated_tracks = new TTree("rk_tracks", "Runge-Kutta simulated tracks");
    TRandom3* rand          = new TRandom3(0);

    X17::TrackRK current_track;

    simulated_tracks->Branch("track", &current_track);

    std::vector<TGraph2D*> track_graphs;

    // The loop for the track simulation.
    for (int i = 0; i < n_tracks; i++)
    {
      // if ((100 * i) % n_tracks == 0) std::cout << "Progress: " << 100*i/n_tracks << " \%\n";
      ReportProgress(i, n_tracks);

      // Generating random initial parameters.
      bool electron;
      X17::Vector origin, orientation;
      double kin_en;

      // X17::GetRandomTrackParams(rand,electron,origin,orientation,kin_en);
      electron    = true;
      origin      = { X17::constants::xmin, 0, X17::constants::win_height * (-1. / 2. + i * 1. / n_tracks) };
      orientation = { 1, 0, 0 };
      kin_en      = 1e+6; // 3e+6 + i * 1.0 / n_tracks * (13e+6 - 3e+6);

      // The actual track simulation.
      X17::RK4<8>* track = GetTrackRK(*magfield, electron, step, kin_en, origin, orientation);
      track->Integrate();

      std::vector<X17::Matrix<8, 1>> results = track->GetResults();
      std::vector<X17::RKPoint> points;

      using namespace X17::constants;
      points.reserve(results.size());
      for (auto r : results)
        points.emplace_back(m2cm * r.at(1, 0), m2cm * r.at(2, 0), m2cm * r.at(3, 0), 1e+9 / c * r.at(0, 0));

      current_track = X17::TrackRK(electron, points, origin, orientation, kin_en);
      simulated_tracks->Fill();

      if ((100 * i) % n_tracks == 0)
        track_graphs.push_back(GetGraphRK(track));
    }

    // Plotting some of the tracks.
    double height = 8;
    using namespace X17::constants;
    TCanvas* c_tracks = new TCanvas("c_tracks", "Example tracks");
    // A histogram for scalling of the axes.
    TH3F* scale
      = new TH3F("scale", "Example tracks;x [cm];y [cm];z [cm]", 1, xmin, xmax, 1, -yhigh, yhigh, 1, -height, height);
    scale->Draw("");
    gStyle->SetOptStat(0);
    scale->GetXaxis()->SetTitleOffset(1.5);
    scale->GetZaxis()->SetTitleOffset(1.5);
    for (auto* g : track_graphs)
    {
      g->SetLineColor(kRed);
      g->Draw("LINE same");
    }

    X17::DefaultLayout::GetDefaultLayout().DrawPads3D(height);

    TFile* outfile = new TFile("../../../data/rk_tracks/rk_tracks_forward2.root", "RECREATE");
    simulated_tracks->Write();
    c_tracks->Write();

    // Freeing memory.
    delete simulated_tracks;
    delete rand;
    delete c_tracks;
    delete scale;
    delete outfile;

    return 0;
  }
} // namespace

int main()
{
  return rk_tracks();
}