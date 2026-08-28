// C++ dependencies
#include <memory>
#include <string>

// ROOT dependencies
#include "TFile.h"
#include "TTree.h"

// X17 dependencies
#include "AppManager.h"
#include "Field.h"
#include "RecoTasks.h"
#include "Track.h"
#include "TrackLoop.h"

int main()
{
  X17::AppManager man("reco_test", 2, "Testing on smaller samples and Runge-Kutta tracks.");

  // Loading the magnetic field data.
  std::unique_ptr<X17::MagField> magfield = X17::AppManager::LoadMagField();

  // Loading the ionization electron drift map.
  std::unique_ptr<X17::DriftMap> map = man.LoadMap("2.0");

  // Loading file microscopic tracks one by one.
  X17::TrackMicro* track_Evf_100 = nullptr;
  TTree track_selection("track_selection", "Track selection");
  track_selection.Branch("track_small", &track_Evf_100);

  std::string single_track_folder = "data/micro_tracks/grid_01/";

  // 8 MeV, min theta, min phi
  TTree* tracks910 = man.LoadTreeFromFile(single_track_folder + "tracks_small910.root", "tracks_small");
  tracks910->SetBranchAddress("track_small", &track_Evf_100);
  tracks910->GetEntry(1);
  track_selection.Fill();

  // 8 MeV, zero theta, zero phi
  TTree* tracks1000 = man.LoadTreeFromFile(single_track_folder + "tracks_small1000.root", "tracks_small");
  tracks1000->SetBranchAddress("track_small", &track_Evf_100);
  tracks1000->GetEntry(4);
  // track_selection.Fill();

  // 3 MeV, zero theta, zero phi
  TTree* tracks92 = man.LoadTreeFromFile(single_track_folder + "tracks_small92.root", "tracks_small");
  tracks92->SetBranchAddress("track_small", &track_Evf_100);
  tracks92->GetEntry(1);
  // track_selection.Fill();

  // TrackLoop for single microscopic track.
  TrackLoop* loop        = new TrackLoop(*map, magfield.get());
  loop->make_track_plots = true;

  RecoPadsTask* t = new RecoPadsTask(); // NOLINT(misc-include-cleaner)
  loop->AddTask(t);
  // loop->AddTask(new MicroCircleAndRKFitTask(t));

  // Loading file with Runge-Kutta tracks.
  TTree* rk_tracks = man.LoadTreeFromFile("data/rk_tracks/rk_tracks_forward2.root", "rk_tracks");
  // rk_tracks->Print();

  // TrackLoop for Runge-Kutta simulated tracks.
  TrackLoop* rk_loop = new TrackLoop(*map, magfield.get());
  // auto t2 = new RKFitCircleTask();
  // rk_loop->AddTask(t2);
  // rk_loop->AddTask(new PlotSelectionTask(t2));
  rk_loop->AddTask(new PlotForwardTask());

  // Processing.
  TFile out_file("track_plots.root", "RECREATE", "Tracks from microscopic simulation");
  loop->ProcessMulti(&track_selection);
  out_file.Close();

  TFile out_file2("../../data/rk_tracks/rk_plots_forward2.root", "RECREATE", "Tracks from Runge-Kutta simulation fit");
  rk_loop->ProcessRK(rk_tracks);
  out_file2.Close();

  return 0;
}