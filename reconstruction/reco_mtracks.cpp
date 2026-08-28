// C++ dependencies
#include <iostream>
#include <memory>
// #include <string>

// ROOT dependencies
#include "TChain.h"
#include "TError.h"
#include "TFile.h"
#include "TTree.h"

// X17 dependencies
#include "AppManager.h"
#include "Field.h"
#include "RecoTasks.h"
#include "TrackLoop.h"
// #include "Utilities.h"

namespace
{
  int reco_mtracks()
  {
    // Setting up the application.
    X17::AppManager man("reco_mtracks", 2, "Reconstruction of microscopic tracks.");

    // Which files to choose
    bool allTracks = true;

    // Loading the magnetic field data.
    std::unique_ptr<X17::MagField> magfield = X17::AppManager::LoadMagField();

    // Loading the ionization electron drift map.
    std::unique_ptr<X17::DriftMap> map = man.LoadMap("2.0");

    // Loading file(s) with microscopic tracks.
    // std::string micro_tracks_folder;
    TChain* micro_tracks = new TChain("tracks_small");

    // if (allTracks)
    // {
    //   micro_tracks_folder = "../../data/micro_tracks/grid_01/";
    //   AddFilesToTChain(micro_tracks, micro_tracks_folder + "tracks_small", ".root", 1, 2000);
    //   micro_tracks_folder = "../../data/micro_tracks/grid_02/";
    //   AddFilesToTChain(micro_tracks, micro_tracks_folder + "tracks_small", ".root", 1, 9702);
    // }
    // else
    // {
    //   micro_tracks_folder = "../../data/micro_tracks/grid_01/";
    //   micro_tracks->Add((micro_tracks_folder + "tracks_small1000.root").c_str());
    // }

    micro_tracks->Add("../../data/mcs_tracks/mcs_tracks3.root");

    std::cout << "Processing " << micro_tracks->GetEntries() << " tracks.\n";

    // TrackLoop for multiple microscopic tracks.
    TrackLoop* multi_loop        = new TrackLoop(*map, magfield.get());
    multi_loop->make_track_plots = !allTracks;

    if (!allTracks)
    {
      // multi_loop->AddTask(new DriftTimeTask());
      multi_loop->AddTask(new XZPlotTask());
      multi_loop->AddTask(new XYPlotTask());
      multi_loop->AddTask(new GraphResTask());
      multi_loop->AddTask(new HistResTask());
    }

    RecoPadsTask* t2 = new RecoPadsTask(-1.6); // NOLINT(misc-include-cleaner)
    multi_loop->AddTask(t2);
    if (allTracks)
      multi_loop->AddTask(new MicroFitAndSaveTask(t2)); // NOLINT(misc-include-cleaner)
    else
      multi_loop->AddTask(new MicroCircleAndRKFitTask(t2)); // NOLINT(misc-include-cleaner)

    gErrorIgnoreLevel = 6001;

    TFile* out_file = nullptr;
    // if (allTracks)
    //   out_file = new TFile((micro_tracks_folder + "../reco_tracks.root").c_str(), "RECREATE",
    //                        "Tracks from microscopic simulation");
    // else
    //   out_file = new TFile((micro_tracks_folder + "track_plots1000.root").c_str(), "RECREATE",
    //                        "Tracks from microscopic simulation");

    out_file = new TFile("../../data/mcs_tracks/reco_tracks3.root", "RECREATE", "Tracks from map-drifted simulation");

    multi_loop->ProcessMulti(micro_tracks);

    out_file->Close();

    delete out_file;

    delete micro_tracks;
    delete multi_loop;

    return 0;
  }
} // namespace

int main()
{
  return reco_mtracks();
}