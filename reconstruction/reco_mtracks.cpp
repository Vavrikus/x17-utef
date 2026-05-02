// C++ dependencies
#include <iostream>
#include <string>

// ROOT dependencies
#include "TChain.h"
#include "TError.h"
#include "TFile.h"
#include "TTree.h"

// X17 dependencies
#include "Field.h"
#include "RecoTasks.h"
#include "TrackLoop.h"
#include "Utilities.h"
#include "Vector.h"

namespace
{
  int reco_track()
  {
    // Which files to choose
    bool allTracks = true;

    // Loading the magnetic field data.
    X17::Field<X17::Vector>* magfield
      = X17::LoadField("../../data/elmag/VecB2.txt", { -20, -30, -30 }, { 20, 30, 30 }, 0.5);

    // Loading the ionization electron drift map.
    TFile* map_input               = new TFile("../../data/ion_map/sample_2.0/map.root");
    const X17::DriftMap* const map = reinterpret_cast<X17::DriftMap*>(map_input->Get("map"));

    // Loading file(s) with microscopic tracks.
    std::string micro_tracks_folder;
    TChain* micro_tracks = new TChain("tracks_small");

    if (allTracks)
    {
      micro_tracks_folder = "../../data/micro_tracks/grid_01/";
      AddFilesToTChain(micro_tracks, micro_tracks_folder + "tracks_small", ".root", 1, 2000);
      micro_tracks_folder = "../../data/micro_tracks/grid_02/";
      AddFilesToTChain(micro_tracks, micro_tracks_folder + "tracks_small", ".root", 1, 9702);
    }
    else
    {
      micro_tracks_folder = "../../data/micro_tracks/grid_01/";
      micro_tracks->Add((micro_tracks_folder + "tracks_small1000.root").c_str());
    }

    std::cout << "Processing " << micro_tracks->GetEntries() << " tracks.\n";

    // TrackLoop for multiple microscopic tracks.
    TrackLoop* multi_loop        = new TrackLoop(*map, magfield);
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
    if (allTracks)
      out_file = new TFile((micro_tracks_folder + "../reco_tracks.root").c_str(), "RECREATE",
                           "Tracks from microscopic simulation");
    else
      out_file = new TFile((micro_tracks_folder + "track_plots1000.root").c_str(), "RECREATE",
                           "Tracks from microscopic simulation");

    multi_loop->ProcessMulti(micro_tracks);

    out_file->Close();
    map_input->Close();

    delete out_file;
    delete map_input;

    delete magfield;
    delete map;
    delete micro_tracks;
    delete multi_loop;

    return 0;
  }
} // namespace

int main()
{
  return reco_track();
}