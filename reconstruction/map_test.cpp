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
#include "RecoTasks.h" // IWYU pragma: keep
#include "TrackLoop.h"
#include "Utilities.h"
#include "Vector.h"

int main()
{
  // Loading the magnetic field data.
  X17::Field<X17::Vector>* magfield
    = X17::LoadField("../../data/elmag/VecB2.txt", { -20, -30, -30 }, { 20, 30, 30 }, 0.5);

  // Loading the ionization electron drift map.
  TFile* map_input               = new TFile("../../data/ion_map/sample_2.0/map.root");
  const X17::DriftMap* const map = reinterpret_cast<X17::DriftMap*>(map_input->Get("map"));

  // Loading file(s) with microscopic tracks.
  std::string micro_tracks_folder;
  TChain* micro_tracks = new TChain("tracks_small");

  micro_tracks_folder = "../../data/micro_tracks/grid_01/";
  AddFilesToTChain(micro_tracks, micro_tracks_folder + "tracks_small", ".root", 1, 2000);
  micro_tracks_folder = "../../data/micro_tracks/grid_02/";
  AddFilesToTChain(micro_tracks, micro_tracks_folder + "tracks_small", ".root", 1, 9702);

  std::cout << "Processing " << micro_tracks->GetEntries() << " tracks.\n";

  // TrackLoop for multiple microscopic tracks.
  TrackLoop* multi_loop        = new TrackLoop(*map, magfield);
  multi_loop->make_track_plots = false;
  // multi_loop->AddTask(new MapRecoCompareTask("c_oldnew_res",true,true));
  multi_loop->AddTask(new RecoPadsTask());                          // NOLINT(misc-include-cleaner)
  multi_loop->AddTask(new MapRecoTask("c_fit_res", true, true));    // NOLINT(misc-include-cleaner)
  multi_loop->AddTask(new MapRecoTask("c_e_fit_res", true, false)); // NOLINT(misc-include-cleaner)
  multi_loop->AddTask(new MapRecoTask("c_p_fit_res", false, true)); // NOLINT(misc-include-cleaner)

  gErrorIgnoreLevel = 6001;

  TFile* out_file
    = new TFile((micro_tracks_folder + "../map_test.root").c_str(), "RECREATE", "Tracks from microscopic simulation");

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