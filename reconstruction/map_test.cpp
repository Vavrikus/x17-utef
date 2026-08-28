// C++ dependencies
#include <iostream>
#include <memory>
#include <string>

// ROOT dependencies
#include "TChain.h"
#include "TError.h"
#include "TFile.h"
#include "TTree.h"

// X17 dependencies
#include "AppManager.h"
#include "Field.h"
#include "RecoTasks.h" // IWYU pragma: keep
#include "TrackLoop.h"
#include "Utilities.h"

int main()
{
  X17::AppManager man("map_test", 2, "Testing the drift map reconstruction.");

  // Loading the magnetic field data.
  std::unique_ptr<X17::MagField> magfield = X17::AppManager::LoadMagField();

  // Loading the ionization electron drift map.
  std::unique_ptr<X17::DriftMap> map = man.LoadMap("2.0");

  // Loading file(s) with microscopic tracks.
  std::string micro_tracks_folder;
  TChain* micro_tracks = new TChain("tracks_small");

  micro_tracks_folder = "data/micro_tracks/grid_01/";
  AddFilesToTChain(micro_tracks, micro_tracks_folder + "tracks_small", ".root", 1, 2000);
  micro_tracks_folder = "data/micro_tracks/grid_02/";
  AddFilesToTChain(micro_tracks, micro_tracks_folder + "tracks_small", ".root", 1, 9702);

  std::cout << "Processing " << micro_tracks->GetEntries() << " tracks.\n";

  // TrackLoop for multiple microscopic tracks.
  TrackLoop* multi_loop        = new TrackLoop(*map, magfield.get());
  multi_loop->make_track_plots = false;
  // multi_loop->AddTask(new MapRecoCompareTask("c_oldnew_res",true,true));
  multi_loop->AddTask(new RecoPadsTask());                          // NOLINT(misc-include-cleaner)
  multi_loop->AddTask(new MapRecoTask("c_fit_res", true, true));    // NOLINT(misc-include-cleaner)
  multi_loop->AddTask(new MapRecoTask("c_e_fit_res", true, false)); // NOLINT(misc-include-cleaner)
  multi_loop->AddTask(new MapRecoTask("c_p_fit_res", false, true)); // NOLINT(misc-include-cleaner)
  multi_loop->AddTask(new EdepTask());

  gErrorIgnoreLevel = 6001;

  TFile* out_file
    = new TFile((micro_tracks_folder + "../map_test.root").c_str(), "RECREATE", "Tracks from microscopic simulation");

  multi_loop->ProcessMulti(micro_tracks);

  out_file->Close();

  delete out_file;

  delete micro_tracks;
  delete multi_loop;

  return 0;
}