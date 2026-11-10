// C++ dependencies
#include <exception>
#include <memory>
#include <string>
#include <vector>

// ROOT dependencies
#include "RtypesCore.h"
#include "TChain.h"
#include "TRandom3.h"
#include "TTree.h"

// X17 dependencies
#include "AppManager.h"
#include "Field.h"
#include "Logger.h"
#include "Points.h" // NOLINT(readability-duplicate-include)
#include "Track.h"
#include "Utilities.h"
#include "Vector.h"
#include "X17Utilities.h"

int main(int /*argc*/, char** /*argv*/)
{
  X17::AppManager man("mcs_tracks", 3, "Drift of MCS tracks using the map.");

  bool driftSecondary = false;

  if (driftSecondary)
    X17_LOGGER.Essential("Drift of secondary ionization electrons (delta rays) included.");
  else
    X17_LOGGER.Essential("Drifting primary ionization electrons only.");

  // std::string filename = "../source/X17geant4/X17geant4/mcs_13MeV.root";
  // TTree* sim_tree = man.LoadTreeFromFile(filename, "gasIonization");
  // TTree* primary_tree = man.LoadTreeFromFile(filename, "primParticle");

  TChain* sim_tree              = new TChain("gasIonization");
  TChain* primary_tree          = new TChain("primParticle");
  std::string mcs_tracks_folder = "/storage/projects/utefx17/martin/data/mcs_tracks/new_tracks/grid01/";

  for (TChain* chain : { sim_tree, primary_tree })
    AddFilesToTChain(chain, mcs_tracks_folder + "mcs_tracks_", ".root", 1, 14);

  // sim_tree->Add(filename.c_str());
  // primary_tree->Add(filename.c_str());

  std::unique_ptr<X17::DriftMap> map = man.LoadMap("2.0");

  X17::StartPoint start;
  int parentID;
  int eventID;
  int prevID = -1;

  sim_tree->SetBranchAddress("X", &start.point.x);
  sim_tree->SetBranchAddress("Y", &start.point.y);
  sim_tree->SetBranchAddress("Z", &start.point.z);
  sim_tree->SetBranchAddress("Time", &start.t);
  sim_tree->SetBranchAddress("ParentID", &parentID);
  sim_tree->SetBranchAddress("Event", &eventID);

  int last_event = -1;
  char particleName[16];
  bool electron;
  X17::Vector orientation;
  X17::Vector origin;
  double kin_energy;

  primary_tree->SetBranchAddress("Particle", particleName);
  primary_tree->SetBranchAddress("VtxX", &origin.x);
  primary_tree->SetBranchAddress("VtxY", &origin.y);
  primary_tree->SetBranchAddress("VtxZ", &origin.z);
  primary_tree->SetBranchAddress("Px", &orientation.x);
  primary_tree->SetBranchAddress("Py", &orientation.y);
  primary_tree->SetBranchAddress("Pz", &orientation.z);
  primary_tree->SetBranchAddress("TotalEnergy", &kin_energy);

  TTree* out_tree
    = man.CreateOutputTree("data/mcs_tracks/mcs_tracks.root", "tracks_small", "Tree of tracks with map-drifted points");

  std::unique_ptr<TRandom3> rand = man.CreateRNG(0);

  X17::TrackMicro track;
  out_tree->Branch("track_small", &track);
  out_tree->SetAutoSave(-1e+9); // Write into file every 1 GB.
  std::vector<X17::MicroPoint> points;

  X17_LOGGER.Essential("Drifting ", sim_tree->GetEntries(), " electrons from ", primary_tree->GetEntries(),
                       " primary particles.");

  // Force TChain to calculate file offsets before looping
  sim_tree->LoadTree(0);
  primary_tree->LoadTree(0);

  long long last_global_event = -1;
  long long prev_global_id    = -1;

  for (Long64_t i = 0; i < sim_tree->GetEntries(); i++)
  {
    X17_LOGGER.ProgressBar(i, sim_tree->GetEntries());

    sim_tree->GetEntry(i);

    if (!driftSecondary && parentID != 1)
      continue;

    // 2. Calculate the unique GLOBAL event ID across all files in the TChain
    int current_tree_idx      = sim_tree->GetTreeNumber();
    long long global_event_id = primary_tree->GetTreeOffset()[current_tree_idx] + eventID;

    if (global_event_id != last_global_event)
    {
      // Fetch using the global index, ensuring it reads from the correct file
      primary_tree->GetEntry(global_event_id);
      last_global_event = global_event_id;

      std::string p_name(particleName);
      if (p_name == "e-")
        electron = true;
      else if (p_name == "e+")
        electron = false;
      else
        X17_LOGGER.Fatal("Unexpected primary particle: ", p_name);

      origin /= 10;                     // mm -> cm
      orientation.Normalize();          // it contains momentum prior
      kin_energy *= 1E+6;               // MeV -> eV
      kin_energy -= X17::constants::E0; // contains total energy prior
    }

    start.point /= 10; // mm -> cm

    X17::EndPoint end;
    try
    {
      if (!X17::IsInTPC(start.point))
      {
        X17_LOGGER.Warning("Ionization outside of TPC: ", start.point.x, ", ", start.point.y, ", ", start.point.z);
        continue;
      }

      end = map->GetField(start.point).GetRandomPoint(rand.get());
    }
    catch (const std::exception& e)
    {
      X17_LOGGER.Fatal(e.what());
    }

    end.t += start.t;

    X17::MicroPoint micro;
    micro.start = start;
    micro.end   = end;

    if (prev_global_id != global_event_id)
    {
      if (global_event_id != prev_global_id + 1)
        X17_LOGGER.Warning("Event IDs out of order: ", prev_global_id, " -> ", global_event_id);
      if (prev_global_id != -1)
      {
        track = X17::TrackMicro(electron, points, origin, orientation, kin_energy, {});
        out_tree->Fill();
        points.clear();
      }
      prev_global_id = global_event_id;
    }

    points.push_back(micro);
  }

  if (!points.empty())
  {
    track = X17::TrackMicro(electron, points, origin, orientation, kin_energy, {});
    out_tree->Fill();
  }

  return 0;
}