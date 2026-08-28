// C++ dependencies
#include <exception>
#include <memory>
#include <string>
#include <vector>

// ROOT dependencies
#include "RtypesCore.h"
#include "TRandom3.h"
#include "TTree.h"

// X17 dependencies
#include "AppManager.h"
#include "Field.h"
#include "Logger.h"
#include "Points.h" // NOLINT(readability-duplicate-include)
#include "Track.h"
#include "Vector.h"
#include "X17Utilities.h"

int main(int /*argc*/, char** /*argv*/)
{
  X17::AppManager man("mcs_tracks", 3, "Drift of MCS tracks using the map.");

  bool driftSecondary = false;

  if (driftSecondary)
    X17::Logger::Get().Essential("Drifting secondary particles (delta rays).");
  else
    X17::Logger::Get().Essential("Drifting primary particles only.");

  TTree* sim_tree = man.LoadTreeFromFile("../source/X17geant4/X17geant4/mcs_3MeV.root", "gasIonization");

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

  TTree* out_tree = man.CreateOutputTree("data/mcs_tracks/mcs_tracks4.root", "tracks_small",
                                         "Tree of tracks with map-drifted points");

  std::unique_ptr<TRandom3> rand = man.CreateRNG(0);

  X17::TrackMicro track;
  out_tree->Branch("track_small", &track);
  std::vector<X17::MicroPoint> points;

  bool electron           = true;
  X17::Vector orientation = { 1, 0, 0 };
  X17::Vector origin      = { 6.51, 0, 0 };
  double kin_energy       = 8.0E+6 - X17::constants::E0;

  X17::Logger::Get().Essential("Drifting " + std::to_string(sim_tree->GetEntries()) + " events.");

  for (Long64_t i = 0; i < sim_tree->GetEntries(); i++)
  {
    X17::Logger::Get().ProgressBar(i, sim_tree->GetEntries());

    sim_tree->GetEntry(i);

    if (!driftSecondary && parentID != 1)
      continue;

    start.point /= 10; // mm -> cm

    X17::EndPoint end;
    try
    {
      end = map->GetField(start.point).GetRandomPoint(rand.get());
    }
    catch (const std::exception& e)
    {
      X17::Logger::Get().Fatal(e.what());
    }

    end.t += start.t;

    X17::MicroPoint micro;
    micro.start = start;
    micro.end   = end;

    if (prevID != eventID)
    {
      if (prevID != -1)
      {
        track = X17::TrackMicro(electron, points, origin, orientation, kin_energy, {});
        out_tree->Fill();
      }
      points.clear();
      prevID = eventID;
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