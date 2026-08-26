// C++ dependencies
#include <iostream>
#include <memory>
#include <random>
#include <vector>

// ROOT dependencies
#include "RtypesCore.h"
#include "TFile.h"
#include "TParameter.h"
#include "TRandom3.h"
#include "TTree.h"

// X17 dependencies
#include "Field.h"
#include "Points.h" // NOLINT(readability-duplicate-include)
#include "Track.h"
#include "Utilities.h"
#include "Vector.h"
#include "X17Utilities.h"

int main(int /*argc*/, char** /*argv*/)
{
  bool driftSecondary = false;

  // Setting the random seed
  UInt_t seed = std::random_device{}();
  std::cout << "Setting random seed: " << seed << '\n';
  std::unique_ptr<TRandom3> rand = std::make_unique<TRandom3>(seed);

  TFile* map_input               = TFile::Open("../../../data/ion_map/sample_2.0/map.root", "READ");
  const X17::DriftMap* const map = reinterpret_cast<X17::DriftMap*>(map_input->Get("map"));

  TFile* sim_input = TFile::Open("../../../../source/X17geant4/X17geant4/mcs_3MeV.root", "READ");
  TTree* sim_tree  = static_cast<TTree*>(sim_input->Get("gasIonization"));

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

  TFile* output   = TFile::Open("../../../data/mcs_tracks/mcs_tracks2.root", "RECREATE");
  TTree* out_tree = new TTree("tracks_small", "Tree of tracks with map-drifted points");

  X17::TrackMicro track;
  out_tree->Branch("track_small", &track);
  std::vector<X17::MicroPoint> points;

  TParameter<UInt_t> pSeed("seed", seed);
  pSeed.Write();

  bool electron           = true;
  X17::Vector orientation = { 1, 0, 0 };
  X17::Vector origin      = { 6.51, 0, 0 };
  double kin_energy       = 8.0E+6 - X17::constants::E0;

  for (Long64_t i = 0; i < sim_tree->GetEntries(); i++)
  {
    ReportProgress(i, sim_tree->GetEntries());

    sim_tree->GetEntry(i);

    if (!driftSecondary && parentID != 1)
      continue;

    start.point /= 10;
    X17::EndPoint end = map->GetField(start.point).GetRandomPoint(rand.get());
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

  out_tree->Write();

  output->Close();
  sim_input->Close();
  map_input->Close();

  return 0;
}