// C++ dependencies
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

// ROOT dependencies
#include "Logger.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"

// Garfield++ dependencies
#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentGrid.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewSignal.hh"

// X17 dependencies
#include "AppManager.h"
#include "JobManager.h"
#include "Points.h"
#include "Track.h"
#include "Utilities.h"
#include "Vector.h"
#include "X17Utilities.h"

using namespace Garfield;
using namespace X17::constants;

// Input parameters: max_id, id, iterations, angle_bins, energy_bins (for grid-like simulation).
// For random simulation of one track run with no parameters.
int main(int argc, char* argv[])
{
  X17::AppManager man("micro_tracks", 3, "Simulation of microscopic tracks.");
  X17::Logger& logger      = X17::Logger::Get();
  X17::JobManager& job_man = X17::JobManager::Get(); // Job manager for parallel jobs.

  bool random_params;

  int iterations = 1;

  if (argc == 1)
  {
    random_params = true;
    logger.Essential("Running random generation of a single track, no extra parameters needed.");
  }
  else if (argc == 2)
  {
    random_params = true;
    iterations    = std::stoi(argv[1]);
    logger.Essential("Running random generation of " + std::to_string(iterations) + " tracks.");
  }
  else
  {
    random_params = false;

    // Check the number of paramaters passed to main function.
    if (argc < 4)
      logger.Fatal("Missing arguments in ion_electrons. Correct arguments: max_id, id, iterations.");

    // Set parameters of this job.
    int max_id = std::stoi(argv[1]);
    int id     = std::stoi(argv[2]);
    iterations = std::stoi(argv[3]);

    // Set parameters.
    int energy_bins = 11;
    job_man.RegisterParameterSteps("Ekin", 3e+6, 13e+6, energy_bins);
    int angle_bins   = 21;
    double theta_max = std::atan((X17::constants::win_height / 2) / X17::constants::xmin);
    job_man.RegisterParameterSteps("theta", -theta_max, theta_max, angle_bins);
    double phi_max = std::atan((X17::constants::win_width / 2) / X17::constants::xmin);
    job_man.RegisterParameterSteps("phi", -phi_max, phi_max, angle_bins);
    job_man.RegisterParameterSteps("charge", -1, 1, 2);

    if (id > max_id)
      logger.Fatal("Parameter id cannot be bigger than max_id.\n");

    job_man.Initialize(max_id, id, iterations);
  }

  // Set the output files.
  std::string folder_path = std::filesystem::current_path().string() + "/build/simulations/micro_tracks";
  std::string outPath     = GetNextFilePath(folder_path, "tracks_full");  // Big file containing the drift lines.
  std::string outPath2    = GetNextFilePath(folder_path, "tracks_small"); // Small file without the drift lines.
  logger.Essential("Full output will be saved to: " + outPath);
  logger.Essential("Small output will be saved to: " + outPath2);

  TTree* tracks
    = man.CreateOutputTree(outPath, "tracks_full", "Tree of simulated microscopic tracks (with drift lines)");

  X17::TrackMicro microtrack;
  tracks->Branch("track_full", &microtrack);

  // Set the gas mixture.
  MediumMagboltz gas;
  gas.SetComposition("ar", 70., "co2", 30.);

  // Add magnetic and electric field.
  ComponentGrid grid;
  const double m2cm = 100.;
  grid.LoadMagneticField("data/elmag/VecB2.txt", "xyz", m2cm);
  grid.LoadElectricField("data/elmag/VecE2.txt", "xyz", false, false, m2cm);
  grid.SetMedium(&gas);

  // Assemble a sensor.
  Sensor sensor;
  sensor.AddComponent(&grid);
  constexpr double space = 1; // Extra space on the sensor so electrons always end on zmax.
  sensor.SetArea(xmin - space, -yhigh - space, zmin - space, xmax + space, yhigh + space, zmax);

  // We use microscopic tracking for the electron avalanche simulation.
  AvalancheMicroscopic aval;
  aval.SetSensor(&sensor);
  aval.EnableSignalCalculation();
  aval.EnableMagneticField();
  aval.EnableDriftLines();

  // Random number generator for the simulation.
  std::unique_ptr<TRandom3> rand = man.CreateRNG();

  logger.Info("Starting simulation.").PushIndent();

  int min_set = random_params ? 1 : job_man.GetMinIndex();
  int max_set = random_params ? 1 : job_man.GetMaxIndex();

  // Loop over the tracks. Generate each track as many times as is the given number of iterations.
  for (int i = min_set; i <= max_set; i++)
    for (int j = 0; j < iterations; j++)
    {
      std::vector<X17::MicroPoint> points;
      std::vector<std::vector<X17::DriftLinePoint>> driftlines;

      // Track parameters.
      bool electron;
      X17::Vector origin, orientation;
      double kin_en, theta, phi;

      // Get initial track parameters (random or grid-like).
      if (random_params)
        X17::GetRandomTrackParams(rand.get(), electron, origin, orientation, kin_en);
      else
      {
        electron = job_man.GetParameterValue("charge", i) < 0;
        theta    = job_man.GetParameterValue("theta", i);
        phi      = job_man.GetParameterValue("phi", i);
        kin_en   = job_man.GetParameterValue("Ekin", i);

        orientation = { std::cos(phi) * std::cos(theta), std::sin(phi) * std::cos(theta), std::sin(theta) };

        // Setting the origin point from the orientation vector (assuming straight line motion from (0,0,0)).
        origin = orientation * X17::constants::xmin / (std::cos(phi) * std::cos(theta));
      }

      logger.Print("TRACK No." + std::to_string(i) + ", iteration " + std::to_string(j + 1) + ":").PushIndent();
      std::string particle = electron ? "electron" : "positron";
      logger.Print("particle:    " + particle);
      logger.Print("origin:      " + origin.ToString());
      logger.Print("orientation: " + orientation.ToString());
      job_man.PrintParameterValues(i).PopIndent();

      // Simulate an ionizing particle using Heed.
      TrackHeed track;
      track.SetParticle(particle);

      track.SetKineticEnergy(kin_en); // Set the particle kinetic energy [eV].
      track.SetSensor(&sensor);
      track.EnableMagneticField();
      track.EnableElectricField();
      track.DisableDeltaElectronTransport(); // This will disable secondary electrons in the track.
      track.EnablePhotonReabsorption(false); // Enable/disable fluorescence reabsorption.

      // Get the default parameters.
      double maxrange = 0., rforstraight = 0., stepstraight = 0., stepcurved = 0.;
      track.GetSteppingLimits(maxrange, rforstraight, stepstraight, stepcurved);

      // Reduce the step size [rad].
      stepcurved = 0.04;
      maxrange   = 0.2;
      track.SetSteppingLimits(maxrange, rforstraight, stepstraight, stepcurved);

      // Set the starting point and momentum vector of the particle.
      double xt = origin.x; // [cm]
      double yt = origin.y; // [cm]
      double zt = origin.z; // [cm]
      double ti = 0;        // [ns]
      double px = orientation.x;
      double py = orientation.y;
      double pz = orientation.z;

      // Simulate the track.
      track.NewTrack(xt, yt, zt, ti, px, py, pz);

      // Loop over the clusters.
      double xc, yc, zc, tc, ec, extra;
      int nc;
      int n_electron = 0;
      while (track.GetCluster(xc, yc, zc, tc, nc, ec, extra))
      {
        for (int j = 0; j < nc; ++j)
        {
          X17::MicroPoint point;
          double xe, ye, ze, te, ee, dxe, dye, dze;
          track.GetElectron(j, xe, ye, ze, te, ee, dxe, dye, dze);

          n_electron++;
          std::cout << "Distance to origin: " << std::sqrt(xe * xe + ye * ye + ze * ze) << "  time " << te
                    << "  number " << n_electron << '\n';

          // Simulate the drift/avalanche of this electron.
          aval.AvalancheElectron(xe, ye, ze, te, ee, dxe, dye, dze);

          // Move electrons that hit the mesh plane into the amplification gap.
          int status;
          aval.GetElectronEndpoint(0, point.start.point.x, point.start.point.y, point.start.point.z, point.start.t,
                                   point.e0, point.end.point.x, point.end.point.y, point.end.point.z, point.end.t,
                                   point.e1, status);
          points.push_back(point);

          // Save driftlines. (Every 10th point)
          std::vector<X17::DriftLinePoint> driftline;

          for (int k = 0; k < aval.GetNumberOfElectronDriftLinePoints(); k += 10)
          {
            X17::DriftLinePoint dl_point;
            aval.GetElectronDriftLinePoint(dl_point.point.x, dl_point.point.y, dl_point.point.z, dl_point.t, k);
            driftline.push_back(dl_point);
          }

          driftlines.push_back(driftline);
        }

        // break; // Only for fast testing (simulates only one electron)!!!
      }

      microtrack = X17::TrackMicro(electron, points, origin, orientation, kin_en, driftlines);

      tracks->Fill();
    }

  logger.PopIndent();

  // 1. Deactivate the driftlines branch.
  tracks->SetBranchStatus("*driftlines*", false);

  // 2. Open the second file manually to avoid triggering AppManager's auto-cleanup.
  TFile outFile2(outPath2.c_str(), "RECREATE");
  outFile2.cd(); // explicitly set gDirectory to the new file

  // 3. Clone the tree. CloneTree() automatically attaches the new tree to gDirectory (outFile2).
  TTree* tracks2 = tracks->CloneTree();
  tracks2->SetName("tracks_small");
  tracks2->SetTitle("Tree of simulated microscopic tracks (no drift lines)");

  // 4. Write and close the small file.
  tracks2->Write();
  outFile2.Close();

  // 5. Re-enable the branch on the original tree so AppManager::~AppManager()
  // writes the final buffer of the full tree correctly when the program exits.
  tracks->SetBranchStatus("*driftlines*", true);

  man.GetOutputFile()->cd();

  return 0;
}