#include <iostream>
#include <cmath>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <TGraph2D.h>
#include <TTree.h>

#include "Garfield/ViewField.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/ComponentGrid.hh"

using namespace Garfield;
using namespace std;

int main(int argc, char * argv[]){

  TApplication app("app", &argc, argv);

  // Set the gas mixture.
  MediumMagboltz gas;
  gas.SetComposition("ar", 90., "co2", 10.);

  // ofstream myfile;
  // myfile.open ("electrons.txt"); //name of the file to save data
  TFile outFile("electrons.root","RECREATE","Electrons from ionization track");
  TTree electrons("electrons","Tree of initial and final points of electrons");

  double x0, y0, z0, t0, e0;
  double x1, y1, z1, t1, e1;
  electrons.Branch("x0",&x0);
  electrons.Branch("y0",&y0);
  electrons.Branch("z0",&z0);
  electrons.Branch("t0",&t0);
  electrons.Branch("e0",&e0);
  electrons.Branch("x1",&x1);
  electrons.Branch("y1",&y1);
  electrons.Branch("z1",&z1);
  electrons.Branch("t1",&t1);
  electrons.Branch("e1",&e1);

  
  ComponentGrid grid;
  const double m2cm = 100.;
  grid.LoadMagneticField("VecB.txt", "xyz", m2cm); 
  grid.LoadElectricField("VecE.txt", "xyz",false,false, m2cm);
  grid.SetMedium(&gas);

  
  // 8 cm drift gap.
  constexpr double yDrift = 8;


  // Assemble a sensor.
  Sensor sensor;
  sensor.AddComponent(&grid); 
  sensor.SetArea(-1, -yDrift, -15, 1, yDrift, 15.);


  //~ // We use microscopic tracking for simulating the electron avalanche.
  AvalancheMicroscopic aval;
  aval.SetSensor(&sensor); 
  // Switch on signal calculation. 
  aval.EnableSignalCalculation(); 
  aval.EnableMagneticField();
  
  //~ // Simulate an ionizing particle (negative pion) using Heed.
  TrackHeed track;
  track.SetParticle("electron"); //Set particle type
  constexpr double momentum = 8.e6; // [eV / c] //Set particle momentum
  track.SetMomentum(momentum);
  track.SetSensor(&sensor);
  track.EnableMagneticField();
  track.EnableElectricField();
  track.DisableDeltaElectronTransport();  //This will disable secondary electrons in the track
  track.EnablePhotonReabsorption(false);  //enable/disable fluorescence reabsorption

  //~ // Construct a viewer to visualise the drift lines.
  ViewDrift driftView;
  track.EnablePlotting(&driftView);
  aval.EnablePlotting(&driftView);
  aval.EnableExcitationMarkers(false);
  aval.EnableIonisationMarkers(false);
  //~ // Get the default parameters.
  double maxrange = 0., rforstraight = 0., stepstraight = 0., stepcurved = 0.;
  track.GetSteppingLimits(maxrange, rforstraight, stepstraight, stepcurved);
  // Reduce the step size [rad].
  
  stepcurved = 0.04;
  maxrange = 0.2;
  track.SetSteppingLimits(maxrange, rforstraight, stepstraight, stepcurved);

  TH1D *h1 = new TH1D("htime","htime",100,0,4000);
  TGraph2D *gr1 = new TGraph2D();

  // Set the starting point and momentum vector of the particle. Randomize in the future
  double xt = 0.0; // [cm]
  double yt = 0.0;
  double zt = 0.0;
  double ti = 0;
  double px = 0;
  double py = 0;
  double pz = 1;
  int k=0;

  // Now simulate a track, with p0 = ẑ

  track.NewTrack(xt, yt, zt, ti, px, py, pz);
  // Loop over the clusters.
  double xc, yc, zc, tc, ec, extra; //parameters
  int nc;
  while (track.GetCluster(xc, yc, zc, tc, nc, ec, extra)) 
  {
    cout << "distance to origin " << sqrt(xc*xc+yc*yc+zc*zc) << "  time " << tc << "  number " << k << "\n";
    for (int j = 0; j < nc; ++j) 
    {
      double xe, ye, ze, te, ee, dxe, dye, dze;
      track.GetElectron(j, xe, ye, ze, te, ee, dxe, dye, dze);
      // Simulate the drift/avalanche of this electron.
      aval.AvalancheElectron(xe, ye, ze, te, 0.1, dxe, dye, dze);
      // Move electrons that hit the mesh plane into the amplification gap.
      int status;
      aval.GetElectronEndpoint(0, x0, y0, z0, t0, e0, 
                                  x1, y1, z1, t1, e1, status);
      electrons.Fill();
      k++;
    }
    // if(k == 50) break;
  }

  std::cout << "total electron generated: "<< k <<std::endl;
  
  // Plotting
  // These canvas will be used to display the drift lines and the field
  if(bool make_plots=true){
  gStyle->SetPadRightMargin(0.15);
  TCanvas* c1 = new TCanvas("c1", "", 800, 800);
  driftView.SetPlane(1,0, 0, 0, 0, 0);
  driftView.SetArea(-15, -8, 15, 8);
  driftView.SetCanvas(c1);
  constexpr bool twod = true;
  driftView.Plot(twod);
  c1->SaveAs("plot.pdf");
  
  TCanvas* c2 = new TCanvas("c2", "", 800, 800);
  driftView.SetPlane(0,1,0, 0, 0, 0);
  driftView.SetArea(-1, -15, 1, 15);
  driftView.SetCanvas(c2);
  constexpr bool twod2 = true;
  driftView.Plot(twod2);
  c2->SaveAs("plot2.pdf");  
  
  TCanvas* c3 = new TCanvas("c3", "", 800, 800);
  driftView.SetPlane(0,0,1, 0, 0, 0);
  driftView.SetArea(-1, -8, 1, 8);
  driftView.SetCanvas(c3);
  constexpr bool twod3 = true;
  driftView.Plot(twod3);
  c3->SaveAs("plot3.pdf");  
  
  TCanvas* c4 = new TCanvas("c4", "", 800, 800);
  h1->Draw();
  c4->SaveAs("plottime.pdf");  
  
  TCanvas* c5 = new TCanvas("c5", "", 800, 800);
  gr1->GetXaxis()->SetTitle("X-axis(cm)");
  gr1->GetYaxis()->SetTitle("Y-axis(cm)");
  gr1->GetZaxis()->SetTitle("Z-axis(cm)");
  gr1->Draw("P");
  c5->SaveAs("trajectory.pdf");  
  }

  electrons.Print();
  outFile.Write();
  outFile.Close();

  // TCanvas* c6 = new TCanvas("c6", "", 800, 800);
  // electrons.Draw("z0:t1","","ap");
  app.Run(true);
  return 0;

}
