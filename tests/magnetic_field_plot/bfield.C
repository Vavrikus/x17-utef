#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>

#include "Garfield/ComponentGrid.hh"
#include "Garfield/ViewField.hh"

using namespace Garfield;

int main(int argc, char *argv[]) {

  TApplication app("app", &argc, argv);

  // Load the field map.
  ComponentGrid grid;
  const double m2cm = 100.;
  grid.LoadMagneticField("../../mag_data/VecB.txt", "xyz", m2cm);
  grid.LoadElectricField("../../mag_data/VecE.txt", "xyz",false,false, m2cm);
  //~ grid.GetElectricField(false);

  ViewField view;
  view.SetComponent(&grid);
  //~ view.SetPlaneYZ();
   view.SetPlaneXZ();
  //view.SetPlane(0,1,0,0,-20,0);

  // Get the mesh parameters.
  unsigned int nx = 0, ny = 0, nz = 0;
  double xMin = 0., yMin = 0., zMin = 0., xMax = 0., yMax = 0., zMax = 0.;
  grid.GetMesh(nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax);
  view.SetArea(xMin, zMin, xMax, zMax);
  TCanvas c1("c1", "", 902, 600);
  view.SetCanvas(&c1);
  view.PlotContour("bmag");
  
  // TCanvas c2("c2", "", 600, 600);
  // view.SetCanvas(&c2);
  // view.PlotProfile(0., 0., zMin, 0., 0., zMax, "by");

  app.Run(true);
}
