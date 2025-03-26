#include <iostream>

#include <TApplication.h>
#include <TArrow.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TLine.h>
#include <TPaletteAxis.h>
#include <TROOT.h>
#include <TStyle.h>

#include "Garfield/ComponentGrid.hh"
#include "Garfield/ViewField.hh"

#include "Field.h"
#include "X17Utilities.h"

using namespace Garfield;

int main(int argc, char *argv[])
{
	TApplication app("app", &argc, argv);
	
	// Load the field map.
	// ComponentGrid grid;
	// const double m2cm = 100.;
	// grid.LoadMagneticField("../../../data/elmag/VecB2.txt", "xyz", m2cm);
	// grid.LoadElectricField("../../../data/elmag/VecE2.txt", "xyz", false, false, m2cm);
	// // grid.GetElectricField(false);
	
	// ViewField view;
	// view.SetComponent(&grid);
	// view.SetPlaneXY();
	// // view.SetPlaneXZ();
	// // view.SetPlane(0,1,0,0,-20,0);
	
	// // Get the mesh parameters.
	// unsigned int nx = 0, ny = 0, nz = 1;
	// double xMin = -0.2, yMin = -0.2, zMin = -0.2, xMax = 0.2, yMax = 0.2, zMax = 0.2;
	// grid.GetMesh(nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax);
	// view.SetArea(xMin, yMin, xMax, yMax);
	
	// TCanvas c1("c1", "", 902, 600);
	// view.SetCanvas(&c1);
	// view.PlotContour("bmag");
	
	// TCanvas c2("c2", "", 600, 600);
	// view.SetCanvas(&c2);
	// view.PlotProfile(-20, 0, 0, 20, 0, 0, "bmag", false);

	// TCanvas c3("c3", "", 600, 600);
	// view.SetCanvas(&c3);
	// view.Plot("bmag");

	// Loading the magnetic field data.
	X17::Field<X17::Vector>* magfield = X17::LoadField("../../../data/elmag/VecB2.txt", {-20,-30,-30}, {20,30,30}, 0.5);
	
	// Grid size
	double xmax = X17::constants::xmax + 1.5;
	double ymax = X17::constants::yhigh + 1.55;
	double zmax = X17::constants::zmax + 1.75;

	double spacing = 0.75;			// Grid spacing
	double scale   = 6.0 * spacing; // Scaling factor for arrows
	double padding = 1.0;
	
	TCanvas *c = new TCanvas("c", "Vector Field", 800, 800);
	c->SetLeftMargin(0.12);
	c->SetRightMargin(0.15);

	TH2F *frame = new TH2F("frame", ";x [cm];y [cm];B [T]", 10, -padding, xmax+padding, 10, -ymax-padding, ymax+padding);
	frame->SetStats(0);
	frame->SetBinContent(1,0);
	frame->SetMinimum(0);
	frame->SetMaximum(255.0/750.0);
	frame->Draw("colz");

	X17::DrawTrapezoid(false);

	for (double x = 0; x <= xmax; x += spacing) {
		for (double y = -ymax; y <= ymax; y += spacing) {
			X17::Vector b = magfield->GetField(x,y,0);
			
			// Draw arrow (x,y) → (x+vx*scale, y+vy*scale)
			if(scale*b.Magnitude() < 2*spacing)
			{
				int colorIndex = gStyle->GetColorPalette((int)(7.5 * b.Magnitude() * 100) % 255);

				TArrow *arrow = new TArrow(x, y, x + scale * b.x, y + scale * b.y, b.Magnitude()/12, "|>");
				arrow->SetLineColor(colorIndex);
				arrow->SetFillColor(colorIndex);
				arrow->Draw();
			}
		}
	}

	c->Draw();

	TCanvas *c2 = new TCanvas("c2", "XZ plane magnitude", 800, 800);
	c2->SetLeftMargin(0.12);
	c2->SetRightMargin(0.16);

	TH2F* h_field = new TH2F("h_field", ";x [cm];z [cm];B [T]", 1000, 0, xmax, 1000, -zmax, zmax);
	h_field->SetStats(0);

	for(int i = 0; i < 1000; i++)
	{
		for(int j = 0; j < 1000; j++)
		{
			double x = h_field->GetXaxis()->GetBinCenter(i);
			double z = h_field->GetYaxis()->GetBinCenter(j);
			X17::Vector b = magfield->GetField(x,0,z);
			h_field->SetBinContent(i,j,b.Magnitude());
		}
	}

	gStyle->SetNumberContours(255);
	h_field->Draw("colz");

	TLine *line = new TLine(X17::constants::xmax, X17::constants::zmax, X17::constants::xmax, -X17::constants::zmax);
	line->SetLineWidth(2);
	line->Draw();
	line = new TLine(X17::constants::xmax, -X17::constants::zmax, X17::constants::xmin, -X17::constants::zmax);
	line->SetLineWidth(2);
	line->Draw();
	line = new TLine(X17::constants::xmin, -X17::constants::zmax, X17::constants::xmin, X17::constants::zmax);
	line->SetLineWidth(2);
	line->Draw();
	line = new TLine(X17::constants::xmin, X17::constants::zmax, X17::constants::xmax, X17::constants::zmax);
	line->SetLineWidth(2);
	line->Draw();

	h_field->GetZaxis()->SetTitleOffset(1.5);

	// TCanvas *c3 = new TCanvas("c3", "Field along z=0.45 in xz plane", 800, 800);
	// c3->SetLeftMargin(0.12);
	// c3->SetRightMargin(0.05);

	// TH1F* h_field_xz = new TH1F("h_field_xz", ";x [cm];B [T]", 1000, 0, xmax);
	// h_field_xz->SetStats(0);

	// for(int i = 0; i < 1000; i++)
	// {
	// 	double x = h_field_xz->GetBinCenter(i);
	// 	X17::Vector b = magfield->GetField(x,0,0.45);
	// 	std::cout << "x = " << x << " b.y = " << b.y << std::endl;
	// 	h_field_xz->SetBinContent(i,b.y);
	// }

	// h_field_xz->Draw();

	app.Run(true);
}
