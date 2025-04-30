#include <iostream>
#include <string>
#include <vector>

#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TImage.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "../../build/X17_dict.cxx"

#include "CircleFit3D.h"
#include "Field.h"
#include "RK4.h"
#include "Track.h"
#include "X17Utilities.h"

X17::TrackMicro ParseOldTrack(TTree* tree, bool newfile)
{
    X17::TrackMicro track;
    X17::MicroPoint* point = nullptr;

    if (newfile)
    {
        tree->SetBranchAddress("point",&point);
    }
    
    else
    {
        point = new X17::MicroPoint();
        tree->SetBranchAddress("x0",&point->start.point.x);
        tree->SetBranchAddress("y0",&point->start.point.y);
        tree->SetBranchAddress("z0",&point->start.point.z);
        tree->SetBranchAddress("t0",&point->start.t);
        tree->SetBranchAddress("e0",&point->e0);
        tree->SetBranchAddress("x1",&point->end.point.x);
        tree->SetBranchAddress("y1",&point->end.point.y);
        tree->SetBranchAddress("z1",&point->end.point.z);
        tree->SetBranchAddress("t1",&point->end.t);
        tree->SetBranchAddress("e1",&point->e1);
    }

    for (int i = 0; i < tree->GetEntries(); i++)
    {
        tree->GetEntry(i);
        track.points.push_back(*point);
    }

    return track;    
}

/// @brief Load track (electron start and endpoints) from a single track file
/// @param folder Name of the folder in the data/single_track directory.
/// @param newcoords True = x,y,z; False = z,x,y
/// @param newfile True = doubles; False = X17::MicroPoint
/// @return The track as a X17::TrackMicro object
X17::TrackMicro LoadTrack(std::string folder, bool newcoords = false, bool newfile = false)
{
    std::string track_path = "../../../data/single_track/" + folder + "/electrons.root";
    TFile* track_file = new TFile(track_path.c_str());
    TTree* track_tree = (TTree*)track_file->Get("electrons");
    return ParseOldTrack(track_tree,newfile);
}

double g_cwidth = 700;
double g_cheight = 500;
double g_xmin = 0;
double g_xmax = 15;
double g_ymin = -0.3;
double g_ymax = 0.3;
double g_zmin = -3;
double g_zmax = 8;

void PlotDriftXZ(X17::TrackMicro track, bool newcoords = false)
{
    TCanvas* c1 = new TCanvas();
        TGraph* g1 = new TGraph();
        std::vector<TLine*> driftlines;
        for (const auto& point : track.points)
        {
            if (newcoords)
            {
                g1->SetPoint(g1->GetN(),point.start.x(),point.start.z());
                driftlines.emplace_back(new TLine(point.start.x(),point.start.z(),point.end.x(),point.end.z()));
            }

            else
            {
                g1->SetPoint(g1->GetN(),point.start.z(),point.start.y());
                driftlines.emplace_back(new TLine(point.start.z(),point.start.y(),point.end.z(),point.end.y()));
            }
        }

        g1->GetXaxis()->SetTitle("x [cm]");
        g1->GetYaxis()->SetTitle("z [cm]");
        g1->GetXaxis()->SetRangeUser(-0.1,15.1);
        g1->GetYaxis()->SetRangeUser(-2.7,8.1);
        g1->Draw("AL");

        for (const auto& line : driftlines)
        {
            line->SetLineColor(X17::Color::GOrange());
            line->SetLineWidth(2);
            line->Draw();
        }

        g1->SetLineWidth(3);
        g1->Draw("same");
}

void PlotDriftXY(X17::TrackMicro track, std::string filename, bool newcoords = false)
{
    TCanvas* c = new TCanvas("","Drift XY",g_cwidth,g_cheight);
    TImage* im_xy = TImage::Open(("../../../tests/thesis_plots/" + filename).c_str());
    
    TH2F* h2 = new TH2F("h2","",10,g_xmin,g_xmax,10,g_ymin,g_ymax);
    h2->SetStats(0);
    h2->GetXaxis()->SetTitle("x [cm]");
    h2->GetYaxis()->SetTitle("y [cm]");
    h2->Draw();

    double l = c->GetLeftMargin();
    double r = c->GetRightMargin();
    double b = c->GetBottomMargin();
    double t = c->GetTopMargin();
    
    // weird factor necessary to get correct size
    im_xy->Scale(1.00*g_cwidth*(1-l-r)*15/(g_xmax-g_xmin),(1414./1478.)*g_cheight*(1-b-t)*0.6/(g_ymax-g_ymin));
    im_xy->Draw();
    
    // c2->cd();
    h2->Draw("AXIS same");
    TLine* ltop   = new TLine(g_xmin,g_ymax,g_xmax,g_ymax);
    TLine* lright = new TLine(g_xmax,g_ymin,g_xmax,g_ymax);
    ltop->Draw();
    lright->Draw();

    TGraph* g2 = new TGraph();
    for (const auto& point : track.points)
    {
        if (newcoords)
        {
            g2->SetPoint(g2->GetN(),point.start.x(),point.start.y());
        }
        
        else
        {
            g2->SetPoint(g2->GetN(),point.start.z(),point.start.x());
        }
    }
    g2->SetLineWidth(3);
    g2->Draw("L same");
}

void PlotDriftYZ(X17::TrackMicro track, std::string filename, bool newcoords = false)
{
    TCanvas* c = new TCanvas("","Drift YZ",g_cwidth,g_cheight);
    TImage* im_yz = TImage::Open(("../../../tests/thesis_plots/" + filename).c_str());

    TH2F* h3 = new TH2F("h3","",10,g_ymin,g_ymax,10,g_zmin,g_zmax);
    h3->SetStats(0);
    h3->GetXaxis()->SetTitle("y [cm]");
    h3->GetYaxis()->SetTitle("z [cm]");
    h3->Draw();

    double l = c->GetLeftMargin();
    double r = c->GetRightMargin();
    double b = c->GetBottomMargin();
    double t = c->GetTopMargin();

    // weird factor necessary to get correct size
    im_yz->Scale(1.004*g_cwidth*(1-l-r)*0.6/(g_ymax-g_ymin),0.957*g_cheight*(1-b-t)*11/(g_zmax-g_zmin));
    im_yz->Draw();
    
    // c3->cd();
    h3->Draw("AXIS same");
    TLine* ltop2   = new TLine(g_ymin,g_zmax,g_ymax,g_zmax);
    TLine* lright2 = new TLine(g_ymax,g_zmin,g_ymax,g_zmax);
    ltop2->Draw();
    lright2->Draw();

    TGraph* g3 = new TGraph();
    for (const auto& point : track.points)
    {
        if (newcoords)
        {
            g3->SetPoint(g3->GetN(),point.start.y(),point.start.z());
        }
        
        else
        {
            g3->SetPoint(g3->GetN(),point.start.x(),point.start.y());
        }
    }
    g3->SetLineWidth(3);
    g3->Draw("L same");
}

int main(int argc, char *argv[])
{
    TApplication app("app", &argc, argv);
    std::cout << "Working directory: " << gSystem->WorkingDirectory() << std::endl;

    int orange_index = TColor::GetFreeColorIndex();
    TColor* orange = new TColor(orange_index,1,0.6,0);

    X17::Field<X17::Vector>* magfield = X17::LoadField("../../../data/elmag/VecB2.txt",{-20,-30,-30},{20,30,30},0.5);

    X17::TrackMicro track1 = LoadTrack("original");
    PlotDriftXZ(track1);
    PlotDriftXY(track1,"drift_xy_9010.png");
    PlotDriftYZ(track1,"drift_yz_9010.png");

    X17::TrackMicro track2 = LoadTrack("new_7030");
    PlotDriftXZ(track2,true);
    PlotDriftXY(track2,"drift_xy_7030.png",true);
    PlotDriftYZ(track2,"drift_yz_7030.png",true);

    // Reconstruct CircleFit3D + RK4
        // std::vector<X17::RecoPoint> reco_points;
        // for (const auto& point : original_track.points)
        // {
        //     if (newcoords) reco_points.emplace_back(point.start.x(),point.start.y(),point.start.z(),1);
        //     else reco_points.emplace_back(point.start.z(),point.start.x(),point.start.y(),1);
        // }

        // X17::CircleFit3D circlefit = X17::CircleFit3D(X17::Vector(0,0,0), X17::Vector(1,0,0));
        // for (const auto& point : reco_points) circlefit.AddPoint(point);
        // circlefit.SetFitter(4,false);
        // circlefit.SetParameters(0,10,1.5,true);
        // circlefit.FitCircle3D();
        // circlefit.PrintFitParams();

        // double cfit3d_E_mid = circlefit.GetEnergy(*magfield,true);
        // double cfit3d_E_avg = circlefit.GetEnergy(*magfield,false);

        // std::cout << "Reconstructed energy (circle fit no pads): " << cfit3d_E_mid << " (middle field), " << cfit3d_E_avg << " (average field).\n";

        // X17::RKFit fit = X17::RKFit(magfield, true, 1E-13, X17::Vector(0,0,0), X17::Vector(1,0,0), reco_points, true);
        // fit.SetFitter();
        // fit.SetEnergy(cfit3d_E_mid);
        // fit.FitRK();
        // fit.PrintFitParams();

    app.Run(true);

    return 0;
}
