// C++ dependencies
#include <iostream>
#include <vector>

// ROOT dependencies
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TTree.h"

// X17 dependencies
#include "../../include/Track.h"

bool compareTimes(X17::MicroPoint p1, X17::MicroPoint p2) 
{
    return p1.t0() < p2.t0();
}

void TrackStyle(TGraph* g)
{
    g->SetLineColor(kBlue);
    g->SetLineWidth(2);
}

void plot_drift()
{
    TFile* input = new TFile("grid_01/tracks_full1000.root");
    TTree* tracks = (TTree*)input->Get("tracks_full");

    X17::TrackMicro* track = nullptr;
    tracks->SetBranchAddress("track_full",&track);

    tracks->GetEntry(4);
    std::sort(track->points.begin(),track->points.end(),compareTimes);

    TMultiGraph* mg_lines_xy = new TMultiGraph("mg_lines_xy","Driftlines");
    TMultiGraph* mg_lines_xz = new TMultiGraph("mg_lines_xz","Driftlines");
    TMultiGraph* mg_lines_yz = new TMultiGraph("mg_lines_yz","Driftlines");

    TGraph* g_track_xy = new TGraph();
    TrackStyle(g_track_xy);
    TGraph* g_track_xz = new TGraph();
    TrackStyle(g_track_xz);
    TGraph* g_track_yz = new TGraph();
    TrackStyle(g_track_yz);

    int i = 1;

    for(X17::MicroPoint p : track->points)
    {
        g_track_xy->AddPoint(p.x0(),p.y0());
        g_track_xz->AddPoint(p.x0(),p.z0());
        g_track_yz->AddPoint(p.y0(),p.z0());
    }

    for(auto v : track->driftlines) 
    {
        std::cout << "Driftline " << i << "/" << track->driftlines.size() << " with " << v.size() << " points." << std::endl;

        TGraph* g_xy = new TGraph();
        g_xy->SetLineColor(kOrange);
        TGraph* g_xz = new TGraph();
        g_xz->SetLineColor(kOrange);
        TGraph* g_yz = new TGraph();
        g_yz->SetLineColor(kOrange);

        for(X17::DriftLinePoint p : v)
        {
            g_xy->AddPoint(p.x(),p.y());
            g_xz->AddPoint(p.x(),p.z());
            g_yz->AddPoint(p.y(),p.z());
        }

        mg_lines_xy->Add(g_xy);
        mg_lines_xz->Add(g_xz);
        mg_lines_yz->Add(g_yz);

        i++;
    }

    gStyle->SetTitleSize(0.05,"XY");
    gStyle->SetLabelSize(0.05,"XY");

    TCanvas* c_xy = new TCanvas("c_xy", "Driftlines XY",640,500);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.13);
    gPad->SetRightMargin(0.05);
    mg_lines_xy->GetXaxis()->SetTitle("x [cm]");
    mg_lines_xy->GetYaxis()->SetTitle("y [cm]");
    mg_lines_xy->Draw("AL");
    g_track_xy->Draw("L");

    TLegend* l_xy = new TLegend(0.8,0.8,1,1);
    l_xy->AddEntry(g_track_xy,"track");
    l_xy->AddEntry(mg_lines_xy,"driftlines");
    // l_xy->Draw();


    TCanvas* c_xz = new TCanvas("c_xz", "Driftlines XZ",640,500);
    mg_lines_xz->GetXaxis()->SetTitle("x [cm]");
    mg_lines_xz->GetYaxis()->SetTitle("z [cm]");
    mg_lines_xz->Draw("AL");
    g_track_xz->Draw("L");

    TLegend* l_xz = new TLegend(0.8,0.8,1,1);
    l_xz->AddEntry(g_track_xz,"track");
    l_xz->AddEntry(mg_lines_xz,"driftlines");
    // l_xz->Draw();


    TCanvas* c_yz = new TCanvas("c_yz", "Driftlines YZ",640,500);
    mg_lines_yz->GetXaxis()->SetTitle("y [cm]");
    mg_lines_yz->GetYaxis()->SetTitle("z [cm]");
    mg_lines_yz->Draw("AL");
    g_track_yz->Draw("L");

    TLegend* l_yz = new TLegend(0.8,0.8,1,1);
    l_yz->AddEntry(g_track_yz,"track");
    l_yz->AddEntry(mg_lines_yz,"driftlines");
    // l_yz->Draw();
}