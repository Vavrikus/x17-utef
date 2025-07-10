// C++ dependencies
#include <iostream>
#include <string>
#include <vector>

// ROOT dependencies
#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TImage.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TPolyLine3D.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "../../build/X17_dict.cxx"

// X17 dependencies
#include "CircleFit2D.h"
#include "CircleFit3D.h"
#include "Color.h"
#include "Field.h"
#include "NSpline.h"
#include "PadLayout.h"
#include "Reconstruction.h"
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

/// @brief Load track (electron start and endpoints) from a single track file.
/// @param folder Name of the folder in the data/single_track directory.
/// @param newcoords True = x,y,z; False = z,x,y
/// @param newfile True = doubles; False = X17::MicroPoint
/// @return The track as a X17::TrackMicro object.
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

TGraph* g_g_res1 = nullptr;
TGraph* g_g_res2 = nullptr;

void PlotDriftXZ(X17::TrackMicro track, bool newcoords = false)
{
    TCanvas* c1 = new TCanvas();
        TGraph* g1 = new TGraph();
        std::vector<TLine*> driftlines;
        for (const auto& point : track.points)
        {
            if (newcoords)
            {
                g1->AddPoint(point.start.x(),point.start.z());
                driftlines.emplace_back(new TLine(point.start.x(),point.start.z(),point.end.x(),point.end.z()));
            }

            else
            {
                g1->AddPoint(point.start.z(),point.start.y());
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
            g2->AddPoint(point.start.x(),point.start.y());
        }
        
        else
        {
            g2->AddPoint(point.start.z(),point.start.x());
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
            g3->AddPoint(point.start.y(),point.start.z());
        }
        
        else
        {
            g3->AddPoint(point.start.x(),point.start.y());
        }
    }
    g3->SetLineWidth(3);
    g3->Draw("L same");
}

void PlotTrackRK(X17::TrackMicro track, X17::Field<X17::Vector>* magfield, bool newcoords = false)
{
    double Ekin = std::sqrt(64E+12 + X17::constants::E0*X17::constants::E0) - X17::constants::E0;
        X17::RK4<8>* trackrk = X17::GetTrackRK(*magfield,true,1E-13,Ekin,X17::Vector(0,0,0),X17::Vector(1,0,0),true);
        trackrk->Integrate();
        std::vector<X17::RKPoint> rk_points;
        for (const auto& vec : trackrk->GetResults())
        {
            using namespace X17::constants;
            rk_points.push_back(X17::RKPoint(m2cm * vec.at(1,0), m2cm * vec.at(2,0), m2cm * vec.at(3,0), 1e+9 / c * vec.at(0,0)));
        }

        TGraph* g_micro = new TGraph();
        TGraph* g_res = new TGraph();
        TGraph* g_diff = new TGraph();
        TGraph* g_diff2 = new TGraph();
        for (const auto& point : track.points)
        {
            X17::Vector start_old = point.start.point;
            X17::Vector start_new = newcoords ? X17::Vector{start_old.x,start_old.y,start_old.z} : X17::Vector{start_old.z,start_old.x,start_old.y};

            X17::Vector closest_point;
            double res2 = X17::GetTrackRKSqDistAndCP(trackrk,start_new,closest_point);
            g_res->AddPoint(start_new.x,10000*std::sqrt(res2));

            X17::Vector exaggerated = closest_point + 1000 * (start_new-closest_point);
            g_micro->AddPoint(exaggerated.x,exaggerated.z);

            X17::Vector diff = start_new - closest_point;
            g_diff->AddPoint(start_new.x,10000*diff.z);
            g_diff2->AddPoint(start_new.x,10000*diff.x);
        }

        TCanvas* c_rk = new TCanvas("","Drift XZ",g_cwidth,g_cheight);
        g_micro->SetMarkerStyle(2);
        g_micro->SetMarkerColor(kRed);
        g_micro->SetLineColor(kNone);
        g_micro->GetXaxis()->SetTitle("x [cm]");
        g_micro->GetYaxis()->SetTitle("z [cm]");
        g_micro->Draw("AP");

        TGraph* g_rk = new TGraph();
        for (const auto& point : rk_points)
            g_rk->AddPoint(point.x(),point.z());
        g_rk->SetLineColor(kBlue);
        g_rk->SetLineWidth(2);
        g_rk->Draw("L same");

        TLegend* l_rk = new TLegend(0.58,0.81,0.93,0.93);
        l_rk->AddEntry(g_rk,"Runge-Kutta track");
        l_rk->AddEntry(g_micro,"HEED ion. electrons");
        l_rk->SetTextSize(0.043);
        l_rk->Draw();

        TCanvas* c_res = new TCanvas("","Residuals",g_cwidth,g_cheight);
        g_res->SetMarkerStyle(2);
        g_res->SetMarkerColor(kRed);
        g_res->SetLineColor(kNone);
        g_res->GetXaxis()->SetTitle("x [cm]");
        g_res->GetYaxis()->SetTitle("residuals [#mum]");
        g_res->Draw("AP");

        TLegend* l_res = new TLegend(0.54,0.82,0.93,0.93);
        l_res->AddEntry(g_res,"#splitline{HEED electron}{residuals to RK4 track}","p");
        l_res->SetTextSize(0.043);
        l_res->Draw();

        //TCanvas* c_diff = new TCanvas("c_diff","Diff",g_cwidth,g_cheight);
        // g_diff->SetMarkerStyle(2);
        // g_diff->SetMarkerColor(kRed);
        // g_diff->Draw("AP");
        // g_diff2->SetMarkerStyle(2);
        // g_diff2->SetMarkerColor(kBlue);
        // g_diff2->Draw("P same");
}

void PlotRASD(X17::TrackMicro track, X17::Field<X17::MapPoint>* map, bool newcoords = false)
{
    TGraph* g_drift_noTPC = new TGraph();
    TGraph* g_drift_TPC = new TGraph();
    for (const auto& point : track.points)
    {
        using namespace X17::constants;
        double z_end = newcoords ? point.end.z() : point.end.y();
        if(z_end > 7.5)
        {
            double x = newcoords ? point.start.x() : point.start.z();
            double z = newcoords ? point.start.z() : point.start.y();
            double t = point.end.t/1000;

            if(x < xmin)
                g_drift_noTPC->AddPoint(8-z,t);
            else
                g_drift_TPC->AddPoint(8-z,t);
        }
    }
    g_drift_noTPC->SetMarkerStyle(2);
    g_drift_noTPC->SetMarkerColor(kBlue);
    g_drift_TPC->SetMarkerStyle(2);
    g_drift_TPC->SetMarkerColor(kRed);

    TF1* f_drift_TPC = new TF1("f_drift_TPC","pol1",7.5,15);
    g_drift_TPC->Fit(f_drift_TPC,"0");
    f_drift_TPC->SetLineColor(kBlack);
    f_drift_TPC->SetLineWidth(2);
    double a0 = f_drift_TPC->GetParameter(0);
    double a1 = f_drift_TPC->GetParameter(1);
    double b0 = -a0/a1; //inverse polynomial param
    double b1 = 1.0/a1; //inverse polynomial param

    std::cout << "Reconstructed v_d: " << b1 << " cm/us, d_0: " << b0 << " cm\n";

    TCanvas* c_drift = new TCanvas("","Drift",g_cwidth,g_cheight);

    TMultiGraph* mg_drift = new TMultiGraph();
    mg_drift->Add(g_drift_noTPC,"P");
    mg_drift->Add(g_drift_TPC,"P");
    mg_drift->GetXaxis()->SetTitle("distance to readout [cm]");
    mg_drift->GetYaxis()->SetTitle("t [#mus]");
    mg_drift->Draw("A");
    f_drift_TPC->Draw("same");

    TLegend* l_drift = new TLegend(0.17,0.745,0.5,0.91);
    l_drift->AddEntry(g_drift_noTPC,"outside TPC","p");
    l_drift->AddEntry(g_drift_TPC,"inside TPC","p");
    l_drift->AddEntry(f_drift_TPC,"linear fit inside TPC","l");
    l_drift->SetTextSize(0.043);
    l_drift->Draw();

    TCanvas* c_rasd = new TCanvas("","RASD",g_cwidth,g_cheight);
    TGraph* xz_original = new TGraph();
    TGraph* xz_reconstructed = new TGraph();
    TGraph* xz_reco_map = new TGraph();
    TGraph* g_res = new TGraph();
    TGraph* g_res_map = new TGraph();
    std::vector<X17::Vector> residuals, oldvsnew;
    double max_residual = -1;
    double max_oldvsnew = -1;
    double x_limit = 2; //newcoords ? -1 : 2;
    for (const auto& point : track.points)
    {
        using namespace X17::constants;
        double z_end = newcoords ? point.end.z() : point.end.y();
        if(z_end > 7.5)
        {
            X17::Vector orig,reco_rasd;
            X17::RecoPoint reco_map, reco_map_old;

            if (newcoords)
            {
                orig = {point.start.x(), point.start.y(),point.start.z()};
                reco_rasd = {point.end.x(), point.end.y(), 8-(b0+b1*point.end.t/1000)};
                if (orig.x > x_limit && orig.x < map->GetXMax())
                {
                    reco_map_old = X17::ReconstructOld(*map, point.end.x(),point.end.y(),point.end.t,1E-9,false);
                    reco_map = X17::Reconstruct(*map,point);
                }
            }
            else
            {
                orig = {point.start.z(), point.start.x(),point.start.y()};
                reco_rasd = {point.end.z(), point.end.x(), 8-(b0+b1*point.end.t/1000)};
                if (orig.x > x_limit && orig.x < map->GetXMax())
                {
                    reco_map_old = X17::ReconstructOld(*map, point.end.z(),point.end.x(),point.end.t,1E-9);
                    reco_map = X17::Reconstruct(*map, X17::EndPoint(point.end.z(),point.end.x(),point.end.y(),point.end.t));
                }
            }
            // std::cout << "Original:            (" << orig.x << ", " << orig.y << ", " << orig.z << ")" << std::endl;
            xz_original->AddPoint(orig.x,orig.z);
            xz_reconstructed->AddPoint(reco_rasd.x,reco_rasd.z);
            if (orig.x > x_limit && orig.x < map->GetXMax())
            {
                xz_reco_map->AddPoint(reco_map.x(),reco_map.z());
                g_res->AddPoint(reco_rasd.x,reco_rasd.Dist(orig));
                g_res_map->AddPoint(reco_map.x(),reco_map.AsVector().Dist(orig));

                X17::Vector residual = reco_map.AsVector() - orig;
                if(residual.Magnitude() > max_residual) max_residual = residual.Magnitude();
                residuals.push_back(residual);

                X17::Vector oldvnew = reco_map_old.AsVector() - reco_map.AsVector();
                if(oldvnew.Magnitude() > max_oldvsnew) max_oldvsnew = oldvnew.Magnitude();
                oldvsnew.push_back(oldvnew);
            }
        }
    }
    xz_original->SetMarkerStyle(7);
    xz_original->SetMarkerColor(kRed);
    xz_reconstructed->SetMarkerStyle(2);
    xz_reconstructed->SetMarkerColor(kBlack);
    
    gStyle->SetTitleYOffset(1.05);
    TMultiGraph* mg_rasd = new TMultiGraph();
    mg_rasd->Add(xz_original,"P");
    mg_rasd->Add(xz_reconstructed,"P");
    mg_rasd->GetXaxis()->SetTitle("x [cm]");
    mg_rasd->GetYaxis()->SetTitle("z [cm]");
    mg_rasd->Draw("A");
    
    TLegend* l_rasd = new TLegend(0.69,0.82,0.93,0.93);
    l_rasd->AddEntry(xz_original,"simulation","p");
    l_rasd->AddEntry(xz_reconstructed,"reconstructed","p");
    l_rasd->SetTextSize(0.043);
    l_rasd->Draw();
    
    TCanvas* c_res = new TCanvas("","Residuals",g_cwidth,g_cheight);
    g_res->SetMarkerStyle(2);
    g_res->SetMarkerColor(kBlack);
    g_res->GetXaxis()->SetTitle("x [cm]");
    g_res->GetYaxis()->SetTitle("residual [cm]");
    g_res->Draw("AP");
    g_res_map->SetMarkerStyle(2);
    g_res_map->SetMarkerColor(kRed);
    g_res_map->Draw("P same");
    
    TLegend* l_res = new TLegend(0.17,0.8,0.5,0.91);
    l_res->AddEntry(g_res,"RASD residuals","p");
    l_res->AddEntry(g_res_map,"Map residuals","p");
    l_res->SetTextSize(0.043);
    l_res->Draw();
    
    TCanvas* c_map = new TCanvas("","Reco map",g_cwidth,g_cheight);
    xz_reco_map->SetMarkerStyle(2);
    xz_reco_map->SetMarkerColor(kBlack);
    TMultiGraph* mg_map = new TMultiGraph();
    mg_map->Add(xz_original,"P");
    mg_map->Add(xz_reco_map,"P");
    mg_map->GetXaxis()->SetTitle("x [cm]");
    mg_map->GetYaxis()->SetTitle("z [cm]");
    mg_map->Draw("A");

    TLegend* l_map = new TLegend(0.69,0.82,0.93,0.93);
    l_map->AddEntry(xz_original,"simulation","p");
    l_map->AddEntry(xz_reco_map,"reconstructed","p");
    l_map->SetTextSize(0.043);
    l_map->Draw();
    
    double mag_sigma;
    X17::Vector sigma = GetStDev(residuals,&mag_sigma);
    double hmin = -1.2*max_residual;
    double hmax =  1.2*max_residual;

    TH1F* h_xres = new TH1F("",";x residual [cm];# of electrons",GetBinsScott(hmin,hmax,sigma.x,residuals.size()),hmin,hmax);
    TH1F* h_yres = new TH1F("",";y residual [cm];# of electrons",GetBinsScott(hmin,hmax,sigma.y,residuals.size()),hmin,hmax);
    TH1F* h_zres = new TH1F("",";z residual [cm];# of electrons",GetBinsScott(hmin,hmax,sigma.z,residuals.size()),hmin,hmax);
    TH1F* h_res = new TH1F("",";residual [cm];# of electrons",GetBinsScott(0,hmax,mag_sigma,residuals.size()),0,hmax);

    for (X17::Vector residual : residuals)
    {
        h_xres->Fill(residual.x);
        h_yres->Fill(residual.y);
        h_zres->Fill(residual.z);
        h_res->Fill(residual.Magnitude());
    }
    
    TCanvas* c2 = new TCanvas("","Residual histograms",2*g_cwidth,2*g_cheight);
    c2->Divide(2,2);
    c2->cd(1); h_xres->Draw("hist");
    c2->cd(2); h_yres->Draw("hist");
    c2->cd(3); h_zres->Draw("hist");
    c2->cd(4); h_res->Draw("hist");
    

    double mag_sigma2;
    X17::Vector sigma2 = GetStDev(oldvsnew,&mag_sigma2);
    double hmin2 = -1.2*max_oldvsnew;
    double hmax2 =  1.2*max_oldvsnew;

    TH1F* h_yres2 = new TH1F("","Y residuals;#Delta y [cm];# of electrons",GetBinsScott(hmin2,hmax2,sigma2.y,oldvsnew.size()),hmin2,hmax2);
    TH1F* h_xres2 = new TH1F("","X residuals;#Delta x [cm];# of electrons",GetBinsScott(hmin2,hmax2,sigma2.x,oldvsnew.size()),hmin2,hmax2);
    TH1F* h_zres2 = new TH1F("","Z residuals;#Delta z [cm];# of electrons",GetBinsScott(hmin2,hmax2,sigma2.z,oldvsnew.size()),hmin2,hmax2);
    TH1F* h_res2 = new TH1F("","Residuals;Deviation [cm];# of electrons",GetBinsScott(0,hmax2,mag_sigma2,oldvsnew.size()),0,hmax2);

    for (X17::Vector residual : oldvsnew)
    {
        h_yres2->Fill(residual.y);
        h_xres2->Fill(residual.x);
        h_zres2->Fill(residual.z);
        h_res2->Fill(residual.Magnitude());
    }
    
    TCanvas* c3 = new TCanvas("","Residual histograms old vs new reco",g_cwidth,g_cheight);
    c3->Divide(2,2);
    c3->cd(2); h_yres2->Draw("hist");
    c3->cd(1); h_xres2->Draw("hist");
    c3->cd(3); h_zres2->Draw("hist");
    c3->cd(4); h_res2->Draw("hist");

    if (!g_g_res1) g_g_res1 = g_res;
    else if (!g_g_res2) g_g_res2 = g_res;
}

void PlotRASDres2()
{
    TCanvas* c = new TCanvas("","Residuals comparison",g_cwidth,g_cheight);
    g_g_res2->SetMarkerColor(kBlue);

    gStyle->SetTitleYOffset(0.9);

    TMultiGraph* mg_res = new TMultiGraph();
    mg_res->Add(g_g_res1,"p");
    mg_res->Add(g_g_res2,"p");
    mg_res->GetXaxis()->SetTitle("x [cm]");
    mg_res->GetYaxis()->SetTitle("residual [cm]");
    mg_res->Draw("A");

    TLegend* l_res = new TLegend(0.17,0.775,0.4,0.91);
    l_res->AddEntry(g_g_res1,"90:10 Ar:CO_{2}","p");
    l_res->AddEntry(g_g_res2,"70:30 Ar:CO_{2}","p");
    l_res->SetTextSize(0.043);
    l_res->Draw();
}

void PlotPadReco(const X17::Field<X17::MapPoint>& map, bool m7030 = true)
{
    int pad_id = 127;//66;
    double x1,y1,x2,y2,xc,yc;

    X17::DefaultLayout::GetDefaultLayout().GetPadCorners(pad_id,x1,y1,x2,y2,true);

    double max_time = m7030 ? 17050 : 4750;
    double time_step = 10;
    int i_max = max_time / time_step;
    double time_step2 = 100;
    int plane_steps = 20;
    int i_max2 = max_time / time_step2;
    X17::Vector corners[4] = {
        {x1, y1, -8},
        {x1, y2, -8},
        {x2, y1, -8},
        {x2, y2, -8}
    };
    std::vector<TPolyLine3D*> pad_lines;
    TGraph2D* g_center = new TGraph2D();
    TGraph2D* g_map_pts = new TGraph2D();
    X17::DefaultLayout::GetDefaultLayout().GetPadCenter(pad_id,xc,yc);

    double x_min = 1e6;
    double x_max = -1e6;
    double y_min = 1e6;
    double y_max = -1e6;

    // Corner lines
    for (auto v_corner : corners)
    {
        TPolyLine3D* corner_line = new TPolyLine3D(i_max + 1);

        for (int i = 0; i < i_max; i++)
        {
            double time = i * time_step;
            X17::EndPoint t_corner (v_corner, time);
            X17::RecoPoint reco = X17::Reconstruct(map,t_corner);
            corner_line->SetPoint(i,reco.x(),reco.y(),reco.z());

            if (reco.x() < x_min) x_min = reco.x();
            if (reco.x() > x_max) x_max = reco.x();
            if (reco.y() < y_min) y_min = reco.y();
            if (reco.y() > y_max) y_max = reco.y();
        }

        pad_lines.push_back(corner_line);
    }

    // Plane lines and centers
    for (int i = 0; i < i_max2; i++)
    {
        double time = i * time_step2;

        X17::RecoPoint creco = X17::Reconstruct(map,{xc,yc,-8,time+time_step2/2.},g_map_pts);
        g_center->AddPoint(creco.x(),creco.y(),creco.z());

        TPolyLine3D* lines[4] = {
            new TPolyLine3D (plane_steps),
            new TPolyLine3D (plane_steps),
            new TPolyLine3D (plane_steps),
            new TPolyLine3D (plane_steps)
        };

        for (int j = 0; j < plane_steps; j++)
        {
            std::vector<X17::RecoPoint> points = {
                X17::Reconstruct(map,{x1 + j*(x2 - x1) / (plane_steps - 1), y1, -8, time}),
                X17::Reconstruct(map,{x2, y1 + j*(y2 - y1) / (plane_steps - 1), -8, time}),
                X17::Reconstruct(map,{x1 + j*(x2 - x1) / (plane_steps - 1), y2, -8, time}),
                X17::Reconstruct(map,{x1, y1 + j*(y2 - y1) / (plane_steps - 1), -8, time})
            };

            for (int k = 0; k < 4; k++)
                lines[k]->SetPoint(j,points[k].x(),points[k].y(),points[k].z());
        }

        for (int l = 0; l < 4; l++)
            pad_lines.push_back(lines[l]);        
    }

    x_min -= 0.1;
    x_max += 0.1;

    double y_avg = (y_min+y_max)/2;
    y_max = y_avg + (x_max-x_min)/2;
    y_min = y_avg - (x_max-x_min)/2;

    TCanvas* c = new TCanvas("","Pad inversion",g_cwidth,g_cheight);
    ApplyThesisStyle(c);
    TH3F* h = new TH3F("",";x [cm];y [cm];z [cm]", 1, x_min, x_max, 1, y_min, y_max, 1, -8, 8);
    ApplyThesisStyle(h);
    h->Draw();
    for (auto line : pad_lines)
        line->Draw("same");
    g_center->SetMarkerStyle(8);
    g_center->SetMarkerColor(kRed);
    g_center->Draw("P same");
    g_map_pts->SetMarkerStyle(8);
    g_map_pts->SetMarkerColor(kBlue);
    g_map_pts->Draw("P same");
}

void PlotSplineAndCirc2D(X17::TrackMicro track, X17::Field<X17::Vector>* magfield, X17::Field<X17::MapPoint>* map, bool newcoords = false)
{
    TGraph* xz_original = new TGraph();
    TGraph* xz_reco_map = new TGraph();
    TGraph* xz_reco_map_copy = new TGraph();
    double x_limit = 2; //newcoords ? -1 : 2;
    for (const auto& point : track.points)
    {
        using namespace X17::constants;
        double z_end = newcoords ? point.end.z() : point.end.y();
        if (z_end > 7.5)
        {
            X17::Vector orig,reco_rasd;
            X17::RecoPoint reco_map, reco_map_old;

            if (newcoords)
            {
                orig = {point.start.x(), point.start.y(),point.start.z()};
                if (orig.x > x_limit && orig.x < map->GetXMax())
                {
                    reco_map_old = X17::ReconstructOld(*map, point.end.x(),point.end.y(),point.end.t,1E-9,false);
                    reco_map = X17::Reconstruct(*map,point);
                }
            }
            else
            {
                orig = {point.start.z(), point.start.x(),point.start.y()};
                if (orig.x > x_limit && orig.x < map->GetXMax())
                {
                    reco_map_old = X17::ReconstructOld(*map, point.end.z(),point.end.x(),point.end.t,1E-9);
                    reco_map = X17::Reconstruct(*map, X17::EndPoint(point.end.z(),point.end.x(),point.end.y(),point.end.t));
                }
            }
            // std::cout << "Original:            (" << orig.x << ", " << orig.y << ", " << orig.z << ")" << std::endl;
            xz_original->AddPoint(orig.x,orig.z);
            if (orig.x > x_limit && orig.x < map->GetXMax())
            {
                xz_reco_map->AddPoint(reco_map.x(),reco_map.z());
                xz_reco_map_copy->AddPoint(reco_map.x(),reco_map.z());
            }
        }
    }
    xz_original->SetMarkerStyle(7);
    xz_original->SetMarkerColor(kRed);
    
    gStyle->SetTitleYOffset(1.05);
    
    TSpline3* spline = X17::FitSplines<4>(xz_reco_map,x_limit,map->GetXMax());
    TF1* f_spline = xz_reco_map->GetFunction("fit");
    f_spline->SetLineColor(kBlue);

    TCanvas* c_map = new TCanvas("","Spline (map)",g_cwidth,g_cheight);
    xz_reco_map->SetMarkerStyle(2);
    xz_reco_map->SetMarkerColor(kBlack);
    xz_reco_map->SetLineColor(kBlue);

    TMultiGraph* mg_map = new TMultiGraph();
    mg_map->Add(xz_original,"P");
    mg_map->Add(xz_reco_map,"P");
    mg_map->GetXaxis()->SetTitle("x [cm]");
    mg_map->GetYaxis()->SetTitle("z [cm]");
    mg_map->Draw("A");

    TLegend* l_map = new TLegend(0.69,0.765,0.93,0.93);
    l_map->AddEntry(xz_original,"simulation","p");
    l_map->AddEntry(xz_reco_map,"reconstructed","p");
    l_map->AddEntry(f_spline,"spline fit","l");
    l_map->SetTextSize(0.043);
    l_map->Draw();

    int xsteps = 200;
    TGraph* g_magfield = new TGraph();
    TGraph* g_radius = new TGraph();
    TGraph* g_energy = new TGraph();

    for (double x = x_limit; x < map->GetXMax(); x += (map->GetXMax()-x_limit)/xsteps)
    {
        double z = spline->Eval(x);
        X17::Vector bfield = magfield->GetField(x,0,z);
        g_magfield->AddPoint(x,bfield.y);

        int i_node = spline->FindX(x);

        double xnode,ynode,b,c,d;
        spline->GetCoeff(i_node,xnode,ynode,b,c,d);
        double dx   = x-xnode;
        double der  = b+dx*(2*c+3*d*dx);
        double der2 = 2*c+6*d*dx;
        double r = std::abs(0.01*pow(1+der*der,1.5)/der2);
        const double clight = 299792458;
        const double E0 = 510998.95;
        double betasq = 1/(1+pow((E0/(clight*r*bfield.y)),2));
        double Ekin = E0*(1/sqrt(1-betasq)-1);
        if (x > 4 && r < 1.5 && r > 0)
        {
            g_radius->AddPoint(x,r*100);
            g_energy->AddPoint(x,Ekin/1e6);
        }
    }
    
    TCanvas* c_magfield = new TCanvas("","Magnetic field",g_cwidth,g_cheight);
    g_magfield->SetMarkerStyle(2);
    g_magfield->SetMarkerColor(kBlack);
    g_magfield->GetXaxis()->SetTitle("x [cm]");
    g_magfield->GetYaxis()->SetTitle("B [T]");
    g_magfield->Draw("AP");

    TCanvas* c_radius = new TCanvas("","Reconstructed Radius",g_cwidth,g_cheight);
    g_radius->SetMarkerStyle(2);
    g_radius->SetMarkerColor(kBlack);
    g_radius->GetXaxis()->SetTitle("x [cm]");
    g_radius->GetYaxis()->SetTitle("r [cm]");
    g_radius->Draw("AP");
    
    TCanvas* c_energy = new TCanvas("","Reconstructed Energy",g_cwidth,g_cheight);
    g_energy->SetMarkerStyle(2);
    g_energy->SetMarkerColor(kBlack);
    g_energy->GetXaxis()->SetTitle("x [cm]");
    g_energy->GetYaxis()->SetTitle("E [MeV]");
    g_energy->Draw("AP");

    TF1* circ2 = X17::FitCircle2(xz_reco_map_copy,x_limit,map->GetXMax());
    TF1* circ = xz_reco_map_copy->GetFunction("circle");
    circ->SetLineColor(kBlue);
    circ->SetLineWidth(2);
    // circ->Draw("same");

    TCanvas* c_map2 = new TCanvas("","Circle 2D (map)",g_cwidth,g_cheight);
    xz_reco_map_copy->SetMarkerStyle(2);
    xz_reco_map_copy->SetMarkerColor(kBlack);
    xz_reco_map_copy->SetLineColor(kBlue);

    TMultiGraph* mg_map2 = new TMultiGraph();
    mg_map2->Add(xz_original,"P");
    mg_map2->Add(xz_reco_map_copy,"P");
    mg_map2->GetXaxis()->SetTitle("x [cm]");
    mg_map2->GetYaxis()->SetTitle("z [cm]");
    mg_map2->Draw("A");

    TLegend* l_map2 = new TLegend(0.69,0.765,0.93,0.93);
    l_map2->AddEntry(xz_original,"simulation","p");
    l_map2->AddEntry(xz_reco_map_copy,"reconstructed","p");
    l_map2->AddEntry(circ,"circle fit","l");
    l_map2->SetTextSize(0.043);
    l_map2->Draw();

    TGraph* g_magfield2 = new TGraph();
    X17::RecoEnergy(circ,*magfield,g_magfield2,x_limit,map->GetXMax(),0.1);

    TCanvas* c_magfield2 = new TCanvas("","Magnetic field",g_cwidth,g_cheight);
    g_magfield2->SetMarkerStyle(2);
    g_magfield2->SetMarkerColor(kBlack);
    g_magfield2->GetXaxis()->SetTitle("x [cm]");
    g_magfield2->GetYaxis()->SetTitle("B [T]");
    g_magfield2->Draw("AP");
}

int main(int argc, char *argv[])
{
    TApplication app("app", &argc, argv);
    std::cout << "Working directory: " << gSystem->WorkingDirectory() << std::endl;

    int orange_index = TColor::GetFreeColorIndex();
    TColor* orange = new TColor(orange_index,1,0.6,0);

    X17::Field<X17::Vector>* magfield = X17::LoadField("../../../data/elmag/VecB2.txt",{-20,-30,-30},{20,30,30},0.5);
    TFile* map_input = new TFile("../../../data/ion_map/sample_1.0/map.root");
    TFile* map_input2 = new TFile("../../../data/ion_map/sample_2.0/map.root");
    X17::Field<X17::MapPoint>* map9010 = (X17::Field<X17::MapPoint>*)map_input->Get("map");
    X17::Field<X17::MapPoint>* map7030 = (X17::Field<X17::MapPoint>*)map_input2->Get("map");

    X17::TrackMicro track1 = LoadTrack("original");
    // PlotDriftXZ(track1);
    // PlotDriftXY(track1,"drift_xy_9010.png");
    // PlotDriftYZ(track1,"drift_yz_9010.png");

    X17::TrackMicro track2 = LoadTrack("new_7030");
    // PlotDriftXZ(track2,true);
    // PlotDriftXY(track2,"drift_xy_7030.png",true);
    // PlotDriftYZ(track2,"drift_yz_7030.png",true);

    gStyle->SetTitleSize(0.06,"XY");
    gStyle->SetLabelSize(0.06,"XY");
    gStyle->SetTitleXOffset(0.9);
    gStyle->SetTitleYOffset(0.9);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadRightMargin(0.07);
    gStyle->SetPadTopMargin(0.07);

    // PlotTrackRK(track1, magfield);
    // PlotTrackRK(track2, magfield, true);
    
    gStyle->SetOptStat(0);

    // PlotRASD(track1,map9010);
    // PlotRASD(track2,map7030,true);
    // PlotRASDres2();

    // PlotPadReco(*map9010,false);
    // PlotPadReco(*map7030,true);

    PlotSplineAndCirc2D(track2,magfield,map7030,true);

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

    app.Run();

    return 0;
}

