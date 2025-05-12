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
#include "TImage.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "../../build/X17_dict.cxx"

// X17 dependencies
#include "CircleFit3D.h"
#include "Color.h"
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

void PlotTrackRK(X17::TrackMicro track, X17::Field<X17::Vector>* magfield)
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
            X17::Vector start_new = {start_old.z,start_old.x,start_old.y};

            X17::Vector closest_point;
            double res2 = X17::GetTrackRKSqDistAndCP(trackrk,start_new,closest_point);
            g_res->AddPoint(start_new.x,10000*std::sqrt(res2));

            X17::Vector exaggerated = closest_point + 1000 * (start_new-closest_point);
            g_micro->AddPoint(exaggerated.x,exaggerated.z);

            X17::Vector diff = start_new - closest_point;
            g_diff->AddPoint(start_new.x,10000*diff.z);
            g_diff2->AddPoint(start_new.x,10000*diff.x);
        }

        TCanvas* c_rk = new TCanvas("c_rk","Drift XZ",g_cwidth,g_cheight);
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

        TCanvas* c_res = new TCanvas("c_res","Residuals",g_cwidth,g_cheight);
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

void PlotRASD(X17::TrackMicro track, bool newcoords = false)
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
    TGraph* g_res = new TGraph();
    for (const auto& point : track.points)
    {
        using namespace X17::constants;
        double z_end = newcoords ? point.end.z() : point.end.y();
        if(z_end > 7.5)
        {
            X17::Vector orig,reco;

            if (newcoords)
            {
                orig = {point.start.x(), point.start.y(),point.start.z()};
                reco = {point.end.x(), point.end.y(), 8-(b0+b1*point.end.t/1000)};
            }
            else
            {
                orig = {point.start.z(), point.start.x(),point.start.y()};
                reco = {point.end.z(), point.end.x(), 8-(b0+b1*point.end.t/1000)};
            }
            xz_original->AddPoint(orig.x,orig.z);
            xz_reconstructed->AddPoint(reco.x,reco.z);
            g_res->AddPoint(reco.x,std::sqrt(reco.SqDist(orig)));//reco.z-orig.z);//reco.x-orig.x);//
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
    g_res->SetMarkerColor(kRed);
    g_res->GetXaxis()->SetTitle("x [cm]");
    g_res->GetYaxis()->SetTitle("residual [cm]");
    g_res->Draw("AP");

    TLegend* l_res = new TLegend(0.17,0.855,0.5,0.91);
    l_res->AddEntry(g_res,"reconstr. residuals","p");
    l_res->SetTextSize(0.043);
    l_res->Draw();

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

int main(int argc, char *argv[])
{
    TApplication app("app", &argc, argv);
    std::cout << "Working directory: " << gSystem->WorkingDirectory() << std::endl;

    int orange_index = TColor::GetFreeColorIndex();
    TColor* orange = new TColor(orange_index,1,0.6,0);

    X17::Field<X17::Vector>* magfield = X17::LoadField("../../../data/elmag/VecB2.txt",{-20,-30,-30},{20,30,30},0.5);

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

    PlotRASD(track1);
    PlotRASD(track2,true);
    PlotRASDres2();

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

