#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TTree.h>

#include <iostream>

#include "../NSplines.h"

using namespace std;

int reco_track()
{
    TFile* inFile = new TFile("build/electrons.root");
    TTree* electrons = (TTree*)inFile->Get("electrons");

    TCanvas* c1 = new TCanvas();
    electrons->Draw("t1:8-y0","y1>7.0");
    TGraph* ty = new TGraph(electrons->GetSelectedRows(), electrons->GetV2(), electrons->GetV1());
    ty->Draw("ap");
    ty->Fit("pol1","","",8.05,11);

    TF1* ty_fit = (TF1*) ty->GetListOfFunctions()->FindObject("pol1");
    double a0 = ty_fit->GetParameter(0);
    double a1 = ty_fit->GetParameter(1);
    double b0 = -a0/a1; //inverse polynomial param
    double b1 = 1.0/a1; //inverse polynomial param

    double x0,y0,z0,t0,x1,y1,z1,t1;
    electrons->SetBranchAddress("x0",&x0);
    electrons->SetBranchAddress("y0",&y0);
    electrons->SetBranchAddress("z0",&z0);
    electrons->SetBranchAddress("t0",&t0);
    electrons->SetBranchAddress("x1",&x1);
    electrons->SetBranchAddress("y1",&y1);
    electrons->SetBranchAddress("z1",&z1);
    electrons->SetBranchAddress("t1",&t1);

    TGraph* zy = new TGraph();
    TGraph* zy_reco = new TGraph();

    for (int i = 0; i < electrons->GetEntries(); ++i)
    {
        electrons->GetEntry(i);
        if (y1 > 7.0) 
        {
            zy->AddPoint(z0,8-y0);
            zy_reco->AddPoint(z1,b0+b1*t1);
        }
    }

    TCanvas* c2 = new TCanvas("c2","Electron track reconstruction");
    
    zy_reco->SetTitle("Electron track reconstruction;z [cm]; distance to readout [cm]");
    zy_reco->SetMarkerStyle(2);
    zy_reco->SetMarkerSize(0.4);
    zy_reco->Draw("ap");

    zy->SetMarkerColor(2);
    zy->SetMarkerStyle(7);
    zy->SetMarkerSize(1.2);
    zy->Draw("p same");

    TLegend* leg = new TLegend(0.129,0.786,0.360,0.887);
    leg->AddEntry(zy,"original");
    leg->AddEntry(zy_reco,"reconstructed");
    leg->Draw("same");

    const int nodes = 5;
    double min = 0;
    double max = 15;
    NSpline<nodes>* s = new NSpline<nodes>(min,max);
    s->SetDer1(0);

    TF1* zy_fit = new TF1("zz_fit",s->GetEval(),min,max,2*nodes+2);
    zy->Fit(zy_fit,"M","",min,max);

    TH1F* residues = new TH1F("h_residues","Residues",50,-0.5,0.5);

    for (int i = 0; i < zy_reco->GetN(); i++)
    {
        double x,y;
        zy_reco->GetPoint(i,x,y);
        residues->Fill(y-zy_fit->Eval(x));
    }

    TCanvas* c3 = new TCanvas("c3","Residues");
    residues->Draw();    

    return 0;
}

int main()
{
    return reco_track();
}

// plotting into 2D histogram--------------------------------------------------------
    // TH2F* h = new TH2F("h","h",350,7.5,11,500,2300,3300);

    // TTree* electrons = (TTree*)inFile->Get("electrons");

    // electrons->Draw("t1:8-y0>>h","y1>7.0","colz");
    // h->Fit("pol1","","",8.1,11);

//spline fit of zz graph-------------------------------------------------------------
    // TCanvas* c2 = new TCanvas();
    // electrons->Draw("z0-z1:z1","y1>7.0");
    // TGraph* zz = new TGraph(electrons->GetSelectedRows(), electrons->GetV2(), electrons->GetV1());
    // zz->Draw("ap");

    // const int nodes = 10;
    // double min = 0;
    // double max = 15;
    // NSpline<nodes>* s = new NSpline<nodes>(min,max);
    // s->SetDer1(0);

    // TF1* zz_fit = new TF1("zz_fit",s->GetEval(),min,max,2*nodes+2);
    // zz->Fit(zz_fit,"M","",min,max);