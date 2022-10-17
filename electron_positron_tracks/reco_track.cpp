#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TTree.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "../NSplines.h"
#include "../VectorField.h"

TF1* FitCircle(TGraph* graph, const double& min, const double& max)
{
    TF1* circle = new TF1("circle","[2]-sqrt([0]^2-(x-[1])^2)",min,max);
    circle->SetParameter(0,5);
    circle->SetParameter(1,(min+max)/2);
    circle->SetParameter(2,13);
    graph->Fit(circle,"M","",min,max);
    return circle;
}

double circle_func(double* x, double* par)
{
    double& xx       = x[0];
    double& radius   = par[0];
    double& node1_x  = par[1];
    double& node2_x  = par[2];
    double& node1_y  = par[3];
    double& node2_y  = par[4];

    double x_mid = (node2_x-node1_x)/2;
    double y_mid = (node2_y-node1_y)/2;
    double r_mid = sqrt(pow(x_mid,2)+pow(y_mid,2));

    double r_c = sqrt(pow(radius,2)-pow(r_mid,2));
    double x0 = node1_x+x_mid-y_mid*r_c/r_mid;
    double y0 = node1_y+y_mid+x_mid*r_c/r_mid;

    double a1 = (node1_x-x0)/sqrt(pow(radius,2)-pow(node1_x-x0,2));
    double a2 = (node2_x-x0)/sqrt(pow(radius,2)-pow(node2_x-x0,2));
    double b1 = node1_y-a1*node1_x;
    double b2 = node2_y-a2*node2_x;

    if (xx < node1_x) return a1*xx+b1;
    if (xx > node2_x) return a2*xx+b2;    
    return y0 - sqrt(pow(radius,2)-pow(xx-x0,2));
}

TF1* FitCircle2(TGraph* graph, const double& min, const double& max)
{
    TF1* circle = new TF1("circle",circle_func,min,max,5);
    circle->SetParameter(0,50);
    circle->SetParameter(1,4);
    circle->SetParameter(2,12);
    circle->SetParameter(3,8);
    circle->SetParameter(4,10);
    // circle->Draw();
    graph->Fit(circle,"M","",min,max);
    return circle;
}

void RecoEnergy(TSpline3* sp_fit, VectorField* magfield, TGraph* energy, TGraph* radius, TGraph* magnetic, double min, double max, double step)
{
    for (double x = min; x <= max; x += step)
    {
        int i_node = sp_fit->FindX(x);

        double xnode,ynode,b,c,d;
        sp_fit->GetCoeff(i_node,xnode,ynode,b,c,d);
        double dx   = x-xnode;
        double der  = b+dx*(2*c+3*d*dx);
        double der2 = 2*c+6*d*dx;
        double r = 0.01*pow(1+der*der,1.5)/der2;
        const double clight = 299792458;
        const double E0 = 510998.95;
        Vector B = magfield->GetField(0,(8-sp_fit->Eval(x))/100,x/100);
        double betasq = 1/(1+pow((E0/(clight*r*B.vx)),2));
        double Ekin = E0*(1/sqrt(1-betasq)-1);
        if (x > 4) energy->AddPoint(x,Ekin/1e6);
        if (r < 1 && r > 0) radius->AddPoint(x,r*100);
        magnetic->AddPoint(x,B.vx);
    }
}

void RecoEnergy(TF1* fit, VectorField* magfield, TGraph* magnetic, double min, double max, double step)
{
    double r = fit->GetParameter(0)/100;
    const double clight = 299792458;
    const double E0 = 510998.95;

    //mean magnetic field
    Vector B = magfield->GetField(0,(8-fit->Eval((max+min)/2))/100,(max+min)/200);

    for (double x = min; x <= max; x += step)
    {
        Vector B2 = magfield->GetField(0,(8-fit->Eval(x))/100,x/100);
        magnetic->AddPoint(x,B2.vx);
    }
    
    double betasq = 1/(1+pow((E0/(clight*r*B.vx)),2));
    double Ekin = E0*(1/sqrt(1-betasq)-1);

    cout << "Kinetic energy: " << Ekin << " eV, Bx: " << B.vx << " T, beta: ";
    cout << sqrt(betasq) << "\n";
}

// estimates difference between two points with given values in map
double Offset(SensorData s, double x1, double z1, double t1)
{
    constexpr double tfact = 0.00327; // time is measured at different scale, it needs weight
    return sqrt(pow(x1-s.x1,2)+pow(z1-s.z1,2)+pow(tfact*(t1-s.t1),2));
}

SensorData RecoPoint(Field<SensorData>* map, double x1, double z1, double t1, double max_err)
{
    // start looking at the same position
    double x = x1;
    double y = (map->ymax+map->ymin)/2;
    double z = z1;
    double step = map->step/10;

    double offset;      // metric of distance between points
    int iterations = 0; // number of iterations should not exceed 100
    double damp = 0.1;  // damping coefficient

    // loop for offset minimization
    do
    {
        // calculate offset gradient
        SensorData xa = map->GetField(x+step,y,z);
        SensorData xb = map->GetField(x-step,y,z);
        SensorData ya = map->GetField(x,y+step,z);
        SensorData yb = map->GetField(x,y-step,z);
        SensorData za = map->GetField(x,y,z+step);
        SensorData zb = map->GetField(x,y,z-step);

        double oxa = Offset(xa,x1,z1,t1);
        double oxb = Offset(xb,x1,z1,t1);
        double oya = Offset(ya,x1,z1,t1);
        double oyb = Offset(yb,x1,z1,t1);
        double oza = Offset(za,x1,z1,t1);
        double ozb = Offset(zb,x1,z1,t1);

        double gradx = (oxa-oxb)/(2*step);
        double grady = (oya-oyb)/(2*step);
        double gradz = (oza-ozb)/(2*step);

        //adjust current guess by minus gradient
        x -= damp*gradx; y -= damp*grady; z -= damp*gradz;

        //check bounds
        if (x < map->xmin) x = map->xmin;
        if (x > map->xmax) x = map->xmax;
        if (y < map->ymin) y = map->ymin;
        if (y > map->ymax) y = map->ymax;
        if (z < map->zmin) z = map->zmin;
        if (z > map->zmax) z = map->zmax;
        //cout << "gradx: " << x << " grady: " << y << " gradz: " << z << "\n";

        // calculate values at current position
        SensorData cur = map->GetField(x,y,z);
        offset = Offset(cur,x1,z1,t1);

        //make sure step isn't too high
        if(offset < 10*step) step /= 10;

        iterations++;
        // cout << "iter: " << iterations << "\n";        
        if (iterations == 1000) cout << "1000 iterations.\n";
    }
    while ((offset > max_err) && (iterations < 1000));

    return {x,y,z,0};
}

int reco_track()
{
    //file with Garfield simulation output
    TFile* inFile = new TFile("build/electrons.root");
    TTree* electrons = (TTree*)inFile->Get("electrons");

    //file with ionization electrons map
    TFile* inFile2 = new TFile("map.root");
    Field<SensorData>* map = (Field<SensorData>*)inFile2->Get("map");

    //plotting drift time vs distance to readout + linear fit
    TCanvas* c_drift = new TCanvas("c_drift","Drift time");
    electrons->Draw("t1:8-y0","y1>7.0");
    TGraph* ty = new TGraph(electrons->GetSelectedRows(), electrons->GetV2(), electrons->GetV1());
    ty->SetTitle("Drift time as function of distance;distance to readout [cm]; time [ns]");
    ty->SetMarkerStyle(2);
    ty->SetMarkerSize(0.4);
    ty->Draw("ap");
    ty->Fit("pol1","","",8.05,11);

    //getting fit parameters
    TF1* ty_fit = (TF1*) ty->GetListOfFunctions()->FindObject("pol1");
    double a0 = ty_fit->GetParameter(0);
    double a1 = ty_fit->GetParameter(1);
    double b0 = -a0/a1; //inverse polynomial param
    double b1 = 1.0/a1; //inverse polynomial param

    //setting variables from TTree
    double x0,y0,z0,t0,x1,y1,z1,t1;
    electrons->SetBranchAddress("x0",&x0);
    electrons->SetBranchAddress("y0",&y0);
    electrons->SetBranchAddress("z0",&z0);
    electrons->SetBranchAddress("t0",&t0);
    electrons->SetBranchAddress("x1",&x1);
    electrons->SetBranchAddress("y1",&y1);
    electrons->SetBranchAddress("z1",&z1);
    electrons->SetBranchAddress("t1",&t1);

    //zy (track) plot from original track and reconstructed from drift time
    TGraph* zy = new TGraph();
    TGraph* zy_reco = new TGraph();

    for (int i = 0; i < electrons->GetEntries(); ++i)
    {
        cout << "i: " << i << " out of " << electrons->GetEntries() << "\n";
        //if ((10000*i)%electrons->GetEntries() == 0) cout << 100*i/electrons->GetEntries() << " \%\n";
        electrons->GetEntry(i);
        if (y1 > 7.0) 
        {
            zy->AddPoint(z0,8-y0);
            SensorData reco = RecoPoint(map,x1,z1,t1,0.001);
            zy_reco->AddPoint(reco.z1,8-reco.y1);
            //zy_reco->AddPoint(z1,b0+b1*t1);
        }
    }

    //setting up track plots (original + reconstructed)
    TCanvas* c_track = new TCanvas("c_track","Electron track reconstruction");
    
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

    double min = 0;//7
    double max = 15;//10.5

    //fitting both tracks with circles
    TF1* circle_fit  = FitCircle2(zy,min,max);
    TF1* circle_fit2 = FitCircle2(zy_reco,min,max);

    //loading magnetic field from txt file (units = meters)
    VectorField* magfield = new VectorField(-0.3,0.3,-0.3,0.3,-0.2,0.2,0.005);
    magfield->LoadField("/home/vavrik/work/X17/electron_positron_tracks/build/VecB.txt");
    
    TGraph* magnetic_x  = new TGraph();
    TGraph* magnetic_x2 = new TGraph();

    double step = 0.1;
    RecoEnergy(circle_fit,magfield,magnetic_x,min,max,step);
    RecoEnergy(circle_fit2,magfield,magnetic_x2,min,max,step);
    
    // TCanvas* c_magnetic = new TCanvas("c_magnetic","Perpendicular magnetic field");
    // magnetic_x->SetTitle("Perpendicular magnetic field;z [cm]; B [T]");
    // magnetic_x->SetMarkerStyle(2);
    // magnetic_x->SetMarkerSize(0.4);
    // magnetic_x->SetMarkerColor(2);
    // magnetic_x->Draw("ap");

    // magnetic_x2->SetMarkerStyle(2);
    // magnetic_x2->SetMarkerSize(0.4);
    // magnetic_x2->Draw("p same");

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

// spline fit of zz graph------------------------------------------------------------
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

// plotting magnetic field magnitude on XZ plane
    // TCanvas* c3 = new TCanvas("c3","Residues");
    // residues->Draw();

    // TH2F* field = new TH2F("h_field","XZ magnetic field plot",121,-0.3,0.3,81,-0.2,0.2);

    // for(int i = 0; i < 121; i++)
    //     for(int j = 0; j < 81; j++)
    //     {
    //         Vector b = magfield->GetField(-0.3+i*0.005,0.0,-0.2+j*0.005);//magfield->field[i][60][j];//
    //         double bmag = sqrt(b.vx*b.vx+b.vy*b.vy+b.vz*b.vz);
    //         field->SetBinContent(i,j,bmag);
    //     }
    
    // TCanvas* c4 = new TCanvas("c4","Magnetic field XZ plane");
    // field->Draw("cont1");

// magnetic field plot No. 2
    // double xbins = 121;
    // double zbins = 81;

    // TH2F* field = new TH2F("h_field","XZ magnetic field plot",xbins,-0.3,0.3,zbins,-0.2,0.2);

    // for(int i = 0; i < xbins; i++)
    //     for(int j = 0; j < zbins; j++)
    //     {
    //         Vector b = magfield->GetField(-0.3+i*0.005,0.0,-0.2+j*0.005);//magfield->field[i][60][j];//
    //         double bmag = sqrt(b.vx*b.vx+b.vy*b.vy+b.vz*b.vz);
    //         field->SetBinContent(i,j,b.vz);
    //     }
    
    // TCanvas* c4 = new TCanvas("c4","Magnetic field XZ plane");
    // gStyle->SetPalette(1.);
    // field->Draw("colz");

// plot residues of reconstructed track
    // TH1F* residues = new TH1F("h_residues","Residues",50,-0.5,0.5);

    // for (int i = 0; i < zy_reco->GetN(); i++)
    // {
    //     double x,y;
    //     zy_reco->GetPoint(i,x,y);
    //     residues->Fill(y-sp_fit->Eval(x));
    // }