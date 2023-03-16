#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPolyLine3D.h>
#include <TTree.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "../NSplines.h"
#include "../VectorField.h"
#include "../X17Utilities.h"

/// @brief Function for simple circlular arc fit (using function for half of circle, curves up)
/// @param graph Graph to be fitted
/// @param min Lower bound of the fit
/// @param max Upper bound of the fit
/// @return The fitted function
TF1* FitCircle(TGraph* graph, const double& min, const double& max)
{
    TF1* circle = new TF1("circle","[2]-sqrt([0]^2-(x-[1])^2)",min,max);
    circle->SetParameter(0,5);
    circle->SetParameter(1,(min+max)/2);
    circle->SetParameter(2,13);
    graph->Fit(circle,"M","",min,max);
    return circle;
}

/// @brief Function for circular arc with smoothly attached lines at endpoints (nodes)
/// @param x Variable x
/// @param par Set of parameters (radius of the circle, 1st and 2nd node x and y coordinates)
/// @return The value at x
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

/// @brief Function for fitting with circular arc with smoothly attached lines at endpoints
/// @param graph Graph to be fitted
/// @param min Lower bound of the fit
/// @param max Upper bound of the fit
/// @return The fitted function
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

/// @brief Function for energy reconstruction from spline fit
/// @param sp_fit Fitted splines
/// @param magfield Magnetic data
/// @param energy Output graph for reconstructed energy as function of coordinate
/// @param radius Output graph for reconstructed radius as function of coordinate
/// @param magnetic Output graph for magnetic field along fitted trajectory
/// @param min Lower bound
/// @param max Upper bound
/// @param step Step between iterations
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
        Vector B = magfield->GetField(x/100,0,(8-sp_fit->Eval(x))/100);
        double betasq = 1/(1+pow((E0/(clight*r*B.vx)),2));
        double Ekin = E0*(1/sqrt(1-betasq)-1);
        if (x > 4) energy->AddPoint(x,Ekin/1e6);
        if (r < 1 && r > 0) radius->AddPoint(x,r*100);
        magnetic->AddPoint(x,B.vy);
    }
}

/// @brief Function for energy reconstruction from circle with lines fit
/// @param fit Fitted circle with lines function
/// @param magfield Magnetic data
/// @param magnetic Output graph for magnetic field along fitted trajectory
/// @param min Lower bound
/// @param max Upper bound
/// @param step Step between iterations
void RecoEnergy(TF1* fit, VectorField* magfield, TGraph* magnetic, double min, double max, double step)
{
    double r = fit->GetParameter(0)/100;
    const double clight = 299792458;
    const double E0 = 510998.95;

    // mean magnetic field
    Vector B = magfield->GetField((max+min)/200,0,(8-fit->Eval((max+min)/2))/100);

    for (double x = min; x <= max; x += step)
    {
        Vector B2 = magfield->GetField(x/100,0,(8-fit->Eval(x))/100);
        magnetic->AddPoint(x,B2.vy);
    }
    
    double betasq = 1/(1+pow((E0/(clight*r*B.vy)),2));
    double Ekin = E0*(1/sqrt(1-betasq)-1);

    cout << "Kinetic energy: " << Ekin << " eV, By: " << B.vy << " T, beta: ";
    cout << sqrt(betasq) << "\n";
}

int reco_track()
{
    // file with Garfield simulation output
    TFile* inFile = new TFile("build/electrons.root");
    TTree* electrons = (TTree*)inFile->Get("electrons");

    // file with ionization electrons map
    TFile* inFile2 = new TFile("map.root");
    Field<SensorData>* map = (Field<SensorData>*)inFile2->Get("map");

    // plotting drift time vs distance to readout + linear fit
    TCanvas* c_drift = new TCanvas("c_drift","Drift time");
    electrons->Draw("t1:8-z0","z1>7.0");
    TGraph* tz = new TGraph(electrons->GetSelectedRows(), electrons->GetV2(), electrons->GetV1());
    tz->SetTitle("Drift time as function of distance;distance to readout [cm]; time [ns]");
    tz->SetMarkerStyle(2);
    tz->SetMarkerSize(0.4);
    tz->Draw("ap");
    tz->Fit("pol1","","",8.05,11);

    // getting fit parameters
    TF1* tz_fit = (TF1*) tz->GetListOfFunctions()->FindObject("pol1");
    double a0 = tz_fit->GetParameter(0);
    double a1 = tz_fit->GetParameter(1);
    double b0 = -a0/a1; //inverse polynomial param
    double b1 = 1.0/a1; //inverse polynomial param

    // setting variables from TTree
    double x0,y0,z0,t0,x1,y1,z1,t1;
    electrons->SetBranchAddress("x0",&x0);
    electrons->SetBranchAddress("y0",&y0);
    electrons->SetBranchAddress("z0",&z0);
    electrons->SetBranchAddress("t0",&t0);
    electrons->SetBranchAddress("x1",&x1);
    electrons->SetBranchAddress("y1",&y1);
    electrons->SetBranchAddress("z1",&z1);
    electrons->SetBranchAddress("t1",&t1);

    // xz (track) plot from original track and reconstructed from drift time
    TGraph* xz      = new TGraph();
    TGraph* xz_reco = new TGraph();

    // xy (track) plot from original and reconstructed tracks
    TGraph* xy      = new TGraph();
    TGraph* xy_reco = new TGraph();

    // reconstruction residuals
    TGraph* gx_res = new TGraph();
    TGraph* gy_res = new TGraph();
    TGraph* gz_res = new TGraph();
    TGraph* gr_res = new TGraph();

    TH1F* hx_res = new TH1F("hx_res","X residuals",25,-0.5,0.5);
    TH1F* hy_res = new TH1F("hx_res","Y residuals",25,-0.5,0.5);
    TH1F* hz_res = new TH1F("hx_res","Z residuals",25,-0.5,0.5);
    TH1F* hr_res = new TH1F("hx_res","Residuals",25,0,1);
    
    // Variables for reconstruction with pads
    constexpr int timebins = 50;
    int padhits[X17::channels][timebins];
    for (int i = 0; i < X17::channels; i++) for (int j = 0; j < timebins; j++) {padhits[i][j] = 0;}

    TGraph2D* g_xyz      = new TGraph2D();
    TGraph2D* g_xyz_reco = new TGraph2D();
    // TH3F* h_xyz_reco = new TH3F("h_xyz_reco","Electron track reconstruction;x [cm];y [cm];z [cm]",)

    // looping through all electrons
    int n_electrons = 0;
    for (int i = 0; i < electrons->GetEntries(); ++i)
    {
        //cout << "\n\ni: " << i << " out of " << electrons->GetEntries() << "\n";
        //if ((10000*i)%electrons->GetEntries() == 0) cout << 100*i/electrons->GetEntries() << " \%\n";
        electrons->GetEntry(i);
        if (X17::IsInSector(x1,y1,0)) 
        {
            n_electrons++;

            // SensorData reco = RecoPoint(map,x1,z1,t1,0.001);
            SensorData reco = map->Invert(x1,y1,t1);
            xz->AddPoint(x0,8-z0);
            xz_reco->AddPoint(reco.x1,8-reco.z1);
            // xz_reco->AddPoint(z1,b0+b1*t1);

            xy->AddPoint(x0,y0);
            xy_reco->AddPoint(reco.x1,reco.y1);

            gx_res->AddPoint(x0,reco.x1-x0);
            gy_res->AddPoint(x0,reco.y1-y0);
            gz_res->AddPoint(x0,reco.z1-z0);
            gr_res->AddPoint(x0,sqrt(pow((reco.x1-x0),2)+pow((reco.y1-y0),2)+pow((reco.z1-z0),2)));

            hx_res->Fill(reco.x1-x0);
            hy_res->Fill(reco.y1-y0);
            hz_res->Fill(reco.z1-z0);
            hr_res->Fill(sqrt(pow((reco.x1-x0),2)+pow((reco.y1-y0),2)+pow((reco.z1-z0),2)));

            // reconstruction with pads
            int channel = X17::GetPad(x1,y1);
            int timebin = t1/100; if(timebin > timebins - 1) cerr << "ERROR: Invalid timebin: " << timebin << endl;
            if(channel == -1) cerr << "ERROR: No pad hit found. Coordinates x,y: " << x1 << ", " << y1 << endl;
            padhits[channel-1][timebin]++;
            g_xyz->AddPoint(x0,y0,z0);
        }

    }
    cout << "\nNumber of electrons in the TPC region: " << n_electrons << "\n";

    // setting up track plots (original + reconstructed)
    TCanvas* c_track_xz = new TCanvas("c_track_xz","Electron track reconstruction");
    
    xz_reco->SetTitle("Electron track reconstruction;x [cm]; distance to readout [cm]");
    xz_reco->SetMarkerStyle(2);
    xz_reco->SetMarkerSize(0.4);
    xz_reco->Draw("ap");

    xz->SetMarkerColor(2);
    xz->SetMarkerStyle(7);
    xz->SetMarkerSize(1.2);
    xz->Draw("p same");

    TLegend* leg_xz = new TLegend(0.129,0.786,0.360,0.887);
    leg_xz->AddEntry(xz,"original");
    leg_xz->AddEntry(xz_reco,"reconstructed");
    leg_xz->Draw("same");

    TCanvas* c_track_xy = new TCanvas("c_track_xy","Electron track reconstruction");
    
    xy_reco->SetTitle("Electron track reconstruction;x [cm]; y [cm]");
    xy_reco->SetMarkerStyle(2);
    xy_reco->SetMarkerSize(0.4);
    xy_reco->Draw("ap");

    xy->SetMarkerColor(2);
    xy->SetMarkerStyle(7);
    xy->SetMarkerSize(1.2);
    xy->Draw("p same");

    TLegend* leg_xy = new TLegend(0.129,0.786,0.360,0.887);
    leg_xy->AddEntry(xy,"original");
    leg_xy->AddEntry(xy_reco,"reconstructed");
    leg_xy->Draw("same");


    TCanvas* c_fit_res = new TCanvas("c_fit_res", "Reconstruction residuals");
    c_fit_res->Divide(2,2);
    
    c_fit_res->cd(1);
    gx_res->SetTitle("X residuals;x [cm];Δx [cm]");
    gx_res->SetMarkerStyle(2);
    gx_res->SetMarkerSize(0.4);
    gx_res->Draw("ap");
    
    c_fit_res->cd(2);
    gy_res->SetTitle("Y residuals;x [cm];Δy [cm]");
    gy_res->SetMarkerStyle(2);
    gy_res->SetMarkerSize(0.4);
    gy_res->Draw("ap");
    
    c_fit_res->cd(3);
    gz_res->SetTitle("Z residuals;x [cm];Δz [cm]");
    gz_res->SetMarkerStyle(2);
    gz_res->SetMarkerSize(0.4);
    gz_res->Draw("ap");
    
    c_fit_res->cd(4);
    gr_res->SetTitle("Residuals;x [cm];distance [cm]");
    gr_res->SetMarkerStyle(2);
    gr_res->SetMarkerSize(0.4);
    gr_res->Draw("ap");


    TCanvas* c_fit_res2 = new TCanvas("c_fit_res2", "Reconstruction residuals");
    c_fit_res2->Divide(2,2);
    c_fit_res2->cd(1); hx_res->Draw();
    c_fit_res2->cd(2); hy_res->Draw();
    c_fit_res2->cd(3); hz_res->Draw();
    c_fit_res2->cd(4); hr_res->Draw();

    double min = 0;//7
    double max = 15;//10.5

    // fitting both tracks with circles
    TF1* circle_fit  = FitCircle2(xz,min,max);
    TF1* circle_fit2 = FitCircle2(xz_reco,min,max);

    // loading magnetic field from txt file (units = meters)
    VectorField* magfield = new VectorField(-0.2,0.2,-0.3,0.3,-0.3,0.3,0.005);
    magfield->LoadField("/home/vavrik/work/X17/electron_positron_tracks/build/VecB.txt");
    
    double minfield,maxfield,minangle,maxangle;
    X17::GetMinMaxField(*magfield,minfield,maxfield);
    X17::GetMinMaxFieldAngle(*magfield,minangle,maxangle);
    cout << "At least 0.0 cm from TPC walls: minimal magnetic field: " << minfield << " maximal: " << maxfield << " minimal angle to electric field: " << minangle << " maximal: " << maxangle << "\n";
    X17::GetMinMaxField(*magfield,minfield,maxfield,0.5);
    X17::GetMinMaxFieldAngle(*magfield,minangle,maxangle,0.5);
    cout << "At least 0.5 cm from TPC walls: Minimal magnetic field: " << minfield << " maximal: " << maxfield << " minimal angle to electric field: " << minangle << " maximal: " << maxangle << "\n";
    X17::GetMinMaxField(*magfield,minfield,maxfield,1);
    X17::GetMinMaxFieldAngle(*magfield,minangle,maxangle,1);
    cout << "At least 1.0 cm from TPC walls: Minimal magnetic field: " << minfield << " maximal: " << maxfield << " minimal angle to electric field: " << minangle << " maximal: " << maxangle << "\n";
    
    TGraph* magnetic_x  = new TGraph();
    TGraph* magnetic_x2 = new TGraph();

    double step = 0.1;
    RecoEnergy(circle_fit,magfield,magnetic_x,min,max,step);
    RecoEnergy(circle_fit2,magfield,magnetic_x2,min,max,step);

    // Reconstruction with pads
    for (int i = 0; i < X17::channels; i++)
    {
        for (int j = 0; j < timebins; j++)
        {
            if(padhits[i][j] != 0)
            {
                double time = 100 * j + 50;
                double xpad,ypad;
                X17::GetPadCenter(i+1,xpad,ypad);

                SensorData reco = map->Invert(xpad,ypad,time);
                g_xyz_reco->AddPoint(reco.x1,reco.y1,reco.z1);
            }
        }        
    }

    TCanvas* c_track_xyz = new TCanvas("c_track_xyz","Electron track reconstruction with pads and time bins");
    
    // h_xyz_reco->Draw("box2");
    g_xyz_reco->SetTitle("Electron track reconstruction;x [cm]; y [cm];z [cm]");
    g_xyz_reco->SetMarkerStyle(2);
    g_xyz_reco->SetMarkerSize(0.4);
    g_xyz_reco->Draw("p");

    g_xyz->SetMarkerColor(2);
    g_xyz->SetMarkerStyle(7);
    g_xyz->SetMarkerSize(1.2);
    g_xyz->Draw("p same");

    TLegend* leg_xyz = new TLegend(0.129,0.786,0.360,0.887);
    leg_xyz->AddEntry(g_xyz,"original");
    leg_xyz->AddEntry(g_xyz_reco,"reconstructed");
    leg_xyz->Draw("same");

    vector<TPolyLine3D*> pad_lines;

    for (int i = 1; i <= X17::channels; i++)
    {
        double x1,y1,x2,y2;
        X17::GetPadCorners(i,x1,y1,x2,y2,true);

        TPolyLine3D* pad = new TPolyLine3D(5);
        pad->SetPoint(0,x1,y1,-2.5);
        pad->SetPoint(1,x1,y2,-2.5);
        pad->SetPoint(2,x2,y2,-2.5);
        pad->SetPoint(3,x2,y1,-2.5);
        pad->SetPoint(4,x1,y1,-2.5);

        pad_lines.push_back(pad);
    }
    
    for(auto l : pad_lines)   l->Draw("AL");
    
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

    // for (int i = 0; i < xz_reco->GetN(); i++)
    // {
    //     double x,y;
    //     xz_reco->GetPoint(i,x,y);
    //     residues->Fill(y-sp_fit->Eval(x));
    // }