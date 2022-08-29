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

using namespace std;

struct Vector
{
    double vx,vy,vz;
};

Vector operator+(const Vector& v1,const Vector& v2)
{
    return Vector{v1.vx+v2.vx,v1.vy+v2.vy,v1.vz+v2.vz};
}

Vector operator-(const Vector& v1,const Vector& v2)
{
    return Vector{v1.vx-v2.vx,v1.vy-v2.vy,v1.vz-v2.vz};
}

Vector operator*(const double& d,const Vector& v)
{
    return Vector{d*v.vx,d*v.vy,d*v.vz};
}

struct VectorField
{
    double xmin,ymin,zmin,xmax,ymax,zmax;
    double step;
    int ximax,yimax,zimax;

    vector<vector<vector<Vector>>> field;

    VectorField(double xmin,double xmax,double ymin, double ymax, double zmin, double zmax, double step)
    : xmin(xmin),xmax(xmax),ymin(ymin),ymax(ymax),zmin(zmin),zmax(zmax),step(step)
    {        
        for (int i = 0; i <= round((xmax-xmin)/step); i++)
        {
            vector<vector<Vector>> v1;
            for (int j = 0; j <= round((ymax-ymin)/step); j++)
            {
                vector<Vector> v2;
                for (int k = 0; k <= round((zmax-zmin)/step); k++)
                {
                    v2.push_back({0,0,0});
                }
                v1.push_back(v2);
            }
            field.push_back(v1);
        }

        this->GetVectorIndexes(xmax,ymax,zmax,ximax,yimax,zimax);
    }

    void GetVectorIndexes(double x, double y, double z, int& xi, int& yi, int& zi)
    {
        if(x<xmin||x>xmax||y<ymin||y>ymax||z<zmin||z>zmax)
            cerr << "Cannot read field out of bounds.\n";
        
        xi = round((x-xmin)/step);
        yi = round((y-ymin)/step);
        zi = round((z-zmin)/step);
    }

    Vector* GetVector(double x, double y, double z)
    {
        int xi,yi,zi;
        this->GetVectorIndexes(x,y,z,xi,yi,zi);

        return &(field[xi][yi][zi]);
    }

    void LoadField(const char* filename)
    {
        std::ifstream inf {filename};

        int lines_read = 0;
        int lines_processed = 0;
        int lines_expected = round((xmax-xmin)/step+1)*round((ymax-ymin)/step+1)*round((zmax-zmin)/step+1);

        while (inf)
        {
            std::string X,Y,Z,VX,VY,VZ;
            inf >> X; inf >> Y; inf >> Z; inf >> VX; inf >> VY; inf >> VZ;
            lines_read++;
            
            if (X != "")
            {
                double x,y,z,vx,vy,vz;
                x = stod(X); y = stod(Y); z = stod(Z);
                vx = stod(VX); vy = stod(VY); vz = stod(VZ);

                *(this->GetVector(x,y,z)) = Vector{vx,vy,vz};
                lines_processed++;
            }
        }

        cout << "Lines read: " << lines_read << " processed: " << lines_processed << " expected: " << lines_expected << "\n";
    }

    Vector GetField(double x, double y, double z)
    {
        int xi,yi,zi;
        this->GetVectorIndexes(x,y,z,xi,yi,zi);
        
        int xi2,yi2,zi2;
        if(x-(xmin+step*xi)<0) xi2 = xi-1; else xi2 = xi + 1; if(xi2>ximax) xi2=xi;
        if(y-(ymin+step*yi)<0) yi2 = yi-1; else yi2 = yi + 1; if(yi2>yimax) yi2=yi;
        if(z-(zmin+step*zi)<0) zi2 = zi-1; else zi2 = zi + 1; if(zi2>zimax) zi2=zi;

        double dx = abs((x-(xmin+step*xi))/step);
        double dy = abs((y-(ymin+step*yi))/step);
        double dz = abs((z-(zmin+step*zi))/step);

        Vector& c000 = field[xi][yi][zi];
        Vector& c001 = field[xi][yi][zi2];
        Vector& c010 = field[xi][yi2][zi];
        Vector& c011 = field[xi][yi2][zi2];
        Vector& c100 = field[xi2][yi][zi];
        Vector& c101 = field[xi2][yi][zi2];
        Vector& c110 = field[xi2][yi2][zi];
        Vector& c111 = field[xi2][yi2][zi2];

        Vector c00 = (1-dx)*c000+dx*c100;
        Vector c01 = (1-dx)*c001+dx*c101;
        Vector c10 = (1-dx)*c010+dx*c110;
        Vector c11 = (1-dx)*c011+dx*c111;

        Vector c0 = (1-dy)*c00+dy*c10;
        Vector c1 = (1-dy)*c01+dy*c11;

        return (1-dz)*c0+dz*c1;
    }
};

template<int nodes>
TSpline3* FitSplines(TGraph* graph, const double min, const double max)
{
    double node_pos[nodes];
    for (int i = 0; i < nodes; i++) node_pos[i] = min+(i/(nodes-1.0))*(max-min);

    NSpline<nodes>* s = new NSpline<nodes>(node_pos);

    TF1* fit = new TF1("zz_fit",s->GetEval(),min,max,2*nodes+2);
    graph->Fit(fit,"M","",min,max);

    double* fit_pars = fit->GetParameters();
    double nodes_y[nodes];
    for(int i = nodes; i < 2*nodes; i++) nodes_y[i-nodes]=fit_pars[i];
    double beg1 = fit_pars[2*nodes];
    double end1 = fit_pars[2*nodes+1];
    
    return new TSpline3("sp3",node_pos,nodes_y,nodes,"b1e1",beg1,end1);
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

int reco_track()
{
    //file with Garfield simulation output
    TFile* inFile = new TFile("build/electrons.root");
    TTree* electrons = (TTree*)inFile->Get("electrons");

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

    //zy (track) plot from original track and reconstructed from drif time
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

    const int nodes = 25;
    const int nodes2 = 5;
    double min = 0;
    double max = 15;

    //fitting both tracks with splines
    TSpline3* sp_fit  = FitSplines<nodes>(zy,min,max);
    TSpline3* sp_fit2 = FitSplines<nodes2>(zy_reco,min,max);

    //loading magnetic field from txt file (units = meters)
    VectorField* magfield = new VectorField(-0.3,0.3,-0.3,0.3,-0.2,0.2,0.005);
    magfield->LoadField("/home/vavrik/work/X17/electron_positron_tracks/build/VecB.txt");
    
    TGraph* energy_x    = new TGraph();
    TGraph* radius_x    = new TGraph();
    TGraph* magnetic_x  = new TGraph();
    TGraph* energy_x2   = new TGraph();
    TGraph* radius_x2   = new TGraph();
    TGraph* magnetic_x2 = new TGraph();

    double step = 0.1;
    RecoEnergy(sp_fit,magfield,energy_x,radius_x,magnetic_x,min,max,step);
    RecoEnergy(sp_fit2,magfield,energy_x2,radius_x2,magnetic_x2,min,max,step);

    TCanvas* c_energy = new TCanvas("c_energy","Reconstructed energy");
    energy_x->SetTitle("Reconstructed energy;z [cm]; Energy [MeV]");
    energy_x->SetMarkerStyle(2);
    energy_x->SetMarkerSize(0.4);
    energy_x->SetMarkerColor(2);
    energy_x->Draw("ap");

    energy_x2->SetMarkerStyle(2);
    energy_x2->SetMarkerSize(0.4);
    energy_x2->Draw("p same");

    TCanvas* c_radius = new TCanvas("c_radius","Reconstructed radius");
    radius_x->SetTitle("Reconstructed radius;z [cm]; Radius [cm]");
    radius_x->SetMarkerStyle(2);
    radius_x->SetMarkerSize(0.4);
    radius_x->SetMarkerColor(2);
    radius_x->Draw("ap");

    radius_x2->SetMarkerStyle(2);
    radius_x2->SetMarkerSize(0.4);
    radius_x2->Draw("p same");
    
    TCanvas* c_magnetic = new TCanvas("c_magnetic","Perpendicular magnetic field");
    magnetic_x->SetTitle("Perpendicular magnetic field;z [cm]; B [T]");
    magnetic_x->SetMarkerStyle(2);
    magnetic_x->SetMarkerSize(0.4);
    magnetic_x->SetMarkerColor(2);
    magnetic_x->Draw("ap");

    magnetic_x2->SetMarkerStyle(2);
    magnetic_x2->SetMarkerSize(0.4);
    magnetic_x2->Draw("p same");

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