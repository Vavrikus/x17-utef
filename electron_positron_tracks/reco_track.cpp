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


struct VectorField
{
    double xmin,ymin,zmin,xmax,ymax,zmax;
    double step;
    int ximax,yimax,zimax;

    vector<vector<vector<Vector>>> field;

    VectorField(double xmin,double xmax,double ymin, double ymax, double zmin, double zmax, double step)
    : xmin(xmin),xmax(xmax),ymin(ymin),ymax(ymax),zmin(zmin),zmax(zmax),step(step)
    {        
        for (int i = 0; i <= floor((xmax-xmin)/step); i++)
        {
            vector<vector<Vector>> v1;
            for (int j = 0; j <= floor((ymax-ymin)/step); j++)
            {
                vector<Vector> v2;
                for (int k = 0; k <= floor((zmax-zmin)/step); k++)
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
        
        xi = floor((x-xmin)/step);
        yi = floor((y-ymin)/step);
        zi = floor((z-zmin)/step);
        // cout << xi << " " << yi << " " << zi << "\n";
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
        int lines_expected = floor((xmax-xmin)/step)*floor((ymax-ymin)/step)*floor((ymax-ymin)/step);

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
        
        double vx = (1-dx)*field[xi][yi][zi].vx+dx*field[xi2][yi][zi].vx;
        double vy = (1-dy)*field[xi][yi][zi].vy+dy*field[xi][yi2][zi].vy;
        double vz = (1-dz)*field[xi][yi][zi].vz+dz*field[xi][yi][zi2].vz;

        return Vector{vx,vy,vz};
    }
};

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

    TH2F* field = new TH2F("h_field","XZ magnetic field plot",121,-0.3,0.3,81,-0.2,0.2);
    VectorField* magfield = new VectorField(-0.3,0.3,-0.3,0.3,-0.2,0.2,0.005);
    magfield->LoadField("/home/vavrik/work/X17/electron_positron_tracks/build/VecB.txt");

    for(int i = 0; i < 121; i++)
        for(int j = 0; j < 81; j++)
        {
            Vector b = magfield->GetField(-0.3+i*0.005,0.0,-0.2+j*0.005);//magfield->field[i][60][j];//
            double bmag = sqrt(b.vx*b.vx+b.vy*b.vy+b.vz*b.vz);
            field->SetBinContent(i,j,bmag);
        }
    
    TCanvas* c4 = new TCanvas("c4","Magnetic field XZ plane");
    field->Draw("cont1");

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