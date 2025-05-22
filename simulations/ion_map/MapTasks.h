#pragma once

// ROOT dependencies
#include "TArrow.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH2F.h"
#include "TPolyLine3D.h"
#include "TStyle.h"
#include "TVectorD.h"

// X17 dependencies
#include "Field.h"
#include "MapTask.h"
#include "PadLayout.h"
#include "Points.h"
#include "Utilities.h"
#include "X17Utilities.h"

using namespace X17;

class Hist_YX_DX : public MapTask
{
    std::vector<TH2F*> v_yx_dx;

public:
    Hist_YX_DX(const X17::Field<X17::MapPoint>& map) : MapTask(map) { }

    void Z_Loop_Start(double z) override
    {
        std::string h_dx_name = "h_yx_dx_" + std::to_string(z);
        std::string h_dx_title = "XY plane dx (z = " + std::to_string(z) + ");x [cm];y [cm]";
        v_yx_dx.push_back(new TH2F(h_dx_name.c_str(),h_dx_title.c_str(),map.GetYCells(),map.GetYMin(),map.GetYMax(),map.GetXCells(),map.GetXMin(),map.GetXMax()));
    }

    void XYZ_Loop(double x, double y, double z, X17::MapPoint current) override
    {
        if(current.x() != 0) v_yx_dx.back()->Fill(y,x,current.x()-x);
    }

    void PostLoop() override
    {
        for(TH2F* h : v_yx_dx) h->Write();
    }
};

class Hist_YX_DY : public MapTask
{
    std::vector<TH2F*> v_yx_dy;

public:
    Hist_YX_DY(const X17::Field<X17::MapPoint>& map) : MapTask(map) { }

    void Z_Loop_Start(double z) override
    {
        std::string h_dy_name = "h_yx_dy_" + std::to_string(z);
        std::string h_dy_title = "XY plane dy (z = " + std::to_string(z) + ");x [cm];y [cm]";
        v_yx_dy.push_back(new TH2F(h_dy_name.c_str(),h_dy_title.c_str(),map.GetYCells(),map.GetYMin(),map.GetYMax(),map.GetXCells(),map.GetXMin(),map.GetXMax()));
    }

    void XYZ_Loop(double x, double y, double z, X17::MapPoint current) override
    {
        if(current.y() != 0) v_yx_dy.back()->Fill(y, x, current.y() - y);
    }

    void PostLoop() override
    {
        for(TH2F* h : v_yx_dy) h->Write();
    }
};


class Hist_YX_T1 : public MapTask
{
    double xmin,xmax,ymin,ymax;
    std::vector<TH2F*> v_yx_t1;

public:
    Hist_YX_T1(const X17::Field<X17::MapPoint>& map, double xmin, double xmax, double ymin, double ymax) : MapTask(map),xmin(xmin),xmax(xmax),ymin(ymin),ymax(ymax) { }

    void Z_Loop_Start(double z) override
    {
        std::string h_t1_name = "h_yx_t1_" + std::to_string(z);
        std::string h_t1_title = "XY plane t1 (z = " + std::to_string(z) + ");x [cm];y [cm]";
        v_yx_t1.push_back(new TH2F(h_t1_name.c_str(),h_t1_title.c_str(),map.GetYCells(),map.GetYMin(),map.GetYMax(),map.GetXCells(),map.GetXMin(),map.GetXMax()));
    }

    void XYZ_Loop(double x, double y, double z, X17::MapPoint current) override
    {
        if(current.t() != 0) v_yx_t1.back()->Fill(y,x,current.t());
    }

    void PostLoop() override
    {
        TCanvas* c_yx_t1 = new TCanvas("c_yx_t1","Map of drift times");
        c_yx_t1->Divide(4,4);
        for (int i = 0; i < v_yx_t1.size()-1; i++)
        {
            c_yx_t1->cd(i+1);
            v_yx_t1[i]->GetXaxis()->SetRangeUser(ymin,ymax);
            v_yx_t1[i]->GetYaxis()->SetRangeUser(xmin,xmax);
            v_yx_t1[i]->Draw("colz");
            double min_percent = 0.8 + 0.13 * (v_yx_t1.size() - 1 - i) / (v_yx_t1.size() - 1); // How high is minimum compared to maximum in percentage.
            v_yx_t1[i]->SetMinimum(min_percent * v_yx_t1[i]->GetMaximum());
            v_yx_t1[i]->SetMaximum(min_percent * v_yx_t1[i]->GetMaximum() + 500);
            X17::DrawTrapezoid(2);
        }
        c_yx_t1->Write();
    }
};

class Graph_YX : public MapTask
{
    double xmin,xmax,ymin,ymax;
    std::vector<TGraph*> v_g_yx;
    std::vector<std::vector<TArrow*>> vv_g_yx_arrows;
    std::vector<TArrow*> arrows;
    std::vector<std::vector<TEllipse*>> vv_g_yx_ellipses;
    std::vector<TEllipse*> ellipses;

    std::unordered_map<double, TGraph*>& v_g_yx_endpts;

public:
    Graph_YX(const X17::Field<X17::MapPoint>& map, double xmin, double xmax, double ymin, double ymax, std::unordered_map<double, TGraph*>& v_g_yx_endpts)
        : MapTask(map),xmin(xmin),xmax(xmax),ymin(ymin),ymax(ymax),v_g_yx_endpts(v_g_yx_endpts) { }

    void Z_Loop_Start(double z) override
    {
        std::string g_yx_name = "g_yx_" + std::to_string(z);
        std::string g_yx_title = "Map of electron readout positions, z = " + std::to_string(z);

        v_g_yx.push_back(new TGraph());

        v_g_yx.back()->SetName(g_yx_name.c_str());
        v_g_yx.back()->SetTitle(g_yx_title.c_str());
        v_g_yx.back()->SetMarkerColor(2);
        v_g_yx.back()->SetMarkerStyle(6);
        v_g_yx.back()->GetXaxis()->SetTitle("y' [cm]");
        v_g_yx.back()->GetXaxis()->SetTitleSize(0.043);
        v_g_yx.back()->GetXaxis()->SetLabelSize(0.04);
        v_g_yx.back()->GetYaxis()->SetTitle("x' [cm]");
        v_g_yx.back()->GetYaxis()->SetTitleSize(0.043);
        v_g_yx.back()->GetYaxis()->SetLabelSize(0.04);
    }

    void XYZ_Loop(double x, double y, double z, X17::MapPoint current) override
    {
        if((current.x() != 0) && (current.y() != 0))
        {
            v_g_yx.back()->AddPoint(current.y(),current.x());
            // v_g_yx.back()->SetPointError(v_g_yx.back()->GetN()-1,current.ydev(),current.xdev());
            TArrow* arrow = new TArrow(y,x,current.y(),current.x());
            arrows.push_back(arrow);
            ellipses.push_back(current.GetErrorEllipse(3,2.0,true));
        }
    }

    void Z_Loop_End() override
    {
        vv_g_yx_arrows.push_back(arrows);
        arrows.clear();
        vv_g_yx_ellipses.push_back(ellipses);
        ellipses.clear();
    }

    void PostLoop() override
    {        
        TCanvas* c_g_yx = new TCanvas("c_g_yx","Map of electron readout positions");
        c_g_yx->Divide(4,4);
        for (int i = 0; i < v_g_yx.size() - 1; i++)
        {
            int j = i / map.GetStep();
            if (j >= v_g_yx.size()) j = v_g_yx.size() - 1;

            c_g_yx->cd(i+1);
            v_g_yx[j]->GetXaxis()->SetRangeUser(ymin,ymax);
            v_g_yx[j]->GetYaxis()->SetRangeUser(xmin,xmax);
            v_g_yx[j]->Draw("AP");
            X17::DrawTrapezoid(2);

            for(TArrow* arr : vv_g_yx_arrows[j]) 
            {
                if(arr->GetX2() > ymin && arr->GetX2() < ymax && arr->GetY2() > xmin && arr->GetY2() < xmax)
                {
                    arr->SetArrowSize(0.002);
                    arr->Draw();
                }
            }

            for(TEllipse* ell : vv_g_yx_ellipses[j])
            {
                if(ell->GetX1() > ymin && ell->GetX1() < ymax && ell->GetY1() > xmin && ell->GetY1() < xmax)
                {
                    ell->SetFillStyle(0);
                    ell->Draw();
                }
            }

            double z = map.GetZMin() + j * map.GetStep();
            v_g_yx_endpts[z]->SetMarkerColor(kRed);
            v_g_yx_endpts[z]->Draw("P same");
        }
        c_g_yx->Write();

        double c_height = 500;
        double c_width = c_height * (ymax - ymin) / (xmax - xmin) * 8./9.;

        TCanvas* c_g_yx_all = new TCanvas("c_g_yx_all","All electron readout positions",c_width,c_height);
        ApplyThesisStyle(v_g_yx[0]);
        v_g_yx[0]->Draw("AP PMC");
        TGraph* bottom = v_g_yx[0];
        v_g_yx[0] = nullptr;
        for (TGraph* g : v_g_yx) if (g) g->Draw("P PMC");
        X17::DrawTrapezoid(2);
        ApplyThesisStyle(c_g_yx_all);
        c_g_yx_all->Write();

        TCanvas* c_g_yx_bottom = new TCanvas("c_g_yx_bottom","Bottom electron readout positions",c_width,c_height);
        bottom->Draw("AP");
        for(TArrow* arr : vv_g_yx_arrows[0]) 
        {
            if(arr->GetX2() > ymin && arr->GetX2() < ymax && arr->GetY2() > xmin && arr->GetY2() < xmax)
            {
                arr->SetArrowSize(0.002);
                arr->Draw();
            }
        }
        for(TEllipse* ell : vv_g_yx_ellipses[0])
        {
            if(ell->GetX1() > ymin && ell->GetX1() < ymax && ell->GetY1() > xmin && ell->GetY1() < xmax)
            {
                ell->SetFillStyle(0);
                ell->Draw();
            }
        }
        v_g_yx_endpts[-8]->Draw("P same");
        X17::DrawTrapezoid(2);
        ApplyThesisStyle(c_g_yx_bottom);
        c_g_yx_bottom->Write();
    }
};

class Graph_ZT : public MapTask
{
    TGraph* g_zt;

public:
    Graph_ZT(const X17::Field<X17::MapPoint>& map) : MapTask(map) { }

    void PreLoop() override
    {
        g_zt = new TGraph();

        g_zt->SetName("g_zt");
        g_zt->SetTitle("Initial height vs drift time");
        g_zt->SetMarkerColor(2);
        g_zt->SetMarkerStyle(6);
        g_zt->GetXaxis()->SetTitle("z [cm]");
        g_zt->GetYaxis()->SetTitle("t [ns]");
    }

    void XYZ_Loop(double x, double y, double z, X17::MapPoint current) override
    {
        if (current.t() != -1) g_zt->AddPoint(z,current.t());
    }

    void PostLoop() override
    {
        TCanvas* c_g_zt = new TCanvas("c_g_zt","Initial height vs drift time");
        ApplyThesisStyle(g_zt);
        g_zt->Draw("AP");
        ApplyThesisStyle(c_g_zt);
        c_g_zt->Write();
    }
};

class Graph_XZ : public MapTask
{
    TGraphErrors* g_xz;
    std::vector<TArrow*> v_g_xz_arrows;

public:
    Graph_XZ(const X17::Field<X17::MapPoint>& map) : MapTask(map) { }

    void PreLoop() override
    {
        g_xz = new TGraphErrors();
        g_xz->SetName("g_xz");
        g_xz->SetTitle("Map of electron displacement (y = 0) 2 sigma");
        g_xz->SetMarkerColor(2);
        g_xz->SetMarkerStyle(6);
        g_xz->GetXaxis()->SetTitle("x' [cm]");
        g_xz->GetYaxis()->SetTitle("z [cm]");
    }

    void ZX_Loop(double z, double x, X17::MapPoint current) override
    {
        if (current.t() != -1)
        {
            g_xz->AddPoint(current.x(),z);
            g_xz->SetPointError(g_xz->GetN() - 1, 2*current.xdev(), 0);
            TArrow* arrow = new TArrow(x,z,current.x(),z);
            v_g_xz_arrows.push_back(arrow);
        }
    }

    void PostLoop() override
    {
        TCanvas* c_g_xz = new TCanvas("c_g_xz","Map of ionization electron displacement");
        ApplyThesisStyle(g_xz);
        g_xz->Draw("AP");
        TLine* lleft  = new TLine(6.51,8,6.51,-8);
        TLine* lright = new TLine(14.61,8,14.61,-8);
        lleft->Draw();lright->Draw();
        // for(TArrow* arr : v_g_xz_arrows) {arr->SetArrowSize(0.004); arr->Draw();}
        ApplyThesisStyle(g_xz);
        ApplyThesisStyle(c_g_xz);
        c_g_xz->Write();
    }
};

class Graph_XT : public MapTask
{
    bool newdata;
    TGraph* g_xt;
    std::vector<TEllipse*> v_g_xt_ellipses;

public:
    Graph_XT(const X17::Field<X17::MapPoint>& map, bool newdata) : MapTask(map), newdata(newdata) { }

    void PreLoop() override
    {
        g_xt = new TGraphErrors();
        g_xt->SetName("g_xt");
        g_xt->SetTitle("Map of drift times (y = 0)");
        g_xt->SetMarkerColor(2);
        g_xt->SetMarkerStyle(2);
        g_xt->GetXaxis()->SetTitle("x' [cm]");
        g_xt->GetYaxis()->SetTitle("t [ns]");
    }

    void ZX_Loop(double z, double x, X17::MapPoint current) override
    {        
        g_xt->AddPoint(current.x(),current.t());
        v_g_xt_ellipses.push_back(current.GetErrorEllipse(1,2.0));
    }

    void PostLoop() override
    {
        TCanvas* c_g_xt = new TCanvas("c_g_xt","Map of ionization electron drift times");
        ApplyThesisStyle(g_xt);
        g_xt->Draw("AP");
        double ymax = newdata ? 18000 : 5000;
        TLine* lleft  = new TLine(6.51,0,6.51,ymax);
        TLine* lright = new TLine(14.61,0,14.61,ymax);
        lleft->Draw();
        lright->Draw();
        for (TEllipse* ell : v_g_xt_ellipses)
        {
            ell->SetFillStyle(0);
            ell->Draw();
        }
        c_g_xt->SetGrid();
        ApplyThesisStyle(c_g_xt);
        c_g_xt->Write();
    }
};

class Hist_XZ_T1 : public MapTask
{
    TH2F* xz_t1;

public:
    Hist_XZ_T1(const X17::Field<X17::MapPoint>& map) : MapTask(map) { }

    void PreLoop() override
    {
        xz_t1 = new TH2F("h_xz_t1","XZ plane t1;x [cm];z [cm]",map.GetXCells(),map.GetXMin(),map.GetXMax(),map.GetZCells(),map.GetZMin(),map.GetZMax());
    }

    void ZX_Loop(double z, double x, X17::MapPoint current) override
    {        
        xz_t1->Fill(x,z,current.t());
    }

    void PostLoop() override
    {
        xz_t1->Write();
    }
};

class GraphXYT : public MapTask
{
    TGraph2D* g_xyt;
    TGraph2D* g_endpts;
    std::vector<TPolyLine3D*> v_lines;

public:
    GraphXYT(const X17::Field<X17::MapPoint>& map, TGraph2D* g_endpts) : MapTask(map), g_endpts(g_endpts) { }

    void PreLoop() override
    {
        g_xyt = new TGraph2D();
    }
    
    void XYZ_Loop(double x, double y, double z, X17::MapPoint current) override
    {
        if (z >= -7.5) return;
        bool stop = (current.t() != -1) && (z < -7);
        bool no_z = true;
        
        g_xyt->AddPoint(current.x(),current.y(),current.t());
        current.Diagonalize(no_z);
        
        for (int i = 0; i < 3; i++)
        {
            X17::Vector err_vec = X17::Vector(current.eigen_vecs(0,i),current.eigen_vecs(1,i),current.eigen_vecs(3,i));
            err_vec *= StdevBiasFactor(current.n)*std::sqrt(current.eigen_vals(i));
            TPolyLine3D* line = new TPolyLine3D(2);
            line->SetPoint(0,current.x()-err_vec.x,current.y()-err_vec.y,current.t()-err_vec.z);
            line->SetPoint(1,current.x()+err_vec.x,current.y()+err_vec.y,current.t()+err_vec.z);
            v_lines.push_back(line);
        }
    }
    
    void PostLoop() override
    {
        TCanvas* c_g_xyt = new TCanvas("c_g_xyt","Map of ionization electron drift times");
        ApplyThesisStyle(g_xyt);
        g_xyt->Draw("AP");
        g_endpts->SetMarkerColor(kRed);
        g_endpts->SetTitle("XYT graph; x' [cm]; y' [cm]; t [ns]");
        ApplyThesisStyle(g_endpts);
        g_endpts->Draw("P");
        for(TPolyLine3D* line : v_lines) line->Draw();
        ApplyThesisStyle(c_g_xyt);
        c_g_xyt->Write();
    }
};