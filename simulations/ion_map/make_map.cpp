// C++ dependencies
#include <iostream>
#include <string>

// ROOT dependencies
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TStyle.h"

// Dictionaries generated by ROOT
#include "../build/X17_dict.cxx"

// X17 dependencies
// #define DEBUG
#include "Color.h"
#include "Field.h"
#include "MapTask.h"
#include "PadLayout.h"
#include "Points.h"
#include "Utilities.h"
#include "X17Utilities.h"

#include "MapTasks.h"

/// @brief Loads ionization electron data from files (named ion(id).root) in a specified folder.
/// @param max_id The maximum ID number of the files to load.
/// @param folder The folder containing the files to load. Default is "../../data/ion_map/sample_2.0/".
/// @return A TChain pointer containing the loaded data.
TChain* LoadData(int max_id, std::string folder = "../../../data/ion_map/sample_2.0/")
{
    TChain* map_data = new TChain("map_data","Data from ionization electrons simulation.");

    for (int i = 1; i <= max_id; i++)
    {
        std::string filepath = folder + "ion" + std::to_string(i) + ".root?#electrons";
        map_data->Add(filepath.c_str());
    }

    return map_data;
}

TGraph* GetQQconfidenceGraph(double sigma, std::vector<std::vector<double>> sq_mahal_sets, TGraph* g_rand, bool no_z = true)
{
    int nsets = sq_mahal_sets.size();
    int rand_size = sq_mahal_sets.back().size();
    int degrees_of_freedom = no_z ? 3 : 4;

    TGraph* g_rand_ci = new TGraph(2*rand_size);
    double confidence = TMath::Erf(sigma/std::sqrt(2));
    for(int i = 0; i < rand_size; i++)
    {
        double chi2_quantile = TMath::ChisquareQuantile((i+1.0)/(rand_size+1.0),degrees_of_freedom);

        std::vector<double> sq_mahal;
        sq_mahal.reserve(nsets);
        for (int j = 0; j < nsets; j++)
        {
            sq_mahal.push_back(sq_mahal_sets[j][i]);
        }
        g_rand->AddPoint(GetQuantile(sq_mahal,0.95),chi2_quantile);
        g_rand_ci->SetPoint(rand_size-i-1,chi2_quantile,GetQuantile(sq_mahal,0.5+confidence/2,true));
        g_rand_ci->SetPoint(rand_size+i,chi2_quantile,GetQuantile(sq_mahal,0.5-confidence/2,true));
    }

    return g_rand_ci;
}

int make_map()
{
    // Load data from all files (results of individual jobs).
    const bool new_data = true;

    std::string data_folder = new_data ? "../../../data/ion_map/sample_2.0/" : "../../../data/ion_map/sample_1.0/";
    int files_count         = new_data ? 1000 : 200;
    double step_size        = new_data ?  0.5 : 1.0;
    double map_xmin         = new_data ? -1.5 : 0.0;

    TChain* map_data_in = LoadData(files_count,data_folder);
    std::cout << "Number of simulated electrons: " << map_data_in->GetEntries() << "\n";

    // Set branches for TChain containing data.
    X17::MicroPoint point; // An object that will hold the information about current ionization electron loaded.
    point.SetTChainBranches(map_data_in,!new_data);

    // Prepare the field that will hold the final values.
    X17::Field<X17::MapPoint> map({map_xmin,-30,-8},{15,30,8},step_size,X17::MapPoint());

    // Variables for checking with the previous position.
    X17::Vector v_prev(0,0,0);        // Vector for the previous position. Used for comparison with current position.
    int same_prev = 1;                // The number of entries with the same position including the current entry.
    X17::EndPoint p_avg(0,0,0,0);     // The object used for calculating of the averages.
    std::vector<X17::EndPoint> p_vec; // Stores the endpoints with the same initial position.

    // Simulation endpoint TGraph
    TGraph2D* g_endpts = new TGraph2D();
    TGraph2D* g_sim  = new TGraph2D();
    TGraph2D* g_reco = new TGraph2D();

    std::unordered_map<double, TGraph*> v_g_yx_endpts;

    TRandom3* rand = new TRandom3(0);
    TGraph* g_qq = new TGraph();
    TGraph* g_rand = new TGraph();

    int nsets = 10000;
    std::vector<std::vector<double>> sq_mahal_sets;
    sq_mahal_sets.reserve(nsets);

    std::vector<double> mardiaB;
    double Bmin = -7.5;
    double Bmax = 7.5;
    // Scott's rule
    double Bbins      = GetBinsScott(Bmin,Bmax,1.0,651*16);
    double Bbins_rand = GetBinsScott(Bmin,Bmax,1.0,nsets);
    TH1F* h_mardia       = new TH1F("h_mardia","Mardia's multivariate normality test",Bbins,Bmin,Bmax);
    TH1F* h_mardia_rand  = new TH1F("h_mardia_rand","Mardia's multivariate normality test",Bbins_rand,Bmin,Bmax);
    
    std::vector<double> mardiaA;
    double Amin = 0;
    double Amax = 35;
    // Fridman-Diaconis rule
    double Abins      = (Amax-Amin)/(2*5.8/std::pow(651*16,1.0/3.0));
    double Abins_rand = (Amax-Amin)/(2*5.8/std::pow(nsets,1.0/3.0));
    TH1F* h_mardia2      = new TH1F("h_mardia2","Mardia's multivariate normality test",Abins,Amin,Amax);
    TH1F* h_mardia_rand2 = new TH1F("h_mardia_rand2","Mardia's multivariate normality test",Abins_rand,Amin,Amax);

    TH1F* h_pvalsA = new TH1F("h_pvalsA","p-values",std::pow(651*16,1.0/3.0),0,1);
    TH1F* h_pvalsB = new TH1F("h_pvalsB","p-values",std::pow(651*16,1.0/3.0),0,1);

    bool filled_mardia = false;
    std::cout << "Calculating the map...\n";

    // Calculate the averages and standard deviations.
    for (int i = 0; i < map_data_in->GetEntries(); i++)
    {
        map_data_in->GetEntry(i);
        
        ReportProgress(i+1,map_data_in->GetEntries());

        if (v_g_yx_endpts.find(point.start.z()) == v_g_yx_endpts.end())
            v_g_yx_endpts[point.start.z()] = new TGraph();
        v_g_yx_endpts[point.start.z()]->AddPoint(point.end.y(),point.end.x());

        if (point.start.z() == -8)
            g_endpts->AddPoint(point.end.x(),point.end.y(),point.end.t);

        if (point.GetInitPos() == v_prev && (i != map_data_in->GetEntries()-1))
        {
            p_avg += point.end;
            p_vec.push_back(point.end);
            same_prev++;
        }

        else
        {
            if (i == map_data_in->GetEntries() - 1) 
            {
                same_prev++;
                p_vec.push_back(point.end);
                p_avg += point.end;
            }

            p_avg /= same_prev;

            if (i != 0) 
            {
                map.SetPoint(v_prev,X17::MapPoint(p_vec.size(),p_avg,CovarianceMatrix(p_vec,p_avg)));
                // std::cout << v_prev.x << "\t" << v_prev.y << "\t" << v_prev.z << "\n";
                g_sim->AddPoint(v_prev.x,v_prev.y,v_prev.z);
                g_reco->AddPoint(p_avg.x(),p_avg.y(),p_avg.t);
                
                if (!filled_mardia && v_prev == X17::Vector(8,3,-8))//15,9,-8))//14.5,7,-8))//
                {
                    std::cout << "\nCalculating Mardia's test Monte Carlo...\n";

                    FillQQplot(g_qq,p_vec);
                    
                    X17::MapPoint current = *map.GetPoint(v_prev);
                    current.Diagonalize(true);
                    
                    for (int n = 0; n < nsets; n++)
                    {
                        ReportProgress(n+1,nsets);

                        std::vector<EndPoint> p_rand;
                        p_rand.reserve(p_vec.size());
                        
                        for (int j = 0; j < p_vec.size(); j++)
                        {
                            X17::EndPoint rand_point = current.GetRandomPoint(rand);
                            // g_endpts->AddPoint(rand_point.x(),rand_point.y(),rand_point.t);
                            p_rand.push_back(rand_point);
                        }                     
                        
                        // FillQQplot(g_rand,p_rand,p_avg);
                        sq_mahal_sets.push_back(SqMahalanobis(p_rand));//CovarianceMatrix(p_rand,p_avg)));//
                        std::sort(sq_mahal_sets.back().begin(),sq_mahal_sets.back().end());

                        double A,B;
                        MardiaTest(p_rand,A,B);
                        h_mardia_rand->Fill(B);
                        h_mardia_rand2->Fill(A);
                        mardiaA.push_back(A);
                        mardiaB.push_back(B);
                    }
                    
                    filled_mardia = true;
                    i = 0;

                    std::cout << "\nMonte Carlo done!\n";
                    continue;
                }
                
                if (filled_mardia)
                {
                    // std::cout << "X: " << v_prev.x << "\tY: " << v_prev.y << "\n"; 
                    X17::MapPoint* current = map.GetPoint(v_prev);
                    
                    double A,B;
                    MardiaTest(p_vec,A,B);

                    double pA = GetPvalue(mardiaA,A);
                    double pB = GetPvalue(mardiaB,B);

                    current->mardia_A = pA;
                    current->mardia_B = pB;
                    
                    if (v_prev.z != 8)
                    {
                        h_mardia->Fill(B);
                        h_mardia2->Fill(A);
                        h_pvalsA->Fill(pA);
                        h_pvalsB->Fill(pB);
                    }
                }
                
                if (p_vec.size() < 100) std::cout << "WARNING: Map point has only " << p_vec.size() << " entries!\n";
                p_vec.clear();
            }

            p_vec.push_back(point.end);
            p_avg  = point.end;
            v_prev = point.GetInitPos();
            same_prev = 1;
        }
    }
    
    // Save the compiled map.
    data_folder += "map.root";
    TFile* outfile = new TFile(data_folder.c_str(),"RECREATE");
    outfile->WriteObjectAny(&map,"X17::Field<X17::MapPoint>","map");
    const X17::Field<X17::MapPoint>& cmap = map; // Map should no longer be modified.
    
    // Plotting.
    bool MakePlots = true;
    
    if(MakePlots)
    {
        // Plots from the map creation.
            TCanvas* c_sim = new TCanvas("c_sim","Simulation");
            g_sim->SetMarkerStyle(6);
            g_sim->SetTitle(";x [cm];y [cm];z [cm]");
            g_sim->GetXaxis()->SetTitleOffset(1.65);
            g_sim->GetYaxis()->SetTitleOffset(1.65);
            g_sim->GetXaxis()->SetTitleSize(0.043);
            g_sim->GetYaxis()->SetTitleSize(0.043);
            g_sim->GetZaxis()->SetTitleSize(0.043);
            g_sim->Draw("P");
            c_sim->Write();
        
            TCanvas* c_reco = new TCanvas("c_reco","Reconstruction");
            g_reco->SetMarkerStyle(6);
            g_reco->SetTitle(";x' [cm];y' [cm];t [ns]");
            g_reco->GetXaxis()->SetTitleOffset(1.65);
            g_reco->GetYaxis()->SetTitleOffset(1.65);
            g_reco->GetZaxis()->SetTitleOffset(1.2);
            g_reco->GetXaxis()->SetTitleSize(0.043);
            g_reco->GetYaxis()->SetTitleSize(0.043);
            g_reco->GetZaxis()->SetTitleSize(0.043);
            g_reco->Draw("P");
            c_reco->Write();
        
            TCanvas* c_qq = new TCanvas("c_qq","QQ plot");
            
            int rand_size = sq_mahal_sets.back().size();
            TGraph* g_rand_ci = GetQQconfidenceGraph(3.0,sq_mahal_sets,g_rand);
            g_rand_ci->SetFillColor(X17::Color::RGB(255,0,0,15));
            g_rand_ci->SetFillStyle(1001);
            g_rand_ci->GetXaxis()->SetTitle("#chi^{2}_{3} quantile");
            g_rand_ci->GetYaxis()->SetTitle("Squared Mahalanobis distance");
            g_rand_ci->SetTitle("");
            g_rand_ci->Draw("AF");
        
            TGraph* g_rand_ci2 = GetQQconfidenceGraph(2.0,sq_mahal_sets,g_rand);
            g_rand_ci2->SetFillColor(X17::Color::RGB(255,0,0,35));
            g_rand_ci2->SetFillStyle(1001);
            g_rand_ci2->Draw("F same");
        
            g_rand->SetMarkerStyle(2);
            g_rand->SetMarkerColor(kRed-9);
            // g_rand->Draw("AP");
        
            g_qq->SetMarkerStyle(4);
            g_qq->Draw("P same");
        
            TF1* line = new TF1("line","x",0,20);
            line->SetLineColor(kRed);
            line->Draw("same");
        
            g_rand_ci->SetLineColor(kRed-4);
            g_rand_ci->Draw("L same");
            c_qq->Write();
            // c_qq->SaveAs("qq.png");
        
            TCanvas* c_mardia = new TCanvas("c_mardia","Mardia test");
            h_mardia->Scale(1./(h_mardia->Integral("width")));
            h_mardia->SetLineColor(kBlack);
            h_mardia->SetLineWidth(2);
            h_mardia->Draw("hist");
            h_mardia_rand->Scale(1./h_mardia_rand->Integral("width"));
            h_mardia_rand->SetLineColor(kRed);
            h_mardia_rand->Draw("hist same");
            TF1* f_gaus = new TF1("gaus", "1/sqrt(2*TMath::Pi()) * exp(-0.5 * x*x)", -5, 5);
            f_gaus->SetLineColor(kBlue);
            f_gaus->Draw("same");
            TLegend* legend = new TLegend(0.7,0.7,0.9,0.9);
            legend->AddEntry(h_mardia,"Data","l");
            legend->AddEntry(h_mardia_rand,"Monte Carlo","l");
            legend->AddEntry(f_gaus,"Asymptotic","l");
            legend->Draw("same");
            c_mardia->Write();
        
            TCanvas* c_mardia2 = new TCanvas("c_mardia2","Mardia test 2");
            h_mardia2->Scale(1./h_mardia2->Integral("width"));
            h_mardia2->SetLineColor(kBlack);
            h_mardia2->SetLineWidth(2);
            h_mardia2->Draw("hist same");
            h_mardia_rand2->Scale(1./h_mardia_rand2->Integral("width"));
            h_mardia_rand2->SetLineColor(kRed);
            h_mardia_rand2->Draw("hist same");
            int k = 3*4*5/6; // Degrees of freedom
            TF1 *chi2_pdf = new TF1("chi2_pdf", "[0] * TMath::Power(x, [1]/2 - 1) * TMath::Exp(-x/2) / (TMath::Gamma([1]/2) * TMath::Power(2, [1]/2))", 0, 30);
            chi2_pdf->SetParameters(1.0, k); // [0] = normalization, [1] = k (dof)
            chi2_pdf->SetLineColor(kBlue);
            chi2_pdf->Draw("same");
            TLegend* legend2 = new TLegend(0.7,0.7,0.9,0.9);
            legend2->AddEntry(h_mardia2,"Data","l");
            legend2->AddEntry(h_mardia_rand2,"Monte Carlo","l");
            legend2->AddEntry(chi2_pdf,"Asymptotic","l");
            legend2->Draw("same");
            c_mardia2->Write();

            TCanvas* c_pvals = new TCanvas("c_pvals","P-values");
            h_pvalsA->Scale(1./h_pvalsA->Integral("width"));
            h_pvalsA->SetLineColor(kRed);
            h_pvalsA->Draw("hist");
            h_pvalsB->Scale(1./h_pvalsB->Integral("width"));
            h_pvalsB->SetLineColor(kBlue);
            h_pvalsB->Draw("hist same");
            TLegend* legend3 = new TLegend(0.7,0.7,0.9,0.9);
            legend3->AddEntry(h_pvalsA,"A statistic","l");
            legend3->AddEntry(h_pvalsB,"B statistic","l");
            legend3->Draw("same");
            h_pvalsB->GetYaxis()->SetRangeUser(0,1.2); 
            c_pvals->Write();
            

        // Plotting limits for some of the plots.
        double ymin = -10; // Minimal plotted y-coordinate.
        double ymax =  10; // Maximal plotted y-coordinate.
        double xmin =   5; // Minimal plotted x-coordinate.
        double xmax =  16; // Maximal plotted x-coordinate.
        
        std::vector<MapTask*> plot_tasks; // The vector containing all plotting tasks to be plotted.
        // plot_tasks.push_back(new Hist_YX_DX(cmap));
        // plot_tasks.push_back(new Hist_YX_DY(cmap));
        // plot_tasks.push_back(new Hist_YX_T1(cmap,xmin,xmax,ymin,ymax));
        plot_tasks.push_back(new Graph_YX(cmap,xmin,xmax,ymin,ymax,v_g_yx_endpts));
        plot_tasks.push_back(new Graph_ZT(cmap));
        plot_tasks.push_back(new Graph_XZ(cmap));
        plot_tasks.push_back(new Graph_XT(cmap,new_data));
        // plot_tasks.push_back(new Hist_XZ_T1(cmap));
        plot_tasks.push_back(new GraphXYT(cmap,g_endpts));

        for(MapTask* t : plot_tasks) t->PreLoop();

        for (int z_xyi = 0; z_xyi < map.GetZCells(); z_xyi++)
        {
            double z = map.GetZMin() + z_xyi * map.GetStep();
            for(MapTask* t : plot_tasks) t->Z_Loop_Start(z);

            for (int xi = 0; xi < map.GetXCells(); xi++)
            {
                double x = map.GetXMin()+xi*map.GetStep();

                for (int yi = 0; yi < map.GetYCells(); yi++)
                {
                    double y = map.GetYMin() + yi * map.GetStep();
                    X17::MapPoint& current = map.at(xi,yi,z_xyi);
                    if (current.t() == -1) continue;

                    for(MapTask* t : plot_tasks) t->XYZ_Loop(x,y,z,current);
                }
            }
            for(MapTask* t : plot_tasks) t->Z_Loop_End();
        }

        double y_xz = 0;                                       // The y-coordinate of xz plots.
        int y_xzi   = round(y_xz-map.GetYMin())/map.GetStep(); // The y-coordinate index of xz plots.

        for (int xi = 0; xi < map.GetXCells(); xi++)
        {
            double x = map.GetXMin() + xi * map.GetStep();
            for (int zi = 0; zi < map.GetZCells(); zi++)
            {
                double z = map.GetZMin() + zi * map.GetStep();
                X17::MapPoint& current = map.at(xi,y_xzi,zi);

                for(MapTask* t : plot_tasks) t->ZX_Loop(z,x,current);
            }
        }

        for(MapTask* t : plot_tasks) t->PostLoop();
        
        // Drawing the distortion of the pads.
        TCanvas* c_pads = new TCanvas("c_pads", "Pads distortion for different times.");
        c_pads->Divide(4,4);
        
        for (int i = 0; i < 16; i++)
        {
            using namespace X17::constants;
            
            c_pads->cd(i+1);
            TGraph* gr = new TGraph();
            gr->AddPoint(-yhigh,xmin);
            gr->AddPoint(yhigh,xmax);
            
            std::string gr_title = "Pads for t = " + std::to_string(i+1) + " #mus";
            gr->SetTitle(gr_title.c_str());
            // gr->GetHistogram()->SetTitleSize(1);
            gr->GetXaxis()->SetLabelSize(0.05);
            gr->GetYaxis()->SetLabelSize(0.05);
            gr->GetXaxis()->SetTitleSize(0.07);
            gr->GetYaxis()->SetTitleSize(0.07);
            gr->GetXaxis()->SetTitleOffset(0.6);
            gr->GetYaxis()->SetTitleOffset(0.6);
            gr->GetXaxis()->SetTitle("x [cm]");
            gr->GetYaxis()->SetTitle("y [cm]");
            gr->Draw("AP");
            
            X17::DefaultLayout& pads = X17::DefaultLayout::GetDefaultLayout();
            pads.DrawPadsDistortion((i + 1) * 16000 / 16, c_pads, &cmap);
            // X17::DrawTrapezoid();
        }
        
        c_pads->Write();
        outfile->Close();
    }

    return 0;
}

int main(int argc, char const *argv[])
{
    return make_map();
}