#pragma once

// ROOT dependencies
#include "TStyle.h"

// X17 dependencies
#include "PadLayout.h"
#include "Reconstruction.h"
#include "TrackLoop.h"
#include "X17Utilities.h"

/// @brief Class for padded reconstruction.
class RecoPadsTask : public X17::RecoTask
{
    friend class MicroCircleAndRKFitTask;
    friend class MicroFitAndSaveTask;

    static constexpr int timebins = 200;
    int padhits[X17::constants::channels][timebins];
    std::vector<X17::RecoPoint> reco_data;
    TGraph2D *g_xyz, *g_xyz_reco;
    TCanvas* c_reco;

    double height; // Height at which the pads will be drawn.

    void PreElectronLoop() override
    {
        for (int i = 0; i < X17::constants::channels; i++)
        for (int j = 0; j < timebins; j++)
            padhits[i][j] = 0;

        if (m_loop->make_track_plots)
        {
            g_xyz      = new TGraph2D();
            g_xyz_reco = new TGraph2D();
        }
    }

    void ElectronLoop() override
    {
        X17::MicroPoint micro = m_loop->curr_micro;

        int channel = X17::DefaultLayout::GetDefaultLayout().GetPad(micro.x1(),micro.y1());
        int timebin = micro.t1() / 100.0;

        if (timebin > timebins - 1) std::cerr << "ERROR: Invalid timebin: " << timebin << std::endl;
        if (channel == -1) std::cerr << "ERROR: No pad hit found. Coordinates x,y: " << micro.x1() << ", " << micro.y1() << std::endl;

        padhits[channel-1][timebin]++;
        
        if (m_loop->make_track_plots) g_xyz->AddPoint(micro.x0(),micro.y0(),micro.z0());
    }

    void PostElectronLoop() override
    {
        using namespace X17::constants;

        // Reconstruction with pads.
        int max_charge = 0;
        X17::Vector vmin = {100,100,100};
        X17::Vector vmax = {-100,-100,-100};
        std::vector<TPolyLine3D*> marker_lines;

        for (int i = 0; i < channels; i++)
            for (int j = 0; j < timebins; j++)
                if (padhits[i][j] > max_charge) max_charge = padhits[i][j];

        std::cout << "Max charge: " << max_charge << std::endl;

        for (int i = 0; i < channels; i++)
        {
            // ReportProgress(i,channels);
            for (int j = 0; j < timebins; j++)
            {
                int charge = padhits[i][j];
                if (charge == 0) continue;

                double time = 100 * j + 50;
                double xpad,ypad;
                X17::DefaultLayout::GetDefaultLayout().GetPadCenter(i+1,xpad,ypad);

                X17::RecoPoint reco = Reconstruct(m_loop->map,X17::EndPoint(xpad,ypad,zmax,time));
                reco_data.emplace_back(reco.x(),reco.y(),reco.z(),charge);
                if (m_loop->make_track_plots)
                    g_xyz_reco->AddPoint(reco.x(),reco.y(),reco.z());

                // Drawing markers.
                bool old_reco = false;
                const int pts = old_reco ? 4 : 100;
                double norm_charge = (double) charge / (double) max_charge;
                int color_index = norm_charge * (gStyle->GetNumberOfColors() - 1);
                int color = gStyle->GetColorPalette(color_index);
                
                double x1,y1,x2,y2;
                X17::DefaultLayout::GetDefaultLayout().GetPadCorners(i+1,x1,y1,x2,y2);
                if (vmin.x > x1) vmin.x = x1;
                if (vmin.y > y1) vmin.y = y1;
                if (vmax.x < x2) vmax.x = x2;
                if (vmax.y < y2) vmax.y = y2;
                
                X17::Vector corners[4] = {
                    {x1, y1, -8},
                    {x1, y2, -8},
                    {x2, y1, -8},
                    {x2, y2, -8}
                };

                // Corner lines
                for (auto v_corner : corners)
                {
                    TPolyLine3D* corner_line = new TPolyLine3D(pts);
                    corner_line->SetLineColor(color);

                    for (int it = 0; it < pts; it++)
                    {
                        double t = time - 50 + 100.0 * it / (pts - 1.0);
                        X17::EndPoint t_corner (v_corner, t);
                        X17::RecoPoint reco2 = old_reco ? X17::ReconstructOld(m_loop->map,t_corner,1e-6,false) : X17::Reconstruct(m_loop->map,t_corner);
                        X17::Vector line_point = reco.AsVector() + norm_charge * (reco2.AsVector() - reco.AsVector());
                        corner_line->SetPoint(it,line_point.x,line_point.y,line_point.z);

                        if (line_point.z < vmin.z) vmin.z = line_point.z;
                        if (line_point.z > vmax.z) vmax.z = line_point.z;
                    }

                    marker_lines.push_back(corner_line);
                }

                // Plane lines
                for (double t : {time - 50, time + 50})
                {
                    TPolyLine3D* lines[4] = {
                        new TPolyLine3D(pts),
                        new TPolyLine3D(pts),
                        new TPolyLine3D(pts),
                        new TPolyLine3D(pts)
                    };

                    for (auto line : lines)
                    {
                        line->SetLineColor(color);
                        marker_lines.push_back(line);
                    }

                    for (int k = 0; k < pts; k++)
                    {
                        std::vector<X17::RecoPoint> points;
                        if (old_reco)
                        {
                            points = {
                                X17::ReconstructOld(m_loop->map,{x1 + k*(x2 - x1) / (pts - 1), y1, -8, t},1e-6,false),
                                X17::ReconstructOld(m_loop->map,{x2, y1 + k*(y2 - y1) / (pts - 1), -8, t},1e-6,false),
                                X17::ReconstructOld(m_loop->map,{x1 + k*(x2 - x1) / (pts - 1), y2, -8, t},1e-6,false),
                                X17::ReconstructOld(m_loop->map,{x1, y1 + k*(y2 - y1) / (pts - 1), -8, t},1e-6,false)
                            };
                        }

                        else
                        {
                            points = {
                                X17::Reconstruct(m_loop->map,{x1 + k*(x2 - x1) / (pts - 1), y1, -8, t}),
                                X17::Reconstruct(m_loop->map,{x2, y1 + k*(y2 - y1) / (pts - 1), -8, t}),
                                X17::Reconstruct(m_loop->map,{x1 + k*(x2 - x1) / (pts - 1), y2, -8, t}),
                                X17::Reconstruct(m_loop->map,{x1, y1 + k*(y2 - y1) / (pts - 1), -8, t})
                            };
                        }

                        for (int l = 0; l < 4; l++)
                        {
                            X17::Vector line_point = reco.AsVector() + norm_charge * (points[l].AsVector() - reco.AsVector());
                            lines[l]->SetPoint(k,line_point.x,line_point.y,line_point.z);
                        }
                    }
                }
            }
        }

        // Margin.
        vmin -= {0.2,0.2,0.2};
        vmax += {0.2,0.2,0.2};

        // Old code for drawing markers.
        // auto reco_markers = GetDataMarkers(reco_data);
        // for (auto m : reco_markers) m->Draw("same");
        
        if (m_loop->make_track_plots)
        {
            std::string c_reco_name = "c_track_xyz_" + std::to_string(m_loop->curr_track_index);
            c_reco = new TCanvas(c_reco_name.c_str(),"Electron track reconstruction with pads and time bins",700*1.33*1.5,500*1.5);
            ApplyThesisStyle(c_reco);
            c_reco->SetRightMargin(0.15);

            // A histogram for scalling of the axes.
            TH3F* scale = new TH3F("scale",";x [cm];y [cm];z [cm]",1,vmin.x,vmax.x,1,vmin.y,vmax.y,1,vmin.z,vmax.z);
            ApplyThesisStyle(scale);
            scale->SetStats(0);
            scale->SetBinContent(1,0);
            scale->SetMinimum(0);
            scale->SetMaximum(max_charge);
            scale->GetXaxis()->SetTitleOffset(1.1);
            scale->GetYaxis()->SetTitleOffset(1.15);
            scale->GetZaxis()->SetTitleOffset(0.8);
            scale->GetYaxis()->SetNdivisions(505);
            scale->Draw("BOX2Z");
            // TPaletteAxis* palette = (TPaletteAxis*)scale->GetListOfFunctions()->FindObject("palette");
            // palette->SetX1NDC(0.9);
            // palette->SetX2NDC(0.92);
            
            g_xyz_reco->SetTitle("Electron track reconstruction;x [cm];y [cm];z [cm]");
            g_xyz_reco->SetMarkerStyle(2);
            g_xyz_reco->SetMarkerSize(0.4);
            g_xyz_reco->Draw("p same");

            g_xyz->SetMarkerColor(2);
            g_xyz->SetMarkerStyle(7);
            g_xyz->SetMarkerSize(1.2);
            g_xyz->Draw("p same");

            for (auto line : marker_lines) line->Draw("same");

            c_reco->Write();
            
            delete c_reco;
            delete scale;
            
            // DefaultLayout::GetDefaultLayout().DrawPads3D(height);
        }

        for (auto line : marker_lines) delete line;
    }

public:
    RecoPadsTask(double pad_height = -2.5) : height(pad_height) {}
};