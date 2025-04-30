#include <iostream>

#include <TApplication.h>
#include <TArrow.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TROOT.h>

#include "PadLayout.h"
#include "X17Utilities.h"

void SetArrowStyle(TArrow* arrow)
{
    arrow->SetLineColor(kRed);
    arrow->SetFillColor(kRed);
    arrow->SetLineWidth(2);
}

void DashedLine(double x1, double y1, double x2, double y2)
{
    TLine *line = new TLine(x1, y1, x2, y2);
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    line->SetLineStyle(2);
    line->Draw();
}

int main(int argc, char *argv[])
{
    using namespace X17::constants;

    TApplication app("app", &argc, argv);

    X17::DefaultLayout& layout = X17::DefaultLayout::GetDefaultLayout();

    const double padding = 1;
    TCanvas* c = new TCanvas("c", "c", 600 * 2 * (yhigh + padding)/(xmax - xmin + 2 * padding), 600);
    c->SetLeftMargin(0.05);
    c->SetRightMargin(0.05);
    
    TH2F* h = new TH2F("h",";y [cm];x [cm]", 1, -yhigh-padding, yhigh+padding, 1, xmin-padding, xmax+padding);
    h->SetStats(0);
    h->Draw("");
    h->GetYaxis()->SetTitleOffset(0.5);
    layout.DrawPads(false,true,c);

    // PAD DIMENSIONS
        double x1_1, y1_1, x1_2, y1_2;
        layout.GetPadCorners(1, x1_1, y1_1, x1_2, y1_2, false);

        DashedLine(y1_1, x1_2, y1_1, x1_2 + 0.5);
        DashedLine(y1_2, x1_2, y1_2, x1_2 + 0.5);

        TArrow *arrow = new TArrow(y1_1, x1_2 + 0.5, y1_2, x1_2 + 0.5, 0.01, "|<|>");
        SetArrowStyle(arrow);
        arrow->Draw();

        // Add text to the right of the arrow indicating height
        TLatex *heightText = new TLatex(y1_2 + 0.2, x1_2 + 0.5, "0.9 cm");
        heightText->SetTextColor(kRed);
        heightText->SetTextSize(0.035);
        heightText->SetTextAlign(12);  // Left alignment
        heightText->Draw();

        double x4_1, y4_1, x4_2, y4_2;
        layout.GetPadCorners(4, x4_1, y4_1, x4_2, y4_2, false);

        DashedLine(y4_1 - 0.5, x4_1, y4_1, x4_1);
        DashedLine(y4_1 - 0.5, x4_2, y4_1, x4_2);

        TArrow *arrow2 = new TArrow(y4_1 - 0.5, x4_1, y4_1 - 0.5, x4_2, 0.01, "|<|>");
        SetArrowStyle(arrow2);
        arrow2->Draw();

        // Add text to the left of the arrow indicating width
        TLatex *widthText = new TLatex(y4_1 - 0.8, (x4_1 + x4_2)/2, "0.6 cm");
        widthText->SetTextColor(kRed);
        widthText->SetTextSize(0.035);
        widthText->SetTextAlign(32);  // Right alignment
        widthText->Draw();

    // SPECIAL PADS
        for (int pad : {102, 124, 127}) {
            double x1, y1, x2, y2;
            layout.GetPadCorners(pad, x1, y1, x2, y2, false);

            DashedLine(y1, x1, y1, x1 - 1.1);
            DashedLine(y2, x1, y2, x1 - 1.1);
            
            TArrow *arrow3 = new TArrow(y1, x1 - 1.1, y2, x1 - 1.1, 0.01, "|<|>");
            SetArrowStyle(arrow3);
            arrow3->Draw();

            TString s_height = pad == 127 ? "0.509 cm" : "0.6 cm";

            // Add text to the right of the arrow indicating height
            TLatex *heightText = new TLatex(y2 + 0.2, x1 - 1.1, s_height);
            heightText->SetTextColor(kRed);
            heightText->SetTextSize(0.035);
            heightText->SetTextAlign(12);  // Left alignment
            heightText->Draw();
        }

    // STAGGERING OFFSET
        // Get pad coordinates for pads 6 and 7
        double x6_1, y6_1, x6_2, y6_2;
        double x7_1, y7_1, x7_2, y7_2;
        layout.GetPadCorners(6, x6_1, y6_1, x6_2, y6_2, false);
        layout.GetPadCorners(7, x7_1, y7_1, x7_2, y7_2, false);

        // Calculate centers
        double cx6 = (x6_1 + x6_2)/2;
        double cy6 = (y6_1 + y6_2)/2;
        double cx7 = (x7_1 + x7_2)/2;
        double cy7 = (y7_1 + y7_2)/2;

        // Create measurement line with two-sided arrow
        DashedLine(y6_1, x6_1, y6_1, x7_1);

        TArrow *stagger_arrow = new TArrow(y6_1, x7_1, y7_1, x7_1, 0.01, "|<|>");
        SetArrowStyle(stagger_arrow);

        // Annotation text (bottom-left of the line)
        TLatex *text = new TLatex();
        text->SetTextColor(kRed);
        text->SetTextSize(0.035);
        text->SetTextAlign(33);  // Bottom-right alignment

        // Calculate text position (midpoint shifted down-left)
        const double text_x = (cy6 + cy7)/2;
        const double text_y = (cx6 + cx7)/2 - 0.9;  // Shift down

        // Draw all elements
        stagger_arrow->Draw();
        text->DrawLatex(text_x, text_y, "Offset: 0.3946 cm");

    // GAP
        double x10_1, y10_1, x10_2, y10_2;
        double x11_1, y11_1, x11_2, y11_2;
        layout.GetPadCorners(10, x10_1, y10_1, x10_2, y10_2, false);
        layout.GetPadCorners(11, x11_1, y11_1, x11_2, y11_2, false);

        DashedLine(y10_1, x10_1, y7_1, x10_1);
        DashedLine(y11_1, x11_2, y7_1, x11_2);

        TArrow *gap_arrow1 = new TArrow(y7_1, x10_1 + 0.5, y7_1, x10_1, 0.01, "|>");
        SetArrowStyle(gap_arrow1);
        gap_arrow1->Draw();

        TArrow *gap_arrow2 = new TArrow(y7_1, x11_2 - 0.5, y7_1, x11_2, 0.01, "|>");
        SetArrowStyle(gap_arrow2);
        gap_arrow2->Draw();
        
        // Text for gap
        TLatex *gap_text = new TLatex();
        gap_text->SetTextSize(0.035);
        gap_text->SetTextColor(kRed);
        gap_text->SetTextAlign(32);
        gap_text->DrawLatex(y7_1 - 0.25, (x10_1 + x11_2)/2, "Gap: 0.08 cm");

        double x49_1, y49_1, x49_2, y49_2;
        double x61_1, y61_1, x61_2, y61_2;
        layout.GetPadCorners(49, x49_1, y49_1, x49_2, y49_2, false);
        layout.GetPadCorners(61, x61_1, y61_1, x61_2, y61_2, false);

        DashedLine(y49_2, x49_2, y49_2, x49_2 + 0.5);
        DashedLine(y61_1, x61_2, y61_1, x61_2 + 0.5);

        TArrow *gap_arrow3 = new TArrow(y49_2 - 0.5, x49_2 + 0.5, y49_2, x49_2 + 0.5, 0.01, "|>");
        SetArrowStyle(gap_arrow3);
        gap_arrow3->Draw();

        TArrow *gap_arrow4 = new TArrow(y61_1 + 0.5, x61_2 + 0.5, y61_1, x61_2 + 0.5, 0.01, "|>");
        SetArrowStyle(gap_arrow4);
        gap_arrow4->Draw();

        TLatex *gap_text2 = new TLatex(y61_1 + 0.7, x61_2 + 0.5, "Gap: 0.08 cm");
        gap_text2->SetTextColor(kRed);
        gap_text2->SetTextSize(0.035);
        gap_text2->SetTextAlign(12);  // Left alignment
        gap_text2->Draw();

    c->Modified();
    c->Update();

    app.Run(true);
}

