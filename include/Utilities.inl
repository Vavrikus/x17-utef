// ROOT dependencies
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMultiGraph.h"

template <typename T>
void ApplyThesisStyle(T* obj)
{
    if constexpr (is_from<T, TCanvas, TPad>)
    {
        double w = obj->GetWw();
        if (w > 800)
        {
            obj->SetLeftMargin(0.1);
            obj->SetRightMargin(0.05);
        }
        else
        {
            obj->SetLeftMargin(0.13);
            obj->SetRightMargin(0.07);
        }
        obj->SetTopMargin(0.07);
        obj->SetBottomMargin(0.13);
    }

    else
    {
        std::vector<TAxis*> axes;
        if constexpr (is_from<T, TH1F, TGraph, TGraphErrors, TMultiGraph>)
            axes = { obj->GetXaxis(), obj->GetYaxis(), nullptr };

        else if constexpr (is_from<T, TH2F, TGraph2D, TH3F>)
            axes = { obj->GetXaxis(), obj->GetYaxis(), obj->GetZaxis() };
        
        else
        {
            std::cout << "Warning: Can't apply thesis style for object of type " << typeid(obj).name() << std::endl;
            return;
        }

        for (TAxis* axis : axes) if (axis)
        {
            axis->SetTitleSize(0.06);
            axis->SetLabelSize(0.055);
        }

        if (axes[2])
        {
            axes[0]->SetTitleOffset(1.2);
            axes[1]->SetTitleOffset(1.2); // previously also 0.9

            double zmax = axes[2]->GetXmax();
            double zmin = axes[2]->GetXmin();
            if (abs(zmin) > zmax) zmax = abs(zmin);
            if (zmax < 9) zmax = 9;
            double z_offset = std::ceil(std::log10(zmax)) * 0.2 + 0.5;
            if (zmin < 0) z_offset += 0.2;
            axes[2]->SetTitleOffset(z_offset); // previously 1.3 for times in 1000s of ns
        }

        else
        {
            axes[0]->SetTitleOffset(0.9);
            axes[1]->SetTitleOffset(0.75);
        }   
    }
}