// C++ dependencies
#include <cmath>

// ROOT dependencies
#include "TF1.h"
#include "TGraph.h"
#include "TSpline.h"

// X17 dependencies
#include "CircleFit2D.h"
#include "Field.h"
#include "Vector.h"
#include "X17Utilities.h"

namespace X17
{
    TF1* FitCircle(TGraph* graph, double min, double max)
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
        double& xx      = x[0];
        double& radius  = par[0];
        double& node1_x = par[1];
        double& node2_x = par[2];
        double& node1_y = par[3];
        double& node2_y = par[4];

        double x_mid = (node2_x - node1_x) / 2;
        double y_mid = (node2_y - node1_y) / 2;
        double r_mid = std::sqrt(pow(x_mid,2) + pow(y_mid,2));

        double r_c = std::sqrt(pow(radius,2) - pow(r_mid,2));
        double x0 = node1_x + x_mid + y_mid * r_c / r_mid;
        double y0 = node1_y + y_mid - x_mid * r_c / r_mid;

        double a1 = (x0 - node1_x) / std::sqrt(pow(radius,2) - pow(node1_x - x0,2));
        double a2 = (x0 - node2_x) / std::sqrt(pow(radius,2) - pow(node2_x - x0,2));
        double b1 = node1_y - a1 * node1_x;
        double b2 = node2_y - a2 * node2_x;

        if (xx < node1_x) return a1 * xx + b1;
        if (xx > node2_x) return a2 * xx + b2;    
        return y0 + std::sqrt(pow(radius,2) - pow(xx - x0,2));
    }

    TF1* FitCircle2(TGraph* graph, double min, double max)
    {
        TF1* circle = new TF1("circle",circle_func,min,max,5);
        circle->SetParameter(0,20);
        circle->SetParLimits(0,5,500);
        circle->SetParameter(1,4);
        circle->SetParLimits(1,0,7);
        circle->SetParameter(2,12);
        circle->SetParLimits(2,7,20);
        circle->SetParameter(3,0);
        circle->SetParameter(4,-2);
        
        graph->Fit(circle,"M","",min,max);
        return circle;
    }

    void RecoEnergy(TF1* fit, const Field<Vector>& magfield, TGraph* magnetic, double min, double max, double step)
    {
        using namespace constants;

        double r = cm2m * fit->GetParameter(0);

        // Mean magnetic field.
        Vector B = magfield.GetField((xmax + xmin) / 2, 0, fit->Eval((xmax + xmin) / 2));

        for (double x = min; x <= max; x += step)
        {
            Vector B2 = magfield.GetField(x, 0, fit->Eval(x));
            magnetic->AddPoint(x,B2.y);
        }
        
        double betasq = 1 / (1 + pow((E0 / (c * r * B.y)), 2));
        double Ekin = E0 * (1 / std::sqrt(1 - betasq) - 1);

        std::cout << "Kinetic energy: " << Ekin << " eV, By: " << B.y << " T, beta: ";
        std::cout << std::sqrt(betasq) << "\n";
    }
}