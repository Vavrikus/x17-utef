// ROOT dependencies
#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1F.h>
#include <TRandom3.h>

// X17 dependencies
#include "Matrix.h"
#include "Points.h"

int main(int argc, char* argv[])
{
    TApplication app("app", &argc, argv);

    int n_sets = 50000;
    int n_points = 100;
    TH1F* hist = new TH1F("hist","hist",150,-7.5,7.5);
    TH1F* hist2 = new TH1F("hist2","hist2",300,0,30);
    TRandom3* rand = new TRandom3(0);

    for (int i = 0; i < n_sets; i++)
    {
        std::vector<X17::EndPoint> points;
        points.reserve(n_points);

        X17::EndPoint average(0,0,0,0);
        X17::Matrix<4,4> cov(1);
        X17::MapPoint current(100,average,cov);

        for (int j = 0; j < n_points; j++)
        {
            points.push_back(current.GetRandomPoint(rand));

            // X17::EndPoint point(rand->Gaus(0,1),rand->Gaus(0,1),rand->Gaus(0,1),rand->Gaus(0,1));
            // points.push_back(point); 
        }

        double A,B;
        X17::MardiaTest(points,A,B);
        hist->Fill(B);
        hist2->Fill(A);
    }

    TCanvas* c = new TCanvas();
    hist->Scale(1./hist->Integral("width"));    
    hist->Draw("hist same");
    TF1* f_gaus = new TF1("gaus", "1/sqrt(2*TMath::Pi()) * exp(-0.5 * x*x)", -5, 5);
    f_gaus->SetLineColor(kBlue);
    f_gaus->Draw("same");

    std::cout << "Mean: " << hist->GetMean() << "\n";
    std::cout << "Sigma: " << hist->GetRMS() << "\n";

    TCanvas* c2 = new TCanvas();
    hist2->Scale(1./hist2->Integral("width"));
    hist2->Draw("hist same");
    int k = 3*4*5/6; // Degrees of freedom
    TF1 *chi2_pdf = new TF1("chi2_pdf", "[0] * TMath::Power(x, [1]/2 - 1) * TMath::Exp(-x/2) / (TMath::Gamma([1]/2) * TMath::Power(2, [1]/2))", 0, 30);
    chi2_pdf->SetParameters(1.0, k); // [0] = normalization, [1] = k (dof)
    chi2_pdf->SetLineColor(kRed);
    chi2_pdf->Draw("same");

    app.Run();
    return 0;
}