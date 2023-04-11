#include "TF1.h"
#include "TGraph.h"
#include "TSpline.h"

#include "VectorField.h"

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