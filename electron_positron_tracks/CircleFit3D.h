#pragma once

#include <algorithm>
#include <vector>

#include "TGraph.h"
#include "TGraph2D.h"
#include "TMarker3DBox.h"
#include "TVirtualFitter.h"

#include "../VectorField.h"
#include "../X17Utilities.h"

using namespace std;

template<typename... T>
double find_min(T... args) {
    double values[] = { args... };
    return *std::min_element(values, values + sizeof...(args));
}

/// @brief Struct for information about reconstructed point for fitting
struct DataPoint
{
    double x,y,z;
    int count; // Current equivalent to charge -- number of electrons in given point

    DataPoint(double x, double y, double z, int c) : x(x),y(y),z(z),count(c) {}
    Vector ToVector() const {return {x,y,z};}
};

vector<TMarker3DBox*> GetDataMarkers(vector<DataPoint> data, double zbin_size = 0.3)
{
    vector<TMarker3DBox*> markers;
    constexpr double max_size = 0.75;

    // find maximal count
    int max_count = 0;
    for (DataPoint p : data) if (p.count > max_count) max_count = p.count;

    // create markers
    for (DataPoint p : data)
    {
        double rel_size = max_size*p.count/max_count;
        double xlen = rel_size*X17::pad_width/2.0;
        double ylen = rel_size*X17::pad_height/2.0;
        double zlen = rel_size*zbin_size/2.0;

        markers.push_back(new TMarker3DBox(p.x,p.y,p.z,xlen,ylen,zlen,0,0));
    }
    
    return markers;
}

/// @brief class for fitting reconstructed track with circular arc with smoothly attached lines
class CircleFit3D
{
    typedef void(*EvalFn)(Int_t&, Double_t*, Double_t&, Double_t*, Int_t);

private:
    vector<DataPoint> fit_data; // Contains recontructed points
    Vector origin,orientation;  // Parameters for the first line (to be fixed in the fit)
    double length;              // Length of the first line (also t_max)
    double alpha;               // Rotation angle of the arc around the first line (if 0 curves towards negative z)
    double radius;              // Radius of the arc
    double phi_max;             // Maximal angle on the arc

    double cos_theta,sin_theta,cos_phi,sin_phi,cos_alpha,sin_alpha; // angle sines and cosines
    Vector origin2,orientation2;  // Parameters calculated for the second line
    Vector originc,center,normal; // Parameters calculated for the circle

    TVirtualFitter* gFitter;
    double l_err,a_err,r_err,phi_err; // Fit error variables

    CircleFit3D() = default;
    CircleFit3D(const CircleFit3D&) = delete;
    CircleFit3D& operator=(const CircleFit3D&) = delete;

    void UpdateCurve()
    {
        cos_theta = orientation.vz;
        sin_theta = sqrt(1-cos_theta*cos_theta);

        cos_phi = orientation.vx/sin_theta;
        sin_phi = orientation.vy/sin_theta;

        cos_alpha = cos(alpha);
        sin_alpha = sin(alpha);

        originc = this->GetLinePoint(length,true);

        normal = Vector {-sin_alpha*cos_theta*cos_phi-cos_alpha*sin_phi,
                         -sin_alpha*cos_theta*sin_phi+cos_alpha*cos_phi,
                         sin_alpha*sin_theta};
        
        center = originc + radius * Vector {cos_alpha*cos_theta*cos_phi-sin_alpha*sin_phi,
                                            cos_alpha*cos_theta*sin_phi+sin_alpha*cos_phi,
                                            -cos_alpha*sin_theta};
        
        origin2 = this->GetCirclePoint(phi_max);

        double cos_varphi = cos(phi_max);
        double sin_varphi = sin(phi_max);
        orientation2 = Vector{sin_varphi*(cos_alpha*cos_theta*cos_phi-sin_alpha*sin_phi)+cos_varphi*sin_theta*cos_phi,
                              sin_varphi*(cos_alpha*cos_theta*sin_phi+sin_alpha*cos_phi)+cos_varphi*sin_theta*sin_phi,
                              -sin_varphi*cos_alpha*sin_theta+cos_varphi*cos_theta};
    }

    double LineSqDist(const DataPoint& point, bool first_line)
    {
        Vector *orig,*orient;
        if (first_line) {orig = &origin;  orient = &orientation; }
        else            {orig = &origin2; orient = &orientation2;}

        double t_close = *orient*(point.ToVector()-*orig)/orient->SqMagnitude();
        if (first_line  && (t_close > length)) t_close = length;
        if (!first_line && (t_close < 0))      t_close = 0;

        Vector line_vector = *orig+t_close*(*orient)-point.ToVector();

        return line_vector.SqMagnitude();
    }

    double CircleSqDist(const DataPoint& point)
    {
        Vector projection = (point.ToVector()-center)-((point.ToVector()-center)*normal)*normal;
        projection.Normalize();

        Vector circle_point   = radius*projection+center;
        Vector circle_vector  = circle_point-point.ToVector();
        Vector arc_beg_vector = originc-point.ToVector();
        Vector arc_end_vector = origin2-point.ToVector();

        return find_min(circle_vector.SqMagnitude(),arc_beg_vector.SqMagnitude(),arc_end_vector.SqMagnitude());
    }

    double SqDistance(const DataPoint& point)
    {
        double sqdist1 = LineSqDist(point,true);
        double sqdist2 = CircleSqDist(point);
        double sqdist3 = LineSqDist(point,false);

        return find_min(sqdist1,sqdist2,sqdist3);
    }

    double SumSq()
    {
        double sum = 0;
        for (DataPoint point : fit_data) sum += point.count*SqDistance(point);
        return sum;
    }

    Vector GetLinePoint(double param, bool first_line)
    {
        if (first_line) return origin+param*orientation;
        else return origin2+param*orientation2;
    }

    Vector GetCirclePoint(double varphi)
    {
        double cos_varphi = cos(varphi);
        double sin_varphi = sin(varphi);

        return originc + radius*Vector{(1-cos_varphi)*(cos_alpha*cos_theta*cos_phi-sin_alpha*sin_phi)+sin_varphi*sin_theta*cos_phi,
                                       (1-cos_varphi)*(cos_alpha*cos_theta*sin_phi+sin_alpha*cos_phi)+sin_varphi*sin_theta*sin_phi,
                                       -(1-cos_varphi)*cos_alpha*sin_theta+sin_varphi*cos_theta};
    }

    void EvalSumSq(int& npar, double* gin, double& sumsq, double* par, int iflag)
    {
        length  = par[0];
        alpha   = par[1];
        radius  = par[2];
        phi_max = par[3];

        UpdateCurve();

        sumsq = SumSq();
    }

    static void Eval(int& npar, double* gin, double& sumsq, double* par, int iflag)
    {
        return GetCircleFit().EvalSumSq(npar,gin,sumsq,par,iflag);
    }

    Vector GetAvgField(VectorField* magfield, double step = 0.1)
    {
        int i = 0;
        double param = 0;
        Vector bfield = {0,0,0};

        while (param/radius < phi_max) 
        {
            Vector cur_point = GetCirclePoint(param/radius);
            if(X17::IsInSector(cur_point))
            {
                bfield += magfield->GetField(cur_point/100.0);
                i++;
            }
            param += step;
        }

        return bfield/i;
    }

    Vector GetMiddleField(VectorField* magfield, double tolerance = 0.0001)
    {
        double xmiddle = (X17::xmax+X17::xmin)/2;
        double low  = 0;
        double high = phi_max;
        double mid  = (low+high)/2;
        Vector vmid = GetCirclePoint(mid);

        while (abs(xmiddle-vmid.vx) > tolerance)
        {
            mid = (low+high)/2;
            vmid  = GetCirclePoint(mid);

            if (vmid.vx < xmiddle) low  = mid;
            else                high = mid;
        }
        
        return magfield->GetField(vmid/100.0);
    }

public:
    static CircleFit3D& GetCircleFit()
    {
        static CircleFit3D instance;
        return instance;
    }

    static CircleFit3D& NewCircleFit(Vector orig, Vector orient)
    {
        CircleFit3D& instance = GetCircleFit();
        instance.fit_data.clear();
        instance.origin      = orig;
        instance.orientation = orient;
        instance.orientation.Normalize();

        return instance;
    }

    vector<DataPoint> GetData() {return fit_data;}

    void AddPoint(double x, double y, double z, int count) {fit_data.emplace_back(x,y,z,count);}

    void AddPoint(DataPoint p) {fit_data.push_back(p);}

    void Prefit()
    {
        // preset free parameters here
        length  = 6.5; 
        alpha   = 0;
        radius  = 20;
        phi_max = 0.3;
        
        
        this->UpdateCurve();
    }
    
    void SetFitter(int parameters = 4, bool print = true)
    {        
        gFitter = TVirtualFitter::Fitter(nullptr,parameters); // the second number is number of parameters
        gFitter->SetFCN(this->Eval);

        if(!print)
        {
            double arg = -1;
            gFitter->ExecuteCommand("SET PRINTOUT",&arg,1);
            gFitter->ExecuteCommand("SET NOW", &arg ,1);
        }
    }

    void FitCircle3D(double max_iter = 500, double toleration = 0.001)
    {
        Prefit();
        gFitter->SetParameter(0,"length",length,0.01,0,20);
        gFitter->SetParameter(1,"alpha",alpha,0.001,0,2*M_PI);
        gFitter->SetParameter(2,"radius",radius,0.01,0,100);
        gFitter->SetParameter(3,"phi_max",phi_max,0.001,0,2*M_PI);

        double arglist[2] = {max_iter,toleration};  // max iterations, step size (toleration)
        gFitter->ExecuteCommand("MIGRAD",arglist,2); // last one num of prints (verbosity)

        length  = gFitter->GetParameter(0);
        l_err   = gFitter->GetParError(0);
        alpha   = gFitter->GetParameter(1);
        a_err   = gFitter->GetParError(1);
        radius  = gFitter->GetParameter(2);
        r_err   = gFitter->GetParError(2);
        phi_max = gFitter->GetParameter(3);
        phi_err = gFitter->GetParError(3);

        UpdateCurve();
    }

    void PrintFitParams()
    {
        cout << "\nCIRCLE FIT PARAMETERS:\n";
        cout << "Length:  " << length  << " +- " << l_err   << "\n";
        cout << "Alpha:   " << alpha   << " +- " << a_err   << "\n";
        cout << "Radius:  " << radius  << " +- " << r_err   << "\n";
        cout << "Phi_max: " << phi_max << " +- " << phi_err << "\n\n";
    }

    TGraph2D* GetGraph(double step = 0.1, double dist = 0)
    {
        TGraph2D* fit_graph = new TGraph2D();
        
        // 1st line
        double param = 0;
        while (param < length)
        {
            Vector cur_pos = GetLinePoint(param,true);
            if(X17::IsInSector(cur_pos,dist)) 
                fit_graph->AddPoint(cur_pos.vx,cur_pos.vy,cur_pos.vz);
            param += step;
        }

        // circle
        param = 0;
        while (param/radius < phi_max)
        {
            Vector cur_pos = GetCirclePoint(param/radius);
            if(X17::IsInSector(cur_pos,dist)) 
                fit_graph->AddPoint(cur_pos.vx,cur_pos.vy,cur_pos.vz);
            param += step;
        }

        // 2nd line
        param = 0;
        while (param < 20)
        {
            Vector cur_pos = GetLinePoint(param,false);
            if(X17::IsInSector(cur_pos,dist)) 
                fit_graph->AddPoint(cur_pos.vx,cur_pos.vy,cur_pos.vz);
            param += step;
        }

        return fit_graph;
    }

    double GetEnergy(VectorField* magfield, bool middle = true)
    {
        Vector bfield;
        if(middle) bfield = GetMiddleField(magfield); //magfield->GetField(GetCirclePoint(phi_max/2.0)/100.0);
        else       bfield = GetAvgField(magfield);

        double b_proj = normal*bfield;
        const double clight = 299792458;
        const double E0 = 510998.95;

        double betasq = 1/(1+pow((E0/(clight*(radius/100.0)*b_proj)),2));
        double Ekin = E0*(1/sqrt(1-betasq)-1);

        return Ekin;
    }

    TGraph* GetEnergyGraph(VectorField* magfield, double step = 0.1)
    {
        TGraph* graph = new TGraph();
        double param = 0;

        while (param/radius < phi_max) 
        {
            Vector cur_point = GetCirclePoint(param/radius);
            if(X17::IsInSector(cur_point))
            {
                Vector bfield = magfield->GetField(cur_point/100.0);
                double b_proj = normal*bfield;
                const double clight = 299792458;
                const double E0 = 510998.95;

                double betasq = 1/(1+pow((E0/(clight*(radius/100.0)*b_proj)),2));
                double Ekin = E0*(1/sqrt(1-betasq)-1);
                graph->AddPoint(param,Ekin);
            }
            param += step;
        }
        
        graph->SetTitle("Energy along the track;Parameter [cm];Energy [eV]");
        return graph;
    }
};
