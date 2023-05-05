#pragma once

#include "TGraph.h"
#include "TGraph2D.h"
#include "TVirtualFitter.h"

#include "Field.h"
#include "Points.h"
#include "Utilities.h"
#include "Vector.h"
#include "X17Utilities.h"

namespace X17
{
    /// @brief class for fitting reconstructed track with circular arc with smoothly attached lines
    class CircleFit3D
    {
        typedef void(*EvalFn)(Int_t&, Double_t*, Double_t&, Double_t*, Int_t);

    private:
        static CircleFit3D* lastfit; // Pointer to the last instance.

        std::vector<RecoPoint> fit_data; // A std::vector of reconstructed points to be fitted.
        Vector origin;                   // The origin of the first line [cm].
        Vector orientation;              // The orientation vector of the first line. Used for theta and phi calculation.
        double length;                   // The length of the first line (also max parameter value) [cm].
        double alpha;                    // The rotation angle of the arc around the first line (if 0 curves towards negative z) [rad].
        double radius;                   // The radius of the arc [rad].
        double phi_max;                  // The maximal angle on the arc [rad].

        double cos_theta; // The cosine of theta.
        double sin_theta; // The sine of theta.
        double cos_phi;   // The cosine of phi.
        double sin_phi;   // The sine of phi.
        double cos_alpha; // The cosine of alpha.
        double sin_alpha; // The sine of alpha.

        Vector origin2;      // The origin of the second line (to be calculated).
        Vector orientation2; // The orientation vector of the second line (to be calculated).

        Vector originc; // The origin of the circular arc (to be calculated).
        Vector center;  // The center of the circular arc (to be calculated).
        Vector normal;  // The normal to the plane of the circular arc (to be calculated).

        TVirtualFitter* gFitter; // Contains TVirtualFitter instance used for the fit.
        double l_err;            // The length fit error.
        double a_err;            // The alpha fit error.
        double r_err;            // The radius fit error.
        double phi_err;          // The phi fit error.

        /// @brief Returns a point on the line.
        /// @param param The parameter that defines the position on the line.
        /// @param first_line Whether to return a point on the first line or the second line.
        /// @return A point on the specified line.
        Vector GetLinePoint(double param, bool first_line)
        {
            if (first_line) return origin  + param * orientation;
            else            return origin2 + param * orientation2;
        }

        /// @brief Calculates a point on the circle described by the current instance.
        /// @param varphi The angle (in radians) at which to calculate the point. Zero coresponds to the beginning of the arc.
        /// @return A vector representing the point on the circle.
        Vector GetCirclePoint(double varphi)
        {
            double cos_varphi = cos(varphi);
            double sin_varphi = sin(varphi);

            return originc + radius*Vector{ (1-cos_varphi)*(cos_alpha*cos_theta*cos_phi - sin_alpha*sin_phi) + sin_varphi*sin_theta*cos_phi,
                                            (1-cos_varphi)*(cos_alpha*cos_theta*sin_phi + sin_alpha*cos_phi) + sin_varphi*sin_theta*sin_phi,
                                           -(1-cos_varphi)*cos_alpha*sin_theta + sin_varphi*cos_theta};
        }

        /// @brief Updates the curve parameters based on the current orientation, length, radius, alpha, and phi_max values.
        /// This function calculates the new values of cos_theta, sin_theta, cos_phi, sin_phi, cos_alpha, sin_alpha, originc,
        /// normal, center, origin2, and orientation2. These values are used to update the properties of the curve.
        /// @note This function assumes that the orientation, length, radius, alpha, and phi_max values have already been set.
        void UpdateCurve()
        {
            cos_theta = orientation.z;
            sin_theta = sqrt(1-cos_theta*cos_theta);

            cos_phi = orientation.x/sin_theta;
            sin_phi = orientation.y/sin_theta;

            cos_alpha = cos(alpha);
            sin_alpha = sin(alpha);

            originc = this->GetLinePoint(length,true);

            normal = Vector {-sin_alpha*cos_theta*cos_phi - cos_alpha*sin_phi,
                             -sin_alpha*cos_theta*sin_phi + cos_alpha*cos_phi,
                              sin_alpha*sin_theta };
            
            center = originc + radius * Vector { cos_alpha*cos_theta*cos_phi - sin_alpha*sin_phi,
                                                 cos_alpha*cos_theta*sin_phi + sin_alpha*cos_phi,
                                                -cos_alpha*sin_theta };
            
            origin2 = this->GetCirclePoint(phi_max);

            double cos_varphi = cos(phi_max);
            double sin_varphi = sin(phi_max);

            orientation2 = Vector{ sin_varphi*(cos_alpha*cos_theta*cos_phi - sin_alpha*sin_phi) + cos_varphi*sin_theta*cos_phi,
                                   sin_varphi*(cos_alpha*cos_theta*sin_phi + sin_alpha*cos_phi) + cos_varphi*sin_theta*sin_phi,
                                  -sin_varphi*cos_alpha*sin_theta + cos_varphi*cos_theta };
        }

        /// @brief Calculates the squared distance between a given point and a line.
        /// @param point The point to calculate the distance from.
        /// @param first_line Flag indicating whether to use the first line (true) or the second line (false).
        /// @return The squared distance between the point and the line.
        double LineSqDist(const RecoPoint& point, bool first_line)
        {
            Vector *orig,*orient;
            if (first_line) { orig = &origin;  orient = &orientation;  }
            else            { orig = &origin2; orient = &orientation2; }

            double t_close = *orient*(point.AsVector()-*orig)/orient->SqMagnitude();
            if (first_line  && (t_close > length)) t_close = length;
            if (!first_line && (t_close < 0))      t_close = 0;

            Vector line_vector = *orig+t_close*(*orient)-point.AsVector();

            return line_vector.SqMagnitude();
        }

        /// @brief Calculate the squared distance between a given point and the circular trajectory.
        /// @param point The point to calculate the distance from.
        /// @return The squared distance between the point and the circular trajectory.
        double CircleSqDist(const RecoPoint& point)
        {
            Vector projection = (point.AsVector()-center)-((point.AsVector()-center)*normal)*normal;
            projection.Normalize();

            Vector circle_point   = radius * projection + center;
            Vector circle_vector  = circle_point - point.AsVector();
            Vector arc_beg_vector = originc - point.AsVector();
            Vector arc_end_vector = origin2 - point.AsVector();

            return find_min(circle_vector.SqMagnitude(),arc_beg_vector.SqMagnitude(),arc_end_vector.SqMagnitude());
        }

        /// @brief Computes the squared distance between the given point and the curve.
        /// @param point The point to compute the squared distance from.
        /// @return The squared distance between the point and the curve.
        double SqDistance(const RecoPoint& point)
        {
            double sqdist1 = LineSqDist(point,true);
            double sqdist2 = CircleSqDist(point);
            double sqdist3 = LineSqDist(point,false);

            return find_min(sqdist1,sqdist2,sqdist3);
        }

        /// @brief Calculates the sum of squared distances from each point in the fit_data vector to the curve. Each distance is scaled by the count/charge.
        /// @return The sum of squared distances from each point in fit_data to the curve.
        double SumSq()
        {
            double sum = 0;
            for (RecoPoint point : fit_data) sum += point.count*SqDistance(point);
            return sum;
        }

        /// @brief Calculate the sum of squared distances for a given set of parameters.
        /// @param npar The number of parameters.
        /// @param gin An array of size npar used to store the first derivatives of the function with respect to the parameters.
        /// @param sumsq The sum of squared distances, which is calculated by this function.
        /// @param par An array of size npar that contains the values of the parameters.
        /// @param iflag An integer that can be used to control the behavior of the function.
        void EvalSumSq(int& npar, double* gin, double& sumsq, double* par, int iflag)
        {
            length  = par[0];
            alpha   = par[1];
            radius  = par[2];
            phi_max = par[3];

            UpdateCurve();

            sumsq = SumSq();
        }

        /// @brief A static wrapper function used to get the correct function pointer type needed for ROOT. Simply calls lastfit->EvalSumSq().
        /// @param npar The number of parameters.
        /// @param gin An array of size npar used to store the first derivatives of the function with respect to the parameters.
        /// @param sumsq The sum of squared distances, which is calculated by this function.
        /// @param par An array of size npar that contains the values of the parameters.
        /// @param iflag An integer that can be used to control the behavior of the function.
        static void Eval(int& npar, double* gin, double& sumsq, double* par, int iflag)
        {
            return lastfit->EvalSumSq(npar,gin,sumsq,par,iflag);
        }

        /// @brief Computes the average magnetic field vector over the circular path defined by the
        ///        current object's parameters, using the provided magnetic field vector field.
        /// @param magfield Pointer to an object containing the magnetic field information.
        /// @param step Step size to use when sampling the circular path. Default value is 0.1.
        /// @return The average magnetic field vector over the circular path.
        Vector GetAvgField(const Field<Vector>& magfield, double step = 0.1)
        {
            int i = 0;
            double param = 0;
            Vector bfield = {0,0,0};

            while (param/radius < phi_max) 
            {
                Vector cur_point = GetCirclePoint(param/radius);
                if(IsInSector(cur_point))
                {
                    bfield += magfield.GetField(cur_point/100.0);
                    i++;
                }
                param += step;
            }

            return bfield/i;
        }

        /// @brief Get the magnetic field at the point on the circle that passes through the middle of the X17 detector in the x direction.
        /// @param magfield The magnetic field vector field to sample from.
        /// @param tolerance The maximum deviation from the middle of the detector in the x direction for which to stop the binary search. Default is 0.0001.
        /// @return The magnetic field vector at the point on the circle that passes through the middle of the X17 detector in the x direction.
        Vector GetMiddleField(const Field<Vector>& magfield, double tolerance = 0.0001)
        {
            using namespace constants;

            double xmiddle = (xmax + xmin) / 2;
            double low  = 0;
            double high = phi_max;
            double mid  = (low + high) / 2;
            Vector vmid = GetCirclePoint(mid);

            while ((abs(xmiddle-vmid.x) > tolerance) && (low*(1+1e-15) < high))
            {
                mid = (low + high) / 2;
                vmid  = GetCirclePoint(mid);

                if (vmid.x < xmiddle) low  = mid;
                else                   high = mid;
            }
            
            return magfield.GetField(vmid/100.0);
        }

    public:
        /// @brief The default constructor for the 3D circle fitting algorithm that smoothly attaches lines.
        CircleFit3D() { lastfit = this; }

        /// @brief Constructor for the 3D circle fitting algorithm that smoothly attaches lines.
        /// @param orig The origin of the first line [cm].
        /// @param orient The orientation vector of the first line.
        CircleFit3D(const Vector& orig, const Vector& orient)
        {
            this->origin      = orig;
            this->orientation = orient;
            this->orientation.Normalize();

            this->gFitter = lastfit->gFitter; // Setting the fitter again every time would take much more time.
            lastfit = this;
        }

        /// @brief Deleted copy constructor. 
        CircleFit3D(const CircleFit3D&) = delete;

        /// @brief Deleted assignment operator. 
        CircleFit3D& operator=(const CircleFit3D&) = delete;

        /// @brief Get the current fit data points.
        /// @return The vector of RecoPoints used for fitting.
        std::vector<RecoPoint> GetData() { return fit_data; }

        /// @brief Add a new RecoPoint to the fit data.
        /// @param x The x-coordinate of the RecoPoint.
        /// @param y The y-coordinate of the RecoPoint.
        /// @param z The z-coordinate of the RecoPoint.
        /// @param count The count associated with the RecoPoint.
        void AddPoint(double x, double y, double z, int count) { fit_data.emplace_back(x,y,z,count); }

        /// @brief Add a new RecoPoint to the fit data.
        /// @param p The RecoPoint to add.
        void AddPoint(RecoPoint p) { fit_data.push_back(p); }

        /// @brief Add a new RecoPoint to the fit data.
        /// @param p The RKPoint with necessary coordinate information.
        void AddPoint(RKPoint p) { fit_data.emplace_back(p.x,p.y,p.z,1); }

        /// @brief Set the alpha angle for the circle fit.
        /// @param electron If true, set alpha to 0; if false, set alpha to pi.
        void SetAlpha(bool electron)
        {
            if (electron) alpha = 0;
            else          alpha = M_PI;
        }

        /// @brief Sets up the TVirtualFitter for the circle fitting with the number of parameters and printout options.
        /// @param parameters Number of fitting parameters.
        /// @param print If true, the fit will printout the fitting status. If false, the fit will be silent.
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

        /// @brief Fits a 3D circle using a maximum of max_iter iterations and the given toleration.
        /// @param max_iter The maximum number of iterations for the fit. Defaults to 500.
        /// @param toleration The toleration for the fit. Defaults to 0.001.
        void FitCircle3D(double max_iter = 500, double toleration = 0.001)
        {
            gFitter->SetParameter(0,"length",length,0.01,-10,5);
            gFitter->SetParameter(1,"alpha",alpha,0.001,-M_PI/2,(3/2)*M_PI);
            gFitter->SetParameter(2,"radius",radius,0.01,10,50);
            gFitter->SetParameter(3,"phi_max",phi_max,0.001,0.15,M_PI/1.5);

            double arglist[2] = {max_iter,toleration};   // max iterations, step size (toleration)
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

        /// @brief Prints the fitting parameters of the circle with their errors.
        void PrintFitParams()
        {
            std::cout << "\nCIRCLE FIT PARAMETERS:\n";
            std::cout << "Length:  " << length  << " +- " << l_err   << "\n";
            std::cout << "Alpha:   " << alpha   << " +- " << a_err   << "\n";
            std::cout << "Radius:  " << radius  << " +- " << r_err   << "\n";
            std::cout << "Phi_max: " << phi_max << " +- " << phi_err << "\n\n";
        }

        /// @brief Generates a TGraph2D representing the fitted circle with attached lines.
        /// @param step The step size between points in the TGraph2D.
        /// @param dist The distance from the TPC walls within which the points are included.
        /// @return A TGraph2D representing the fitted circle with attached lines.
        TGraph2D* GetGraph(double step = 0.1, double dist = 0)
        {
            TGraph2D* fit_graph = new TGraph2D();
            
            // 1st line
            double param = 0;
            while (param < length)
            {
                Vector cur_pos = GetLinePoint(param,true);
                if(IsInSector(cur_pos,dist)) 
                    fit_graph->AddPoint(cur_pos.x,cur_pos.y,cur_pos.z);
                param += step;
            }

            // circle
            param = 0;
            while (param/radius < phi_max)
            {
                Vector cur_pos = GetCirclePoint(param/radius);
                if(IsInSector(cur_pos,dist)) 
                    fit_graph->AddPoint(cur_pos.x,cur_pos.y,cur_pos.z);
                param += step;
            }

            // 2nd line
            param = 0;
            while (param < 20)
            {
                Vector cur_pos = GetLinePoint(param,false);
                if(IsInSector(cur_pos,dist)) 
                    fit_graph->AddPoint(cur_pos.x,cur_pos.y,cur_pos.z);
                param += step;
            }

            return fit_graph;
        }

        /// @brief Calculates the energy of a charged particle moving along the 3D circle fit.
        /// @param magfield The magnetic field data used for the calculation.
        /// @param middle Boolean flag indicating whether to use the middle field of the circle (true) or the average field (false).
        /// @return The kinetic energy of the particle.
        double GetEnergy(const Field<Vector>& magfield, bool middle = true)
        {
            using namespace constants;

            Vector bfield;
            if(middle) bfield = GetMiddleField(magfield);
            else       bfield = GetAvgField(magfield);

            double b_proj = normal*bfield;

            double betasq = 1/(1+pow((E0/(c*(radius/100.0)*b_proj)),2));
            double Ekin = E0*(1/sqrt(1-betasq)-1);

            return Ekin;
        }

        /// @brief Returns a graph representing the energy along the track in the given magnetic field.
        /// @param magfield The magnetic field data used for the calculation.
        /// @param step The step size for the parameter along the track.
        /// @return A TGraph object representing the energy as a function of the parameter.
        TGraph* GetEnergyGraph(const Field<Vector>& magfield, double step = 0.1)
        {
            TGraph* graph = new TGraph();
            double param = 0;

            while (param/radius < phi_max) 
            {
                Vector cur_point = GetCirclePoint(param/radius);
                if(IsInSector(cur_point))
                {
                    Vector bfield = magfield.GetField(cur_point/100.0);
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

    CircleFit3D* CircleFit3D::lastfit = nullptr;
} // namespace X17


////////////////////////////////////////////
// FITTING CIRCLE ONLY, 6 FREE PARAMETERS //
////////////////////////////////////////////
//
// class CircleFit3D
// {
//     typedef void(*EvalFn)(Int_t&, Double_t*, Double_t&, Double_t*, Int_t);
//
// private:
//     vector<DataPoint> fit_data; // Contains recontructed points
//     Vector origin,orientation;  // Parameters for the first line (to be fixed in the fit)
//     double alpha;               // Rotation angle of the arc around the first line (if 0 curves towards negative z)
//     double radius;              // Radius of the arc
//     double y0;                  // Coordinate y on enterance to TPC
//     double z0;                  // Coordinate z on enterance to TPC
//     double theta;               // Angle from z-axis on enterance to TPC
//     double phi;                 // XY plane angle from x-axis on enterance to TPC
//
//     double cos_theta,sin_theta,cos_phi,sin_phi,cos_alpha,sin_alpha; // angle sines and cosines
//     Vector originc,center,normal; // Parameters calculated for the circle
//
//     TVirtualFitter* gFitter;
//     double a_err,r_err,y0_err,z0_err,theta_err,phi_err; // Fit error variables
//
//     CircleFit3D() = default;
//     CircleFit3D(const CircleFit3D&) = delete;
//     CircleFit3D& operator=(const CircleFit3D&) = delete;
//
//     void UpdateCurve()
//     {
//         cos_theta = cos(theta);
//         sin_theta = sin(theta);
//
//         cos_phi = cos(phi);
//         sin_phi = sin(phi);
//
//         cos_alpha = cos(alpha);
//         sin_alpha = sin(alpha);
//
//         originc = {X17::xmin,y0,z0};
//
//         normal = Vector {-sin_alpha*cos_theta*cos_phi-cos_alpha*sin_phi,
//                          -sin_alpha*cos_theta*sin_phi+cos_alpha*cos_phi,
//                          sin_alpha*sin_theta};
//      
//         center = originc + radius * Vector {cos_alpha*cos_theta*cos_phi-sin_alpha*sin_phi,
//                                             cos_alpha*cos_theta*sin_phi+sin_alpha*cos_phi,
//                                             -cos_alpha*sin_theta};
//     }
//
//     double SqDistance(const DataPoint& point)
//     {
//         Vector projection = (point.ToVector()-center)-((point.ToVector()-center)*normal)*normal;
//         projection.Normalize();
//
//         Vector circle_point   = radius*projection+center;
//         Vector circle_vector  = circle_point-point.ToVector();
//
//         return circle_vector.SqMagnitude();
//     }
//
//     double SumSq()
//     {
//         double sum = 0;
//         for (DataPoint point : fit_data) sum += point.count*SqDistance(point);
//         return sum;
//     }
//
//     Vector GetCirclePoint(double varphi)
//     {
//         double cos_varphi = cos(varphi);
//         double sin_varphi = sin(varphi);
//
//         return originc + radius*Vector{(1-cos_varphi)*(cos_alpha*cos_theta*cos_phi-sin_alpha*sin_phi)+sin_varphi*sin_theta*cos_phi,
//                                        (1-cos_varphi)*(cos_alpha*cos_theta*sin_phi+sin_alpha*cos_phi)+sin_varphi*sin_theta*sin_phi,
//                                        -(1-cos_varphi)*cos_alpha*sin_theta+sin_varphi*cos_theta};
//     }
//
//     double GetMaxPhi()
//     {
//         // A*cos + B = C*sin
//         double A = cos_alpha*cos_theta*cos_phi-sin_alpha*sin_phi;
//         double B = (X17::xmax - originc.x)/radius - A;
//         double C = sin_theta*cos_phi;
//
//         // (A^2+C^2) cos^2 + AB cos + (B^2-C^2) = 0
//         double a = A*A+C*C;
//         double b = A*B;
//         double c = B*B-C*C;
//
//         return acos((-b+sqrt(b*b-4*a*c))/(2*a));
//     }
//
//     void EvalSumSq(int& npar, double* gin, double& sumsq, double* par, int iflag)
//     {
//         alpha   = par[0];
//         radius  = par[1];
//         y0      = par[2];
//         z0      = par[3];
//
//         UpdateCurve();
//
//         sumsq = SumSq();
//     }
//
//     static void Eval(int& npar, double* gin, double& sumsq, double* par, int iflag)
//     {
//         return GetCircleFit().EvalSumSq(npar,gin,sumsq,par,iflag);
//     }
//
//     Vector GetAvgField(VectorField* magfield, double step = 0.1)
//     {
//         int i = 0;
//         double param = 0;
//         Vector bfield = {0,0,0};
//         Vector cur_point = GetCirclePoint(param/radius);
//
//         while (X17::IsInSector(cur_point,-0.01)) 
//         {
//             bfield += magfield->GetField(cur_point/100.0);
//             i++;
//             param += step;
//         }
//
//         return bfield/i;
//     }
//
//     Vector GetMiddleField(VectorField* magfield, double tolerance = 0.0001)
//     {
//         double xmiddle = (X17::xmax+X17::xmin)/2;
//         double low  = 0;
//         double high = GetMaxPhi();
//         double mid  = (low+high)/2;
//         Vector vmid = GetCirclePoint(mid);
//
//         while (abs(xmiddle-vmid.x) > tolerance)
//         {
//             mid = (low+high)/2;
//             vmid  = GetCirclePoint(mid);
//
//             if (vmid.x < xmiddle) low  = mid;
//             else                high = mid;
//         }
//      
//         return magfield->GetField(vmid/100.0);
//     }
//
// public:
//     static CircleFit3D& GetCircleFit()
//     {
//         static CircleFit3D instance;
//         return instance;
//     }
//
//     static CircleFit3D& NewCircleFit(Vector orig, Vector orient)
//     {
//         CircleFit3D& instance = GetCircleFit();
//         instance.fit_data.clear();
//         instance.origin      = orig;
//         instance.orientation = orient;
//         instance.orientation.Normalize();
//
//         return instance;
//     }
//
//     vector<DataPoint> GetData() {return fit_data;}
//
//     void AddPoint(double x, double y, double z, int count) {fit_data.emplace_back(x,y,z,count);}
//
//     void AddPoint(DataPoint p) {fit_data.push_back(p);}
//
//     void Prefit()
//     {
//         double param = (X17::xmin - origin.x)/orientation.x;
//         Vector v0 = origin+param*orientation;
//
//         // preset free parameters here
//         alpha  = 0;
//         radius = 20;
//         y0     = v0.y;
//         z0     = v0.z;
//         theta  = acos(orientation.z);
//         phi    = asin(orientation.y/sin(theta));
//
//         this->UpdateCurve();
//     }
//  
//     void SetFitter(int parameters = 6, bool print = true)
//     {        
//         gFitter = TVirtualFitter::Fitter(nullptr,parameters); // the second number is number of parameters
//         gFitter->SetFCN(this->Eval);
//
//         if(!print)
//         {
//             double arg = -1;
//             gFitter->ExecuteCommand("SET PRINTOUT",&arg,1);
//             gFitter->ExecuteCommand("SET NOW", &arg ,1);
//         }
//     }
//
//     void FitCircle3D(double max_iter = 500, double toleration = 0.001)
//     {
//         Prefit();
//         gFitter->SetParameter(0,"alpha",  alpha,  0.001,  0,2*M_PI);
//         gFitter->SetParameter(1,"radius", radius, 0.001,  0,100);
//         gFitter->SetParameter(2,"y0",     y0,     0.01,  -X17::win_width/2,X17::win_width/2);
//         gFitter->SetParameter(3,"z0",     z0,     0.01,  -X17::win_height/2,X17::win_height/2);
//         gFitter->SetParameter(4,"theta",  theta,  0.001,  0,M_PI);
//         gFitter->SetParameter(5,"phi",    phi,    0.001,  0,2*M_PI);
//
//         double arglist[2] = {max_iter,toleration};  // max iterations, step size (toleration)
//         gFitter->ExecuteCommand("MIGRAD",arglist,2); // last one num of prints (verbosity)
//
//         alpha  = gFitter->GetParameter(0);
//         radius = gFitter->GetParameter(1);
//         y0     = gFitter->GetParameter(2);
//         z0     = gFitter->GetParameter(3);
//         theta  = gFitter->GetParameter(4);
//         phi    = gFitter->GetParameter(5);
//
//         a_err     = gFitter->GetParError(0);
//         r_err     = gFitter->GetParError(1);
//         y0_err    = gFitter->GetParError(2);
//         z0_err    = gFitter->GetParError(3);
//         theta_err = gFitter->GetParError(4);
//         phi_err   = gFitter->GetParError(5);
//
//         UpdateCurve();
//     }
//
//     void PrintFitParams()
//     {
//         cout << "\nCIRCLE FIT PARAMETERS:\n";
//         cout << "Alpha:   " << alpha   << " +- " << a_err     << "\n";
//         cout << "Radius:  " << radius  << " +- " << r_err     << "\n";
//         cout << "Y0:      " << y0      << " +- " << y0_err    << "\n";
//         cout << "Z0:      " << z0      << " +- " << z0_err    << "\n";
//         cout << "Theta:   " << theta   << " +- " << theta_err << "\n";
//         cout << "Phi:     " << phi     << " +- " << phi_err   << "\n\n";
//     }
//
//     TGraph2D* GetGraph(double step = 0.1, double dist = 0)
//     {
//         TGraph2D* fit_graph = new TGraph2D();
//      
//         double param = 0;
//         while (param/radius < GetMaxPhi())
//         {
//             Vector cur_pos = GetCirclePoint(param/radius);
//             if(X17::IsInSector(cur_pos,dist)) 
//                 fit_graph->AddPoint(cur_pos.x,cur_pos.y,cur_pos.z);
//             param += step;
//         }
//
//         return fit_graph;
//     }
//
//     double GetEnergy(VectorField* magfield, bool middle = true)
//     {
//         Vector bfield;
//         if(middle) bfield = GetMiddleField(magfield); //magfield->GetField(GetCirclePoint(phi_max/2.0)/100.0);
//         else       bfield = GetAvgField(magfield);
//
//         double b_proj = normal*bfield;
//         const double clight = 299792458;
//         const double E0 = 510998.95;
//
//         double betasq = 1/(1+pow((E0/(clight*(radius/100.0)*b_proj)),2));
//         double Ekin = E0*(1/sqrt(1-betasq)-1);
//
//         return Ekin;
//     }
//
//     TGraph* GetEnergyGraph(VectorField* magfield, double step = 0.1)
//     {
//         TGraph* graph = new TGraph();
//         double param = 0;
//
//         while (param/radius < GetMaxPhi()) 
//         {
//             Vector cur_point = GetCirclePoint(param/radius);
//             if(X17::IsInSector(cur_point))
//             {
//                 Vector bfield = magfield->GetField(cur_point/100.0);
//                 double b_proj = normal*bfield;
//                 const double clight = 299792458;
//                 const double E0 = 510998.95;
//
//                 double betasq = 1/(1+pow((E0/(clight*(radius/100.0)*b_proj)),2));
//                 double Ekin = E0*(1/sqrt(1-betasq)-1);
//                 graph->AddPoint(param,Ekin);
//             }
//             param += step;
//         }
//      
//         graph->SetTitle("Energy along the track;Parameter [cm];Energy [eV]");
//         return graph;
//     }
// };