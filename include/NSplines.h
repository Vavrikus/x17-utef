#pragma once

#include <TSpline.h>

namespace X17
{
	/// @brief Represents a cubic spline that interpolates a set of N nodes. Defined by the positions of the nodes, derivatives at the start and end of the spline.
	/// @tparam N The number of nodes.
	template<int N>
	class NSpline 
	{
		typedef std::function<double(double*,double*)> EvalFn;

	private:
		double nodes_x[N];        // The x-coordinates of spline nodes.
		bool FixDerStart = false; // Is the derivative at the start fixed in fit?
		bool FixDerEnd   = false; // Is the derivative at the end fixed in fit?
		double der_start;         // Derivative at the start of the spline.
		double der_end;           // Derivative at the end of the spline.

	public:
		/// @brief Constructor for NSpline.
		/// @param nodes_x An array of N doubles representing the x-coordinates of the nodes.
		NSpline(const double nodes_x[N])
		{
			for (int i = 0; i < N; ++i) this->nodes_x[i] = nodes_x[i];
		}

		/// @brief Constructor for NSpline that generates a set of evenly spaced nodes between min and max.
		/// @param min The minimum x-value of the spline.
		/// @param max The maximum x-value of the spline.
		NSpline(double min, double max)
		{
			for (int i = 0; i < N; i++) this->nodes_x[i] = min+(i/(N-1.0))*(max-min);
		}
		
		/// @brief Set the starting derivative for the spline and fix it in fit.
		/// @param der_start The value of the starting derivative.
		void SetDerStart(double der_start) {this->der_start = der_start; FixDerStart = true;}

		/// @brief Set the ending derivative for the spline and fix it in fit.
		/// @param der_start The value of the ending derivative.
		void SetDerEnd(double der_end) {this->der_end = der_end; FixDerEnd = true;}

		/// @brief Evaluates the NSpline at a given x-value.
		/// @param x A pointer to an array of doubles representing the x-value at which to evaluate the spline.
		/// @param par A pointer to an array of doubles representing the parameters of the spline.
		/// @return The value of the NSpline at the given x-value.
		double Eval(double* x, double* par)
		{
			/*Fit parameters:
			par[0-N-1]=X of nodes (to be fixed in the fit!)
			par[N-2N-1]=Y of nodes
			par[2N-2N+1]=first derivative at begin and end (to be fixed in the fit!)
			*/

			for (int i = 0; i < N; ++i) par[i] = nodes_x[i];

			Double_t xx = x[0];

			Double_t xn[N];
			Double_t yn[N];

			for (int i = 0; i < N; ++i) xn[i] = par[i];
			for (int i = 0; i < N; ++i) yn[i] = par[N+i];

			if(FixDerStart) par[2*N] = der_start;
			if(FixDerEnd) par[2*N+1] = der_end;

			Double_t b1 = par[2*N];
			Double_t e1 = par[2*N+1];

			TSpline3 sp3("sp3", xn, yn, N, "b1e1", b1, e1);

			return sp3.Eval(xx);		
		}

		/// @brief Returns a std::function object that can be used to evaluate the NSpline.
		/// @return A std::function object that takes two pointers to arrays of doubles and returns a double.
		EvalFn GetEval()
		{
			return [this](double* a, double* b)
			{
				return this->Eval(a,b);
			};
		}

	};

	/// @brief Fits cubic splines to a TGraph object.
	/// @tparam nodes The number of nodes for the spline interpolation.
	/// @param graph The TGraph object to fit splines to.
	/// @param min The minimum x-coordinate value for the spline.
	/// @param max The maximum x-coordinate value for the spline.
	/// @return A pointer to the TSpline3 object representing the fitted spline.
	template<int nodes>
	TSpline3* FitSplines(TGraph* graph, const double& min, const double& max)
	{
		double nodes_x[nodes];
		for (int i = 0; i < nodes; i++) nodes_x[i] = min+(i/(nodes-1.0))*(max-min);

		NSpline<nodes>* s = new NSpline<nodes>(nodes_x);

		TF1* fit = new TF1("fit",s->GetEval(),min,max,2*nodes+2);
		graph->Fit(fit,"M","",min,max);

		double* fit_pars = fit->GetParameters();
		double nodes_y[nodes];
		for(int i = nodes; i < 2*nodes; i++) nodes_y[i-nodes]=fit_pars[i];
		double beg1 = fit_pars[2*nodes];
		double end1 = fit_pars[2*nodes+1];
		
		return new TSpline3("sp3",nodes_x,nodes_y,nodes,"b1e1",beg1,end1);
	}
} // namespace X17