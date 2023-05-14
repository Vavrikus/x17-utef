#pragma once

// ROOT dependencies
#include "TSpline.h"

namespace X17
{
	/// @brief Represents a cubic spline that interpolates a set of N nodes. Defined by the positions of the nodes, derivatives at the start and end of the spline.
	/// @tparam N The number of nodes.
	template<int N>
	class NSpline 
	{
		typedef std::function<double(double*,double*)> EvalFn;

	private:
		double m_nodes_x[N];          // The x-coordinates of spline nodes.
		bool m_fix_der_start = false; // Is the derivative at the start fixed in fit?
		bool m_fix_der_end   = false; // Is the derivative at the end fixed in fit?
		double m_der_start;           // Derivative at the start of the spline.
		double m_der_end;             // Derivative at the end of the spline.

	public:
		/// @brief Constructor for NSpline.
		/// @param nodes_x An array of N doubles representing the x-coordinates of the nodes.
		NSpline(const double nodes_x[N])
		{
			for (int i = 0; i < N; ++i) this->m_nodes_x[i] = nodes_x[i];
		}

		/// @brief Constructor for NSpline that generates a set of evenly spaced nodes between min and max.
		/// @param min The minimum x-value of the spline.
		/// @param max The maximum x-value of the spline.
		NSpline(double min, double max)
		{
			for (int i = 0; i < N; i++) this->m_nodes_x[i] = min + (i / (N - 1.0)) * (max - min);
		}
		
		/// @brief Set the starting derivative for the spline and fix it in fit.
		/// @param der_start The value of the starting derivative.
		void SetDerStart(double der_start) {this->m_der_start = der_start; m_fix_der_start = true;}

		/// @brief Set the ending derivative for the spline and fix it in fit.
		/// @param der_start The value of the ending derivative.
		void SetDerEnd(double der_end) {this->m_der_end = der_end; m_fix_der_end = true;}

		/// @brief Evaluates the NSpline at a given x-value.
		/// @param x A pointer to an array of doubles representing the x-value at which to evaluate the spline.
		/// @param par A pointer to an array of doubles representing the parameters of the spline.
		/// @return The value of the NSpline at the given x-value.
		double Eval(double* x, double* par) const;

		/// @brief Returns a std::function object that can be used to evaluate the NSpline.
		/// @return A std::function object that takes two pointers to arrays of doubles and returns a double.
		EvalFn GetEval() const;
	};

	/// @brief Fits cubic splines to a TGraph object.
	/// @tparam nodes The number of nodes for the spline interpolation.
	/// @param graph The TGraph object to fit splines to.
	/// @param min The minimum x-coordinate value for the spline.
	/// @param max The maximum x-coordinate value for the spline.
	/// @return A pointer to the TSpline3 object representing the fitted spline.
	template<int nodes>
	TSpline3* FitSplines(TGraph* graph, double min, double max);
} // namespace X17

// Templated function definitions.
#include "NSpline.inl"