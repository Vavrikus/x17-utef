// X17 dependencies
#include "NSpline.h"

namespace X17
{
	//// Public methods.

	template<int N>
	double NSpline<N>::Eval(double* x, double* par)
	{
		/*Fit parameters:
		par[0-N-1]=X of nodes (to be fixed in the fit!)
		par[N-2N-1]=Y of nodes
		par[2N-2N+1]=first derivative at begin and end (to be fixed in the fit!)
		*/

		for (int i = 0; i < N; ++i) par[i] = m_nodes_x[i];

		Double_t xx = x[0];

		Double_t xn[N];
		Double_t yn[N];

		for (int i = 0; i < N; ++i) xn[i] = par[i];
		for (int i = 0; i < N; ++i) yn[i] = par[N + i];

		if(m_fix_der_start) par[2 * N] = m_der_start;
		if(m_fix_der_end) par[2 * N + 1] = m_der_end;

		Double_t b1 = par[2 * N];
		Double_t e1 = par[2 * N + 1];

		TSpline3 sp3("sp3", xn, yn, N, "b1e1", b1, e1);

		return sp3.Eval(xx);		
	}

	template<int N>
	typename NSpline<N>::EvalFn NSpline<N>::GetEval()
	{
		return [this](double* a, double* b)
		{
			return this->Eval(a,b);
		};
	}





	//// Functions related to the NSpline class.

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