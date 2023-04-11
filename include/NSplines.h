#pragma once
#include <TSpline.h>

template<int N>
class NSpline 
{
	typedef std::function<double(double*,double*)> EvalFn;

private:
	double node_pos[N];
	bool FixDer1 = false;
	bool FixDer2 = false;
	double der1,der2;

public:
	NSpline(const double node_pos[N])
	{
		for (int i = 0; i < N; ++i) this->node_pos[i] = node_pos[i];
	}

    NSpline(double min, double max)
    {
        for (int i = 0; i < N; i++) this->node_pos[i] = min+(i/(N-1.0))*(max-min);
    }

	void SetDer1(double der1) {this->der1 = der1; FixDer1 = true;}
	void SetDer2(double der2) {this->der2 = der2; FixDer2 = true;}

	double Eval(double* x, double* par)
	{
		/*Fit parameters:
		par[0-N-1]=X of nodes (to be fixed in the fit!)
		par[N-2N-1]=Y of nodes
		par[2N-2N+1]=first derivative at begin and end (to be fixed in the fit!)
		*/

		for (int i = 0; i < N; ++i) par[i] = node_pos[i];

		Double_t xx = x[0];

		Double_t xn[N];
		Double_t yn[N];

		for (int i = 0; i < N; ++i) xn[i] = par[i];
		for (int i = 0; i < N; ++i) yn[i] = par[N+i];

		if(FixDer1) par[2*N]   = der1;
		if(FixDer2) par[2*N+1] = der2;

		Double_t b1 = par[2*N];
		Double_t e1 = par[2*N+1];

		TSpline3 sp3("sp3", xn, yn, N, "b1e1", b1, e1);

		return sp3.Eval(xx);		
	}

	EvalFn GetEval()
	{
		return [this](double* a, double* b)
		{
			return this->Eval(a,b);
		};
	}

};

template<int nodes>
TSpline3* FitSplines(TGraph* graph, const double& min, const double& max)
{
    double node_pos[nodes];
    for (int i = 0; i < nodes; i++) node_pos[i] = min+(i/(nodes-1.0))*(max-min);

    NSpline<nodes>* s = new NSpline<nodes>(node_pos);

    TF1* fit = new TF1("fit",s->GetEval(),min,max,2*nodes+2);
    graph->Fit(fit,"M","",min,max);

    double* fit_pars = fit->GetParameters();
    double nodes_y[nodes];
    for(int i = nodes; i < 2*nodes; i++) nodes_y[i-nodes]=fit_pars[i];
    double beg1 = fit_pars[2*nodes];
    double end1 = fit_pars[2*nodes+1];
    
    return new TSpline3("sp3",node_pos,nodes_y,nodes,"b1e1",beg1,end1);
}