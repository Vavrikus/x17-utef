// X17 dependencies
#include "RK4.h"

namespace X17
{
    template<int N>        
    RK4<N>::RK4(double start, double step, VecFn equation, VectorN initial, EndFn end_condition)
        : m_param(start), m_step(step), m_dif_eq(equation), m_end_fn(end_condition), m_current(initial)
    {
        results.push_back(m_current);
    }

    template<int N>
    void RK4<N>::Integrate()
    {
        while(!m_end_fn(m_step, m_current))
        {
            VectorN k1,k2,k3,k4;
            m_dif_eq(m_param,          m_current,      k1);
            k1 *= m_step;
            m_dif_eq(m_param+m_step/2, m_current+k1/2, k2);
            k2 *= m_step;
            m_dif_eq(m_param+m_step/2, m_current+k2/2, k3);
            k3 *= m_step;
            m_dif_eq(m_param+m_step,   m_current+k3,   k4);
            k4 *= m_step;

            VectorN diff = (1.0/6.0) * (k1 + 2*k2 + 2*k3 + k4);
            m_current += diff;
            m_param   += m_step;
            results.push_back(m_current);
        }
    }
} // namespace X17