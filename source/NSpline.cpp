// ROOT dependencies
#include "TSpline.h"

// X17 dependencies
#include "Field.h"
#include "Vector.h"
#include "X17Utilities.h"

namespace X17
{
    void RecoEnergy(TSpline3* sp_fit, const Field<Vector>& magfield, TGraph* energy, TGraph* radius, TGraph* magnetic, double min, double max, double step)
    {
        using namespace constants;

        for (double x = min; x <= max; x += step)
        {
            int i_node = sp_fit->FindX(x);

            double xnode, ynode, b, c, d;
            sp_fit->GetCoeff(i_node, xnode, ynode, b, c, d);

            // Calculate the radius [m] of the track from the fitted spline.
            double dx   = x - xnode;
            double der  = b + dx * (2 * c + 3 * d * dx);
            double der2 = 2 * c + 6 * d * dx;
            double r    = cm2m * pow(1 + der * der, 1.5) / der2;
            
            Vector B      = magfield.GetField(x, 0, 8 - sp_fit->Eval(x));
            double betasq = 1 / (1 + pow((E0 / (c * r * B.x)), 2));
            double Ekin   = E0 * (1 / std::sqrt(1 - betasq) - 1);

            if (x > 4)
                energy->AddPoint(x, Ekin / 1e6);
            if (r < 1 && r > 0)
                radius->AddPoint(x, r * m2cm);

            magnetic->AddPoint(x, B.y);
        }
    }
}