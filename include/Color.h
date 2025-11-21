#pragma once

// C++ dependencies
#include <iostream>

// ROOT dependencies
#include "TColor.h"
#include "TROOT.h"

namespace X17::Color
{
    /// @brief CIELAB color space.
    struct LAB
    {
        double L, a, b;
    };

    /// @brief sRGB to linear RGB conversion.
    /// @param c Color value (0-1).
    /// @return Linear RGB value.
    inline double inverse_gamma(double c)
    {
        return (c <= 0.04045) ? c / 12.92 : std::pow((c + 0.055) / 1.055, 2.4);
    }

    // RGB (NDC: 0–1) → XYZ

    /// @brief Converts RGB to CIE XYZ (D65).
    /// @param r Red component (0-1).
    /// @param g Green component (0-1).
    /// @param b Blue component (0-1).
    /// @param X Output X component.
    /// @param Y Output Y component.
    /// @param Z Output Z component.
    inline void RGBtoXYZ(float r, float g, float b, double& X, double& Y, double& Z)
    {
        r = inverse_gamma(r);
        g = inverse_gamma(g);
        b = inverse_gamma(b);

        X = r * 0.4124 + g * 0.3576 + b * 0.1805;
        Y = r * 0.2126 + g * 0.7152 + b * 0.0722;
        Z = r * 0.0193 + g * 0.1192 + b * 0.9505;
    }

    /// @brief Converts CIE XYZ to CIELAB.
    /// @param X X component.
    /// @param Y Y component.
    /// @param Z Z component.
    /// @return CIELAB color.
    inline LAB XYZtoLAB(double X, double Y, double Z)
    {
        const double Xn = 0.95047, Yn = 1.00000, Zn = 1.08883;
        auto f = [](double t) {
            return (t > 0.008856) ? std::cbrt(t) : (7.787 * t + 16.0 / 116.0);
        };
        double fx = f(X / Xn), fy = f(Y / Yn), fz = f(Z / Zn);
        return {116.0 * fy - 16.0, 500.0 * (fx - fy), 200.0 * (fy - fz)};
    }

    /// @brief Converts RGB to CIELAB.
    /// @param r Red component (0-1).
    /// @param g Green component (0-1).
    /// @param b Blue component (0-1).
    /// @return CIELAB color.
    inline LAB RGBtoLAB(float r, float g, float b)
    {
        double X, Y, Z;
        RGBtoXYZ(r, g, b, X, Y, Z);
        return XYZtoLAB(X, Y, Z);
    }

    /// @brief ΔE76 distance between two LAB colors.
    /// @param c1 First color.
    /// @param c2 Second color.
    /// @return ΔE76 distance.
    inline double DeltaE76(const LAB& c1, const LAB& c2)
    {
        return std::sqrt(
            std::pow(c1.L - c2.L, 2) +
            std::pow(c1.a - c2.a, 2) +
            std::pow(c1.b - c2.b, 2)
        );
    }

    /// @brief Returns ΔE76 distance between an RGB color and a ROOT color.
    /// @param r Red component (0-1).
    /// @param g Green component (0-1).
    /// @param b Blue component (0-1).
    /// @param col ROOT color.
    /// @return ΔE76 distance.
    inline double ColorDistanceCIELAB(float r, float g, float b, const TColor* col)
    {
        LAB ref = RGBtoLAB(r, g, b);
        LAB tgt = RGBtoLAB(col->GetRed(), col->GetGreen(), col->GetBlue());
        return DeltaE76(ref, tgt);
    }

    /// @brief Creates a new RGB color and returns its index.
    /// @param R Red value (0-255).
    /// @param G Green value (0-255).
    /// @param B Blue value (0-255).
    /// @param saturation Saturation value (0-100), defaults to 100.
    /// @param find_closest If true, the closest predefined color will be found and returned, defaults to true.
    /// @return Index of the created ROOT color.
    inline int RGB(int R, int G, int B, float saturation = 100.0, bool find_closest = true)
    {
        float red_NDC   = 1 + (R/255.0 - 1) * saturation/100.0;
        float green_NDC = 1 + (G/255.0 - 1) * saturation/100.0;
        float blue_NDC  = 1 + (B/255.0 - 1) * saturation/100.0; 

        if(find_closest)
        {
            double closest_dist = std::numeric_limits<double>::max();
            int closest_index = -1;
            for (int i = 0; i < gROOT->GetListOfColors()->GetSize(); i++)
            {
                TColor* color = (TColor*)gROOT->GetListOfColors()->At(i);
                if(!color) continue;
                double dist = ColorDistanceCIELAB(red_NDC, green_NDC, blue_NDC, color);
                if (dist < closest_dist)
                {
                    closest_dist = dist;
                    closest_index = color->GetNumber();
                }
            }
            
            // std::cout << "X17::Color::RGB(): Closest color index: " << closest_index << std::endl;
            return closest_index;
        }
        else
        {
            int color_index = TColor::GetColor(red_NDC, green_NDC, blue_NDC);
            // std::cout << "X17::Color::RGB(): Color index: " << color_index << std::endl;
            return color_index;
        }
    }

    /// @brief Creates a new orange color (from Garfield) and returns its index.
    /// @return Index of the created ROOT color.
    inline int GOrange() { return RGB(255, 153, 0); }
}