// C++ dependencies
#include <algorithm>

// X17 dependencies
#include "Vector.h"

namespace X17
{
    //// Public methods.

    void Vector::operator+=(Vector v)
    {
        this->x += v.x;
        this->y += v.y;
        this->z += v.z;
    }

    void Vector::operator-=(Vector v)
    {
        this->x -= v.x;
        this->y -= v.y;
        this->z -= v.z;
    }

    void Vector::operator*=(double d)
    {
        x *= d;
        y *= d;
        z *= d;
    }

    void Vector::operator/=(double d)
    {
        x /= d;
        y /= d;
        z /= d;
    }





    //// Functions related to the Vector class.

    double LineSqDist(Vector origin, Vector orientation, double max_param, Vector point)
    {
        // Normalize the orientation vector to ensure correct computations.
        Vector norm_orientation = orientation;
        norm_orientation.Normalize();
        
        // Compute the parameter t_close and clamp it to [0, max_param].
        double t_close = norm_orientation * (point - origin);
        t_close = std::clamp(t_close, 0.0, max_param);

        // Compute the vector from the line to the point and return its squared magnitude.
        Vector line_vector = origin + t_close * norm_orientation - point;
        return line_vector.SqMagnitude();
    }

    double LineSqDistAndCP(Vector origin, Vector orientation, double max_param, Vector point, Vector& closest_point)
    {
        // Normalize the orientation vector to ensure correct computations.
        Vector norm_orientation = orientation;
        norm_orientation.Normalize();
        
        // Compute the parameter t_close and clamp it to [0, max_param].
        double t_close = norm_orientation * (point - origin);
        t_close = std::clamp(t_close, 0.0, max_param);

        // Compute the vector from the line to the point and return its squared magnitude.
        closest_point = origin + t_close * norm_orientation;
        Vector line_vector = closest_point - point;
        return line_vector.SqMagnitude();
    }
} // namespace X17