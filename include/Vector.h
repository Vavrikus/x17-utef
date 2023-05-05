#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace X17
{
    /// @brief A 3D vector class with x, y, and z components.
    struct Vector
    {
        double x; // The x component of the vector.
        double y; // The y component of the vector.
        double z; // The z component of the vector.

        /// @brief Default constructor that initializes all components to 0.
        Vector() : x(0), y(0), z(0) { }

        /// @brief Constructor that takes three doubles to initialize the x, y, and z components.
        /// @param x The x component of the vector.
        /// @param y The y component of the vector.
        /// @param z The z component of the vector.
        Vector(double x, double y, double z) : x(x), y(y), z(z) { }

        /// @brief Add another vector to this vector.
        /// @param[in] v The vector to add.
        void operator+=(const Vector& v)
        {
            this->x += v.x;
            this->y += v.y;
            this->z += v.z;
        }
        /// @brief Divide each component of the vector by a scalar.
        /// @param d The scalar to divide by.
        void operator/=(const double& d)
        {
            x /= d;
            y /= d;
            z /= d;
        }

        /// @brief Equality operator for comparing two Vector objects. All components must be equal.
        /// @param v The vector to compare with.
        /// @return True if the vectors are equal, false otherwise.
        bool operator==(const Vector& v) const { return (x == v.x) && (y == v.y) && (z == v.z); }

        /// @brief Compute the square of the magnitude of this vector.
        /// @return The square of the magnitude.
        double SqMagnitude() const { return x * x + y * y + z * z; }

        /// @brief Compute the magnitude of this vector.
        /// @return The magnitude.
        double Magnitude() const { return sqrt(this->SqMagnitude()); }

        /// @brief Normalize this vector.
        void Normalize();

        /// @brief Compute the angle in radians between this vector and another vector.
        /// @param[in] other The other vector.
        /// @return The angle between the two vectors in radians.
        double Angle(const Vector& other) const
        {
            return acos((x * other.x + y * other.y + z * other.z) / (this->Magnitude() * other.Magnitude()));
        }

        /// @brief Compute the square of the distance between this vector and another vector.
        /// @param[in] other The other vector.
        /// @return The square of the distance.
        double SqDist(const Vector& other) const { return pow(x - other.x, 2) + pow(y - other.y, 2) + pow(z - other.z, 2); }
    };

    /// @brief Adds two vectors component-wise and returns the result.
    /// @param v1 The first vector.
    /// @param v2 The second vector.
    /// @return The sum of the two vectors.
    Vector operator+(const Vector& v1, const Vector& v2)
    {
        return Vector{v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
    }

    /// @brief Subtracts two vectors component-wise.
    /// @param v1 The first vector.
    /// @param v2 The second vector.
    /// @return Vector The difference between the two vectors.
    Vector operator-(const Vector& v1, const Vector& v2)
    {
        return Vector{v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
    }

    /// @brief Scalar multiplication of a vector.
    /// @param d The scalar to multiply with.
    /// @param v The vector to multiply.
    /// @return The scaled vector.
    Vector operator*(const double& d, const Vector& v)
    {
        return Vector{d * v.x, d * v.y, d * v.z};
    }

    /// @brief Scalar multiplication of a vector.
    /// @param v The vector to multiply.
    /// @param d The scalar to multiply with.
    /// @return The scaled vector.
    Vector operator*(const Vector& v, const double& d)
    {
        return Vector{d * v.x, d * v.y, d * v.z};
    }

    /// @brief Divide a vector by a scalar.
    /// @param v The vector to be divided.
    /// @param d The scalar to divide by (cannot be zero).
    /// @return The resulting vector.
    /// @throw std::invalid_argument if d is zero (only if DEBUG defined)
    Vector operator/(const Vector& v, const double& d)
    {
    #ifdef DEBUG
        if (d == 0) throw std::invalid_argument("division by zero");
    #endif

        return Vector{v.x/d, v.y/d, v.z/d};
    }
    
    /// @brief Compute the dot product of two vectors.
    /// @param v1 The first vector.
    /// @param v2 The second vector.
    /// @return The dot product.
    double operator*(const Vector& v1, const Vector& v2)
    {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }

    /// @brief Normalize a vector. Divides the vector by its magnitude to produce a unit vector in the same direction.
    void Vector::Normalize() {*this = (*this) / (this->Magnitude());}

    /// @brief Computes the squared distance between a line and a point in 3D space.
    /// @param origin The origin of the line.
    /// @param orientation The orientation vector of the line.
    /// @param max_param The maximum allowed parameter of the line.
    /// @param point The point to compute the distance from.
    /// @return The squared distance between the line and the point.
    double LineSqDist(const Vector& origin, const Vector& orientation, double max_param, const Vector& point)
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
} // namespace X17