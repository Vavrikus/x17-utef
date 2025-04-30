#pragma once

// C++ dependencies
#include <cmath>

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
        double Angle(Vector other) const
        {
            return acos((x * other.x + y * other.y + z * other.z) / (this->Magnitude() * other.Magnitude()));
        }

        /// @brief Compute the square of the distance between this vector and another vector.
        /// @param[in] other The other vector.
        /// @return The square of the distance.
        double SqDist(Vector other) const { return pow(x - other.x, 2) + pow(y - other.y, 2) + pow(z - other.z, 2); }

        /// @brief Add another vector to this vector.
        /// @param[in] v The vector to add.
        void operator+=(Vector v);

        /// @brief Divide each component of the vector by a scalar.
        /// @param d The scalar to divide by.
        void operator/=(double d);

        /// @brief Equality operator for comparing two Vector objects. All components must be equal.
        /// @param v The vector to compare with.
        /// @return True if the vectors are equal, false otherwise.
        bool operator==(Vector v) const { return (x == v.x) && (y == v.y) && (z == v.z); }

        /// @brief Adds two vectors component-wise and returns the result.
        /// @param v2 The second vector.
        /// @return The sum of the two vectors.
        Vector operator+(Vector v2) const
        {
            return Vector{x + v2.x, y + v2.y, z + v2.z};
        }

        /// @brief Subtracts two vectors component-wise.
        /// @param v2 The second vector.
        /// @return Vector The difference between the two vectors.
        Vector operator-(Vector v2) const
        {
            return Vector{x - v2.x, y - v2.y, z - v2.z};
        }

        /// @brief Scalar multiplication of a vector.
        /// @param d The scalar to multiply with.
        /// @return The scaled vector.
        Vector operator*(double d) const
        {
            return Vector{d * x, d * y, d * z};
        }

        /// @brief Divide a vector by a scalar.
        /// @param d The scalar to divide by (cannot be zero).
        /// @return The resulting vector.
        /// @throw std::invalid_argument if d is zero (only if DEBUG defined).
        Vector operator/(double d) const
        {
        #ifdef DEBUG
            if (d == 0) throw std::invalid_argument("division by zero");
        #endif

            return Vector{x / d, y / d, z / d};
        }
        
        /// @brief Compute the dot product of two vectors.
        /// @param v2 The second vector.
        /// @return The dot product.
        double operator*(Vector v2) const
        {
            return x * v2.x + y * v2.y + z * v2.z;
        }
    };

    /// @brief Scalar multiplication of a vector.
    /// @param d The scalar to multiply with.
    /// @param v The vector to multiply.
    /// @return The scaled vector.
    inline Vector operator*(double d, Vector v)
    {
        return Vector{d * v.x, d * v.y, d * v.z};
    }

    /// @brief Normalize a vector. Divides the vector by its magnitude to produce a unit vector in the same direction.
    inline void Vector::Normalize() {*this = (*this) / (this->Magnitude());}

    /// @brief Computes the squared distance between a line and a point in 3D space.
    /// @param origin The origin of the line.
    /// @param orientation The orientation vector of the line.
    /// @param max_param The maximum allowed parameter of the line.
    /// @param point The point to compute the distance from.
    /// @return The squared distance between the line and the point.
    double LineSqDist(Vector origin, Vector orientation, double max_param, Vector point);

    /// @brief Computes the squared distance between a line and a point in 3D space and returns the closest point on the line.
    /// @param origin The origin of the line.
    /// @param orientation The orientation vector of the line.
    /// @param max_param The maximum allowed parameter of the line.
    /// @param point The point to compute the distance from.
    /// @param closest_point Will be set to the closest point on the line to the given point.
    /// @return The squared distance between the line and the point.
    double LineSqDistAndCP(Vector origin, Vector orientation, double max_param, Vector point, Vector& closest_point);
} // namespace X17