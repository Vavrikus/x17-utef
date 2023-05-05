#pragma once

#include <cmath>
#include <vector>

#include "TRandom3.h"

/// @brief Find the minimum value among the given parameters.
/// @tparam T The type of the parameters, must be a numeric type.
/// @param args The parameters to find the minimum value from.
/// @return The minimum value among the parameters.
template<typename... T>
double find_min(T... args) 
{
    double values[] = { args... };
    return *std::min_element(values, values + sizeof...(args));
}

/// @brief Calculates the standard deviation of a vector of double values assuming a Gaussian distribution.
/// @tparam T Numeric type that has operators += and /. Must be possible to cast it to zero.
/// @param values The vector of double values.
/// @param average The average of the vector of double values.
/// @return The standard deviation of the vector of double values.
template<typename T>
T stdev(std::vector<T> values, T average)
{
    T sqdev_sum = static_cast<T>(0);
    for (T d : values) sqdev_sum += pow((d - average), 2);

    return sqrt(sqdev_sum / values.size());
}

/// @brief Returns random number between given minimal and maximal value.
/// @param rand The TRandom3 instance that should be used to generate the number.
/// @param min The minimal value.
/// @param max The maximal value.
/// @return A random number between given minimal and maximal value.
double RandomMinMax(TRandom3* rand, double min, double max) { return min + (max-min)*rand->Rndm(); }