#pragma once

// C++ dependencies
#include <string>
#include <vector>

// ROOT dependencies
#include "TH1F.h"
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

/// @brief Returns random number between given minimal and maximal value.
/// @param rand The TRandom3 instance that should be used to generate the number.
/// @param min The minimal value.
/// @param max The maximal value.
/// @return A random number between given minimal and maximal value.
inline double RandomMinMax(TRandom3* rand, double min, double max) { return min + (max - min) * rand->Rndm(); }

/// @brief Function for finding next available filename in given folder (such as name1.root).
/// @param folder_path Folder where the search should happen.
/// @param prefix Name prefix of the files.
/// @param suffix Suffix of the file (default is .root).
/// @return Next available filename with given patern.
std::string GetNextFilePath(std::string folder_path, std::string prefix, std::string suffix = ".root");

/// @brief Function for determining the sign of a number.
/// @param x The number sign of which will be determined.
/// @return Returns 1 for positive, 0 for zero and -1 for negative x.
int sign(double x);

/// @brief Function for adding multiple files (with enumerated names) to a TChain. Example file name file1.root
/// @param chain TChain to which the files will be added.
/// @param prefix The prefix part of the filename (file1.root --> file).
/// @param suffix The suffix part of the filename (file1.root --> .root).
/// @param start The first index of the filename.
/// @param end The last index of the filename.
void AddFilesToTChain(TChain* chain, std::string prefix, std::string suffix, int start, int end);

/// @brief Uses the histogram maximum and linear interpolation to approximate FWHM.
/// @param histogram The histogram to be used.
/// @param draw Should the result get drawn as a line?
/// @return The FWHM approximation of the histogram.
double GetFWHM(TH1F* histogram, bool draw);

/// @brief Calculates the bias factor for the common normal distribution standard deviation estimator.
/// @param N Number of samples.
/// @return The bias factor.
double StdevBiasFactor(int N);

/// @brief Calculates the average value of a vector.
/// @tparam T The type of the values in the vector.
/// @param values The vector of values.
/// @return The average value of the vector.
template <typename T>
T GetAverage(const std::vector<T>& values)
{
    T sum = 0;
    for (T value : values) sum += value;
    return sum/values.size();
}

/// @brief Calculates the linearly interpolated quantile of a vector.
/// @tparam T The type of the values in the vector.
/// @param values The vector of values.
/// @param quantile The quantile to calculate.
/// @param sorted Is the vector sorted?
/// @return The linearly interpolated quantile of the vector.
template <typename T>
T GetQuantile(std::vector<T>& values, double quantile, bool sorted = false)
{
    if (!sorted) std::sort(values.begin(),values.end());
    double index = values.size()*quantile;
    int floor_index = std::floor(index);
    double mantissa = index - floor_index;
    return (1-mantissa)*values[floor_index] + mantissa*values[floor_index+1];
}