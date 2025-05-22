#pragma once

// C++ dependencies
#include <iomanip>
#include <string>
#include <vector>

// ROOT dependencies
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
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

/// @brief Calculates the bias factor (1/c_4(N)) for the common normal distribution standard deviation estimator.
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
/// @param values The vector of values.
/// @param quantile The quantile to calculate.
/// @param sorted Is the vector sorted?
/// @return The linearly interpolated quantile of the vector.
double GetQuantile(std::vector<double>& values, double quantile, bool sorted = false);

/// @brief Calculates the symmetric p-value of a value in a vector.
/// @param values The vector of values.
/// @param value The value to calculate the p-value for.
/// @param sorted Is the vector sorted?
/// @return The p-value of the value in the vector.
double GetPvalue(std::vector<double>& values, double value, bool sorted = false);

/// @brief Prints a progress bar to the console.
/// @param current The current value of the progress.
/// @param total The total value of the progress.
/// @param width The width of the progress bar.
void ReportProgress(int current, int total, int width = 50);

/// @brief Checks if a type is a member of a parameter pack.
/// @tparam T The type to check.
/// @tparam ...Ts The types in the parameter pack.
template <typename T, typename... Ts>
constexpr bool is_from = (std::is_same_v<T, Ts> || ...);

template <typename T>
void ApplyThesisStyle(T* obj)
{
    if constexpr (is_from<T, TCanvas, TPad>)
    {
        double w = obj->GetWw();
        if (w > 800)
        {
            obj->SetLeftMargin(0.1);
            obj->SetRightMargin(0.05);
        }
        else
        {
            obj->SetLeftMargin(0.13);
            obj->SetRightMargin(0.07);
        }
        obj->SetTopMargin(0.07);
        obj->SetBottomMargin(0.13);
    }

    else
    {
        std::vector<TAxis*> axes;
        if constexpr (is_from<T, TH1F, TGraph, TGraphErrors>)
            axes = { obj->GetXaxis(), obj->GetYaxis(), nullptr };

        else if constexpr (is_from<T, TH2F, TGraph2D>)
            axes = { obj->GetXaxis(), obj->GetYaxis(), obj->GetZaxis() };
        
        else
        {
            std::cout << "Warning: Can't apply thesis style for object of type " << typeid(obj).name() << std::endl;
            return;
        }

        for (TAxis* axis : axes) if (axis)
        {
            axis->SetTitleSize(0.06);
            axis->SetLabelSize(0.055);
        }

        if (axes[2])
        {
            axes[0]->SetTitleOffset(1.2);
            axes[1]->SetTitleOffset(0.9);
            axes[2]->SetTitleOffset(1.3);
        }
        else
        {
            axes[0]->SetTitleOffset(0.9);
            axes[1]->SetTitleOffset(0.75);
        }
        
    }
}