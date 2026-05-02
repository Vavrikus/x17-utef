// C++ dependencies
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <regex>
#include <string>
#include <vector>

// ROOT dependencies
#include "Rtypes.h"
#include "TChain.h"
#include "TF1.h"
#include "TGraph2D.h"
#include "TH1.h"
#include "TLine.h"
#include "TPolyLine3D.h"

// X17 dependencies
#include "Utilities.h"
#include "Vector.h"

std::string GetNextFilePath(const std::string& folder_path, const std::string& prefix, const std::string& suffix)
{
  std::string regex_patern = prefix + "(\\d+)\\" + suffix;
  std::regex filePattern(regex_patern); // Regular expression to match the file pattern.
  int next_file_index = 1;              // Start with the default next file index of 1.

  // Iterate over files in the folder.
  for (const auto& entry : std::filesystem::directory_iterator(folder_path))
  {
    // Check if it's a regular file.
    if (std::filesystem::is_regular_file(entry.path()))
    {
      std::smatch match;
      std::string file_name = entry.path().filename().string(); // Get the file name.

      // Check if the file name matches the pattern.
      if (std::regex_match(file_name, match, filePattern))
      {
        int file_index = std::stoi(match[1].str()); // Extract the file index from the matched group.
        if (file_index >= next_file_index)
          next_file_index = file_index + 1; // Update the next file index if a higher index is found.
      }
    }
  }

  std::string next_filename = prefix + std::to_string(next_file_index) + suffix; // Generate the next file name.
  return folder_path + "/" + next_filename;
}

int sign(double x)
{
  if (x > 0)
    return 1;
  if (x < 0)
    return -1;
  return 0;
}

void AddFilesToTChain(TChain* chain, const std::string& prefix, const std::string& suffix, int start, int end)
{
  for (int i = start; i <= end; ++i)
  {
    std::string filename = prefix + std::to_string(i);
    filename += suffix;

    // Check if the file exists, if it doesn't print a warning.
    if (std::filesystem::exists(filename))
    {
      chain->Add(filename.c_str());
    }
    else
      std::cerr << "Warning: File " << filename << " does not exist. Skipping." << '\n';
  }
}

double GetFWHM(TH1F* h, bool draw)
{
  double max    = h->GetMaximum();
  int left_bin  = h->FindFirstBinAbove(max / 2.0);
  int right_bin = h->FindLastBinAbove(max / 2.0);

  double left_slope  = (h->GetBinContent(left_bin) - h->GetBinContent(left_bin - 1)) / h->GetBinWidth(1);
  double right_slope = (h->GetBinContent(right_bin + 1) - h->GetBinContent(right_bin)) / h->GetBinWidth(1);

  double left_delta  = (h->GetBinContent(left_bin) - (max / 2)) / left_slope;
  double right_delta = -(h->GetBinContent(right_bin) - (max / 2)) / right_slope;

  double left  = h->GetBinCenter(left_bin) - left_delta;
  double right = h->GetBinCenter(right_bin) + right_delta;

  if (draw)
  {
    TLine* l = new TLine(left, max / 2, right, max / 2);
    l->SetLineColor(kRed);
    l->Draw();
  }

  return right - left;
}

double StdevBiasFactor(int N)
{
  // Using logarithmic version of the gamma function to avoid overflow
  double log_g1 = std::lgamma(0.5 * N);
  double log_g2 = std::lgamma(0.5 * (N - 1));
  double root   = std::sqrt(2.0 / (N - 1));

  double val = 1.0 / (root * std::exp(log_g1 - log_g2));
  return val;
}

X17::Vector GetStDev(const std::vector<X17::Vector>& values, double* mag_sigma)
{
  X17::Vector average = GetAverage(values);
  X17::Vector sum_sq(0, 0, 0);
  for (X17::Vector value : values)
  {
    X17::Vector diff = value - average;
    sum_sq += X17::Vector(diff.x * diff.x, diff.y * diff.y, diff.z * diff.z);
  }
  sum_sq /= static_cast<double>(values.size() - 1);

  double mag = std::sqrt(sum_sq.x + sum_sq.y + sum_sq.z);
  *mag_sigma = mag * StdevBiasFactor(static_cast<int>(values.size()));

  X17::Vector stdev(std::sqrt(sum_sq.x), std::sqrt(sum_sq.y), std::sqrt(sum_sq.z));
  return stdev * StdevBiasFactor(static_cast<int>(values.size()));
}

double GetQuantile(std::vector<double>& values, double quantile, bool sorted)
{
  if (!sorted)
    std::ranges::sort(values);
  double index    = static_cast<double>(values.size()) * quantile;
  int floor_index = std::floor(index);
  double mantissa = index - floor_index;
  return (1 - mantissa) * values[floor_index] + mantissa * values[floor_index + 1];
}

double GetPvalue(std::vector<double>& values, double value, bool sorted)
{
  if (!sorted)
    std::ranges::sort(values);

  // Handle edge cases
  if (value < values.front())
    return 0.0; // Extremely low  (p ≈ 0)
  if (value > values.back())
    return 0.0; // Extremely high (p ≈ 0)

  double vals_size = static_cast<double>(values.size());

  // Fraction of values ≤ input (left tail)
  double left_p = static_cast<double>(std::ranges::lower_bound(values, value) - values.begin()) / vals_size;

  // Fraction of values ≥ input (right tail)
  double right_p = static_cast<double>(values.end() - std::ranges::upper_bound(values, value)) / vals_size;

  // Two-tailed p-value: min(left, right) × 2
  return 2.0 * std::min(left_p, right_p);
}

void ReportProgress(int current, int total, int width)
{
  current++;
  float progress = static_cast<float>(current) / static_cast<float>(total);
  int filled     = ifloor(progress * static_cast<float>(width));

  std::cout << "\r["; // Carriage return to overwrite the line
  for (int i = 0; i < width; ++i)
    std::cout << (i < filled ? '=' : ' ');

  std::cout << "] " << std::setw(3) << ifloor(progress * 100) << "%";
  std::cout.flush();

  if (current == total)
    std::cout << '\n';
}

TPolyLine3D* GetLine3D(TGraph2D* graph)
{
  if (graph->GetN() > 1)
  {
    TPolyLine3D* line = new TPolyLine3D(graph->GetN());
    for (int i = 0; i < graph->GetN(); ++i)
      line->SetPoint(i, graph->GetX()[i], graph->GetY()[i], graph->GetZ()[i]);
    return line;
  }
  return nullptr;
}

void GetSkewGausStats(TF1* skew_gaus, double& mean, double& sigma, double& skew, double& fwhm)
{
  double p1 = skew_gaus->GetParameter(1);
  double p2 = skew_gaus->GetParameter(2);
  double p3 = skew_gaus->GetParameter(3);

  // Actual moments differ: (see wikipedia Skew normal distribution)
  double delta = p3 / (std::sqrt(1 + p3 * p3));
  mean         = p1 + delta * p2 * std::sqrt(2 / std::numbers::pi);
  sigma        = p2 * std::sqrt(1 - delta * delta * 2 / std::numbers::pi);
  skew         = (4 - std::numbers::pi) / 2 * std::pow((delta * std::sqrt(2 / std::numbers::pi)), 3)
                 / std::pow((1 - delta * delta * 2 / std::numbers::pi), 1.5);

  // Numerical FWHM
  double x_max = skew_gaus->GetMaximumX();
  double y_max = skew_gaus->GetMaximum();
  double x1    = skew_gaus->GetX(y_max / 2., skew_gaus->GetXmin(), x_max); // left half
  double x2    = skew_gaus->GetX(y_max / 2., x_max, skew_gaus->GetXmax()); // right half
  fwhm         = x2 - x1;
}
