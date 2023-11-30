// C++ dependencies
#include <filesystem>
#include <iostream>
#include <regex>
#include <string>

// ROOT dependencies
#include "TChain.h"
#include "TLine.h"

// X17 dependencies
#include "Utilities.h"

std::string GetNextFilePath(std::string folder_path, std::string prefix, std::string suffix)
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
                int file_index = std::stoi(match[1].str());                          // Extract the file index from the matched group.
                if (file_index >= next_file_index) next_file_index = file_index + 1; // Update the next file index if a higher index is found.
            }
        }
    }

    std::string next_filename = prefix + std::to_string(next_file_index) + suffix; // Generate the next file name.
    return folder_path + "/" + next_filename;
}

int sign(double x) 
{
    if (x > 0)      return  1;
    else if (x < 0) return -1;
    else            return  0;
}

void AddFilesToTChain(TChain* chain, std::string prefix, std::string suffix, int start, int end)
{
    for (int i = start; i <= end; ++i) 
    {
        std::string filename = prefix + std::to_string(i) + suffix;
        
        // Check if the file exists, if it doesn't print a warning.
        if (std::filesystem::exists(filename))
        {
            chain->Add(filename.c_str());
        } 
        else std::cerr << "Warning: File " << filename << " does not exist. Skipping." << std::endl;
    }
}

double GetFWHM(TH1F* h, bool draw)
{
    double max = h->GetMaximum();
    int left_bin   = h->FindFirstBinAbove(max/2.0);
    int right_bin  = h->FindLastBinAbove(max/2.0);

    double left_slope  = (h->GetBinContent(left_bin)    - h->GetBinContent(left_bin-1)) / h->GetBinWidth(1);
    double right_slope = (h->GetBinContent(right_bin+1) - h->GetBinContent(right_bin))  / h->GetBinWidth(1);

    double left_delta  =  (h->GetBinContent(left_bin)  - (max/2)) / left_slope;
    double right_delta = -(h->GetBinContent(right_bin) - (max/2)) / right_slope;

    double left  = h->GetBinCenter(left_bin)  - left_delta;
    double right = h->GetBinCenter(right_bin) + right_delta;

    if(draw)
    {
        TLine* l = new TLine(left,max/2,right,max/2);
        l->SetLineColor(kRed);
        l->Draw();
    }

    return right - left;
}