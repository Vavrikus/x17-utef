// C++ dependencies
#include <filesystem>
#include <regex>

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