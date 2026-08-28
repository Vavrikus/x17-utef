#pragma once

// C++ dependencies
#include <memory>
#include <string>
#include <vector>

// ROOT dependencies
#include "TFile.h"
#include "TTree.h"

// X17 dependencies
#include "Field.h"

namespace X17
{
  /// @brief Class for managing the application and I/O.
  class AppManager
  {
  public:
    /// @brief Constructor. Prints a welcome message and sets the working directory.
    /// @param app_name The name of the application.
    /// @param filepath_depth The depth of the binary filepath from the project root.
    /// @param app_description The description of the application.
    explicit AppManager(std::string app_name, int filepath_depth = 2, std::string app_description = "");

    /// @brief Destructor. Writes the output tree and closes the ROOT files.
    ~AppManager();

    static std::unique_ptr<MagField> LoadMagField();

    std::unique_ptr<DriftMap> LoadMap(const std::string& version);

    TTree* LoadTreeFromFile(const std::string& filename, const std::string& tree_name);

    TTree* CreateOutputTree(const std::string& filename, const std::string& tree_name,
                            const std::string& tree_title = "");

    /// @brief Creates a ROOT random number generator.
    /// @param seed The seed for the random number generator. If not provided or set to 0, a random seed will be used.
    /// @return A unique pointer to the random number generator.
    std::unique_ptr<TRandom3> CreateRNG(UInt_t seed = 0);

  private:
    /// @brief Sets the working directory to the root of the project.
    void ToRootDir_() const;

  private:
    std::string m_app_name;
    std::string m_app_description;
    int m_filepath_depth; ///< The depth of the binary filepath from the project root.

    std::vector<std::unique_ptr<TFile>> m_open_files;
    std::unique_ptr<TFile> m_output_file;
    TTree* m_output_tree;
  };
} // namespace X17