// C++ dependencies
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>

#ifdef __linux__
#include <unistd.h>

#include <linux/limits.h>
#include <sys/types.h>
#elif _WIN32
#include <windows.h>
#elif __APPLE__
#include <limits.h>

#include <mach-o/dyld.h>
#endif

// ROOT dependencies
#include "RtypesCore.h"
#include "TFile.h"
#include "TParameter.h"
#include "TRandom3.h"
#include "TTree.h"

// X17 dependencies
#include "AppManager.h"
#include "Field.h"
#include "Logger.h"

namespace X17
{
  using std::string;

  AppManager::AppManager(string app_name, int filepath_depth, string app_description)
    : m_app_name(std::move(app_name)), m_filepath_depth(filepath_depth), m_app_description(std::move(app_description))
  {
    auto print_line = [](const std::string& text)
    {
      constexpr int WIDTH = 72;
      std::cout << std::left << std::setw(WIDTH) << text << "🭵\n";
    };

#ifndef NDEBUG
    string build_type = "Debug";
#else
    string build_type = "RelWithDebInfo";
#endif

    std::cout << "\n –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––\n";

    print_line(" 🭰 Welcome to the X17 OFTPC simulation and reconstruction framework.");
    print_line(" 🭰 Built on " + string(__DATE__) + " at " + string(__TIME__) + " as " + build_type);

    std::cout << " 🭰 Current time: ";
    Logger::PrintTime(false);
    std::cout << std::setw(41) << "" << "🭵\n";

    print_line(" 🭰 " + m_app_name + " - " + m_app_description);

    std::cout << " –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––\n\n";

    ToRootDir_();
  }

  AppManager::~AppManager()
  {
    m_output_tree->Write();
    m_output_file->Close();
    for (auto& file : m_open_files)
      file->Close();
  }
  std::unique_ptr<MagField> AppManager::LoadMagField()
  {
    return std::unique_ptr<MagField>(X17::LoadField("data/elmag/VecB2.txt", { -20, -30, -30 }, { 20, 30, 30 }, 0.5));
  }

  void AppManager::ToRootDir_() const
  {
    namespace fs = std::filesystem;

    fs::path bin_dir;
    try
    {
#ifdef __linux__
      char result[PATH_MAX];
      ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
      if (count == -1)
        throw std::runtime_error("Failed to read /proc/self/exe");
      bin_dir = fs::path(string(result, count));

#elif _WIN32
      char result[MAX_PATH];
      GetModuleFileNameA(NULL, result, MAX_PATH);
      bin_dir = fs::path(result);

#elif __APPLE__
      char path[PATH_MAX];
      uint32_t size = sizeof(path);
      if (_NSGetExecutablePath(path, &size) == 0)
      {
        // macOS does not always return a clean canonical path here,
        // so we canonicalize it just to be safe.
        return fs::canonical(fs::path(path));
      }
      throw std::runtime_error("Failed to get Mac executable path");

#else
      throw std::runtime_error("OS not supported");
#endif

      fs::path target_dir = bin_dir.parent_path();

      for (int i = 0; i < m_filepath_depth; i++)
        target_dir = target_dir.parent_path();

      fs::current_path(fs::weakly_canonical(target_dir));

      Logger::Get().Info("Setting working directory to " + target_dir.string());
    }
    catch (const std::exception& e)
    {
      std::string msg = "Failed to set working directory: " + string(e.what());
      Logger::Get().Error(msg);
    }
  }

  std::unique_ptr<DriftMap> AppManager::LoadMap(const string& version)
  {
    Logger::Get().Info("Loading map: " + version + ".");
    TFile* map_input = TFile::Open(("data/ion_map/sample_" + version + "/map.root").c_str(), "READ");
    m_open_files.emplace_back(map_input);
    return std::unique_ptr<DriftMap>(reinterpret_cast<X17::DriftMap*>(map_input->Get("map")));
  }

  TTree* AppManager::LoadTreeFromFile(const string& filename, const string& tree_name)
  {
    Logger::Get().Essential("Loading tree: " + tree_name + " from " + filename + ".");
    TFile* tree_input = TFile::Open(filename.c_str(), "READ");
    m_open_files.emplace_back(tree_input);
    return static_cast<TTree*>(tree_input->Get(tree_name.c_str()));
  }

  TTree* AppManager::CreateOutputTree(const string& filename, const string& tree_name, const string& tree_title)
  {
    Logger::Get().Info("Creating output tree: " + tree_name + " in " + filename + ".");
    m_output_file = std::make_unique<TFile>(filename.c_str(), "RECREATE");
    m_output_tree = new TTree(tree_name.c_str(), tree_title.c_str());
    m_output_tree->SetDirectory(m_output_file.get());
    return m_output_tree;
  }

  std::unique_ptr<TRandom3> AppManager::CreateRNG(UInt_t seed)
  {
    if (seed == 0)
      seed = std::random_device{}();
    X17::Logger::Get().Info("Setting TRandom3 seed: " + std::to_string(seed));

    if (!m_output_file)
      Logger::Get().Fatal("No output file set, cannot save seed.");

    TParameter<UInt_t> pSeed("seed", seed);
    pSeed.Write();

    return std::make_unique<TRandom3>(seed);
  }
} // namespace X17