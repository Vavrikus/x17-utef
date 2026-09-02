#pragma once

// C++ dependencies
#include <cmath>
#include <string>
#include <vector>

// X17 dependencies
#include "Logger.h"

namespace X17
{
  /// @brief Struct representing a parameter for a parallel job.
  struct JobParameter
  {
    std::string name;
    double min_bound;
    double max_bound;
    double step_size;
    int steps;
  };

  /// @brief Singleton class for managing the parameters in parallel jobs.
  class JobManager
  {
  public:
    static JobManager& Get()
    {
      static JobManager instance;
      return instance;
    }

    JobManager(const JobManager&)            = delete;
    JobManager& operator=(const JobManager&) = delete;

    [[nodiscard]] int GetMinIndex() const { return m_min_index; }
    [[nodiscard]] int GetMaxIndex() const { return m_max_index; }
    [[nodiscard]] int GetIterations() const { return m_iterations; }

    /// @brief Registers a parameter for a parallel job using its step size.
    void RegisterParameterStepSize(const std::string& name, double min_bound, double max_bound, double step_size);

    /// @brief Registers a parameter for a parallel job using its number of steps.
    void RegisterParameterSteps(const std::string& name, double min_bound, double max_bound, int steps);

    /// @brief Initializes the manager for a specific job.
    void Initialize(int max_id, int id, int iter);

    Logger& PrintParameters();

    Logger& PrintParameterValues(int item_index);

    /// @brief Returns the value of a job parameter for a specific item.
    /// @param name The name of the parameter.
    /// @param item_index The index of the processed item.
    /// @return The value of the parameter for the specified item.
    [[nodiscard]] double GetParameterValue(const std::string& name, int item_index) const;

  private:
    JobManager() = default;

  private:
    bool m_initialized = false;
    bool m_use_cost    = false;

    int m_id         = -1; ///< Current unique job ID.
    int m_max_id     = -1; ///< Number of parallel jobs.
    int m_iterations = -1; ///< Number of iterations for each job.

    int m_min_index   = -1; ///< Index of the first processed item in this job (e.g., a track).
    int m_max_index   = -1; ///< Index of the last processed item in this job (e.g., a track).
    int m_total_items = -1; ///< Total number of items to process in all jobs (e.g., number of tracks).

    std::vector<JobParameter> m_parameters;
    std::size_t m_max_name_length = 0;
  };
} // namespace X17