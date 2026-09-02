// C++ dependencies
#include <algorithm>
#include <utility>

// X17 dependencies
#include "JobManager.h"

namespace X17
{
  void JobManager::RegisterParameterStepSize(const std::string& name, double min_bound, double max_bound,
                                             double step_size)
  {
    if (m_initialized)
      Logger::Get().Fatal("Registering job parameter with step size after manager initialization.");
    m_parameters.emplace_back(name, min_bound, max_bound, step_size, std::ceil((max_bound - min_bound) / step_size));
    m_max_name_length = std::max<std::size_t>(name.length(), m_max_name_length);
  }

  void JobManager::RegisterParameterSteps(const std::string& name, double min_bound, double max_bound, int steps)
  {
    if (m_initialized)
      Logger::Get().Fatal("Registering job parameter with number of steps after manager initialization.");
    m_parameters.emplace_back(name, min_bound, max_bound, (max_bound - min_bound) / (steps - 1), steps);
    m_max_name_length = std::max<std::size_t>(name.length(), m_max_name_length);
  }

  void JobManager::Initialize(int max_id, int id, int iter)
  {
    m_max_id     = max_id;
    m_id         = id;
    m_iterations = iter;

    m_total_items = 1;
    for (auto& parameter : m_parameters)
      m_total_items *= parameter.steps;

    if (!m_use_cost)
    {
      // Integer division to avoid rounding errors.
      m_min_index = (m_total_items * (m_id - 1)) / m_max_id + 1;
      m_max_index = (m_total_items * m_id) / m_max_id;
    }
    else
    {
      Logger::Get().Fatal("JobManager::Initialize: Cost function not implemented.");
    }

    Logger::Get()
      .Essential(std::string("Initializing job ") + std::to_string(m_id) + " of " + std::to_string(m_max_id) + ".")
      .PushIndent();
    Logger::Get().Print("Processing items " + std::to_string(m_min_index) + " to " + std::to_string(m_max_index)
                        + " out of " + std::to_string(m_total_items) + ".");
    Logger::Get().Print("Parameters:").PushIndent();
    PrintParameters().PopIndent(2);

    m_initialized = true;
  }

  Logger& JobManager::PrintParameters()
  {
    Logger& logger = Logger::Get();
    for (const auto& parameter : m_parameters)
    {
      const std::string padding(m_max_name_length - parameter.name.size() + 1, ' ');
      logger.Print(parameter.name + ":" + padding + std::to_string(parameter.steps) + " steps from "
                   + std::to_string(parameter.min_bound) + " to " + std::to_string(parameter.max_bound) + " (step size "
                   + std::to_string(parameter.step_size) + ").");
    }
    return logger;
  }

  Logger& JobManager::PrintParameterValues(int item_index)
  {
    Logger& logger = Logger::Get();
    for (const auto& parameter : m_parameters)
    {
      const std::string padding(m_max_name_length - parameter.name.size() + 1, ' ');
      logger.Print(parameter.name + ":" + padding + std::to_string(GetParameterValue(parameter.name, item_index)));
    }
    return logger;
  }

  [[nodiscard]] double JobManager::GetParameterValue(const std::string& name, int item_index) const
  {
    if (!m_initialized)
      Logger::Get().Fatal("Accessing job parameter before manager initialization.");

    int i_param     = -1;
    int steps_below = 1;

    for (int i = static_cast<int>(m_parameters.size()) - 1; i >= 0; i--)
    {
      if (m_parameters[i].name == name)
      {
        i_param = i;
        break;
      }

      steps_below *= m_parameters[i].steps;
    }

    if (std::cmp_equal(i_param, -1))
      Logger::Get().Fatal("Unknown job parameter: " + name);

    int steps = m_parameters[i_param].steps;
    int index = (item_index - 1) / steps_below % steps;

    return m_parameters[i_param].min_bound + index * m_parameters[i_param].step_size;
  }
} // namespace X17