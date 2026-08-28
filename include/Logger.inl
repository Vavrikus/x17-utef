#define LOGGER_INL

// C++ dependencies
#include <iomanip>
#include <iostream>

// X17 dependencies
#include "Logger.h" // NOLINT(misc-header-include-cycle)
#include "Utilities.h"

namespace X17
{
  template <Number T>
  void Logger::ProgressBar(T current, T total, int width)
  {
    current++;

    if (current % (total / 200) != 0 && current != total)
      return;

    float progress = static_cast<float>(current) / static_cast<float>(total);
    int filled     = ifloor(progress * static_cast<float>(width));

    for (int i = 0; i < m_indent; i++)
      std::cout << "  ";

    std::cout << "\r["; // Carriage return to overwrite the line
    for (int i = 0; i < width; ++i)
      std::cout << (i < filled ? '=' : ' ');

    std::cout << "] " << std::setw(3) << ifloor(progress * 100) << "%";
    std::cout.flush();

    if (current == total)
      std::cout << '\n';
  }
} // namespace X17