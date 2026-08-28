// C++ dependencies
#include <chrono>
#include <cstdlib>
#include <format>
#include <iostream>
#include <string>

// X17 dependencies
#include "Logger.h"

namespace X17
{
  void Logger::PrintTime(bool brackets)
  {
    auto now                 = std::chrono::system_clock::now();
    const char left_bracket  = brackets ? '[' : ' ';
    const char right_bracket = brackets ? ']' : ' ';
    std::string now_str      = std::format("{:%H:%M:%S}", std::chrono::zoned_time{ std::chrono::current_zone(), now });
    now_str                  = now_str.substr(0, 8);
    std::cout << left_bracket << now_str << right_bracket << " ";
  }

  void Logger::Log(LogLevel level, const std::string& message) const
  {
    for (int i = 0; i < m_indent; i++)
      std::cout << "  ";

    if (level != LogLevel::Plain)
    {
      PrintTime(true);
    }

    switch (level)
    {
      case LogLevel::Debug:
        std::cout << "\033[1;37mDEBUG:\033[0m ";
        break;
      case LogLevel::Info:
        std::cout << "\033[1;34mINFO:\033[0m ";
        break;
      case LogLevel::Essential:
        std::cout << "\033[1;32mESSENTIAL:\033[0m ";
        break;
      case LogLevel::Warning:
        std::cout << "\033[1;33mWARNING:\033[0m ";
        break;
      case LogLevel::Error:
        std::cout << "\033[1;31mERROR:\033[0m ";
        break;
      case LogLevel::Fatal:
        std::cout << "\033[1;35mFATAL:\033[0m ";
        break;
      case LogLevel::Plain:
        break;
    }

    std::cout << message << "\n";

    if (level >= LogLevel::Error)
      std::cout << std::flush;

    if (level == LogLevel::Fatal)
      std::exit(1);
  }
} // namespace X17