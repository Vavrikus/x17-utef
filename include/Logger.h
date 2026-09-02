#pragma once

// C++ dependencies
#include <concepts>
#include <cstdint>
#include <string>

namespace X17
{
  template <typename T>
  concept Number = std::integral<T> || std::floating_point<T>;

  /// @brief Enumeration of log levels.
  enum class LogLevel : std::uint8_t
  {
    Debug,
    Info,
    Essential,
    Warning,
    Error,
    Fatal,
    Plain,
  };

  /// @brief Singleton class for logging messages.
  class Logger
  {
  public:
    static Logger& Get()
    {
      static Logger instance;
      return instance;
    }

    Logger(const Logger&)            = delete;
    Logger& operator=(const Logger&) = delete;

    void SetLogLevel(LogLevel level) { m_log_level = level; }

    Logger& Log(LogLevel level, const std::string& message);
    Logger& Debug(const std::string& message) { return Log(LogLevel::Debug, message); }
    Logger& Info(const std::string& message) { return Log(LogLevel::Info, message); }
    Logger& Essential(const std::string& message) { return Log(LogLevel::Essential, message); }
    Logger& Warning(const std::string& message) { return Log(LogLevel::Warning, message); }
    Logger& Error(const std::string& message) { return Log(LogLevel::Error, message); }
    Logger& Fatal(const std::string& message) { return Log(LogLevel::Fatal, message); }
    Logger& Print(const std::string& message) { return Log(LogLevel::Plain, message); }

    void PushIndent(int increment = 1) { m_indent += increment; }
    void PopIndent(int decrement = 1) { m_indent -= decrement; }
    void ResetIndent() { m_indent = 0; }

    static void PrintTime(bool brackets = false);

    /// @brief Prints a progress bar to the console.
    /// @param current The current value of the progress.
    /// @param total The total value of the progress.
    /// @param width The width of the progress bar.
    template <Number T>
    void ProgressBar(T current, T total, int width = 50); // NOLINT(readability-redundant-declaration)

  private:
    Logger() = default;

  private:
    LogLevel m_log_level = LogLevel::Info;
    int m_indent         = 0;

    LogLevel m_last_level = LogLevel::Debug; ///< The last non-plain log level.
  };
} // namespace X17

// Templated function definitions.
#ifndef LOGGER_INL
#include "Logger.inl"
#endif
