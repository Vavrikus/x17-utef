#pragma once

// C++ dependencies
#include <concepts>
#include <cstdint>
#include <sstream>
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

  template <typename... Args>
  std::string build_string(Args&&... args)
  {
    std::ostringstream oss;
    (oss << ... << std::forward<Args>(args));
    return oss.str();
  }

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

    /// @brief Logs a message with a variable number of arguments.
    template <typename... Args>
    requires (sizeof...(Args) > 1)
    Logger& Log(LogLevel level, Args&&... args)
    {
      return Log(level, build_string(std::forward<Args>(args)...));
    }

    Logger& Debug(const std::string& message) { return Log(LogLevel::Debug, message); }
    Logger& Info(const std::string& message) { return Log(LogLevel::Info, message); }
    Logger& Essential(const std::string& message) { return Log(LogLevel::Essential, message); }
    Logger& Warning(const std::string& message) { return Log(LogLevel::Warning, message); }
    Logger& Error(const std::string& message) { return Log(LogLevel::Error, message); }
    Logger& Fatal(const std::string& message) { return Log(LogLevel::Fatal, message); }
    Logger& Print(const std::string& message) { return Log(LogLevel::Plain, message); }

    template <typename... Args>
    requires (sizeof...(Args) > 1)
    Logger& Debug(Args&&... args)
    {
      return Log(LogLevel::Debug, std::forward<Args>(args)...);
    }
    template <typename... Args>
    requires (sizeof...(Args) > 1)
    Logger& Info(Args&&... args)
    {
      return Log(LogLevel::Info, std::forward<Args>(args)...);
    }
    template <typename... Args>
    requires (sizeof...(Args) > 1)
    Logger& Essential(Args&&... args)
    {
      return Log(LogLevel::Essential, std::forward<Args>(args)...);
    }
    template <typename... Args>
    requires (sizeof...(Args) > 1)
    Logger& Warning(Args&&... args)
    {
      return Log(LogLevel::Warning, std::forward<Args>(args)...);
    }
    template <typename... Args>
    requires (sizeof...(Args) > 1)
    Logger& Error(Args&&... args)
    {
      return Log(LogLevel::Error, std::forward<Args>(args)...);
    }
    template <typename... Args>
    requires (sizeof...(Args) > 1)
    Logger& Fatal(Args&&... args)
    {
      return Log(LogLevel::Fatal, std::forward<Args>(args)...);
    }
    template <typename... Args>
    requires (sizeof...(Args) > 1)
    Logger& Print(Args&&... args)
    {
      return Log(LogLevel::Plain, std::forward<Args>(args)...);
    }

    void PushIndent(int increment = 1) { m_indent += increment; }
    void PopIndent(int decrement = 1) { m_indent -= decrement; }
    void ResetIndent() { m_indent = 0; }

    static void PrintTime(bool brackets = false);

    /// @brief Prints a progress bar to the console.
    /// @param current The current value of the progress.
    /// @param total The total value of the progress.
    /// @param width The width of the progress bar.
    template <Number T>
    void ProgressBar(T current, T total,
                     int width = 50); // NOLINT(readability-redundant-declaration)

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
