#pragma once

#include <chrono>
#include <iostream>
#include <stack>
#include <string>
#include <utility>

#define PROFILE_SCOPE(name) Timer timer##__LINE__(name)
#define PROFILE_FUNCTION() PROFILE_SCOPE(__PRETTY_FUNCTION__)

class Timer
{
private:
  std::string m_name;
  std::chrono::high_resolution_clock::time_point m_start;

  static std::stack<Timer*>& ActiveTimers_()
  {
    static std::stack<Timer*> timers;
    return timers;
  }

public:
  explicit Timer(std::string name, bool addToStack = true)
    : m_name(std::move(name))
  {
    m_start = std::chrono::high_resolution_clock::now();
    if (addToStack)
      ActiveTimers_().push(this);
  }

  ~Timer()
  {
    auto endTime  = std::chrono::high_resolution_clock::now();
    auto duration = endTime - m_start;
    std::cout << m_name << " (real time): " << std::chrono::duration<double>(duration).count() << " s" << '\n';

    if (ActiveTimers_().top() == this)
      ActiveTimers_().pop();
  }

  static double GetTime()
  {
    if (ActiveTimers_().empty())
      return 0.0;

    auto endTime  = std::chrono::high_resolution_clock::now();
    auto duration = endTime - ActiveTimers_().top()->m_start;
    return std::chrono::duration<double>(duration).count();
  }
};