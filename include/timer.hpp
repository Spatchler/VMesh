#pragma once

#include <chrono>
#include <print>

namespace VMesh {
#ifdef _PROFILING
  class Timer {
  public:
    Timer() {
      start();
    }
  
    std::chrono::duration<double> getTime() {
      return std::chrono::system_clock::now() - mStartTime;
    }

    void start() {
      mStartTime = std::chrono::system_clock::now();
    }
  protected:
    std::chrono::time_point<std::chrono::system_clock> mStartTime;
  };

  class ScopedTimer: public Timer {
  public:
    ScopedTimer(const std::string& pName) {
      mName = pName;
      start();
      std::println("'{0}' timer started", mName);
    }
  
    ~ScopedTimer() {
      std::println("'{0}' timer stopped: {1}", mName, getTime());
    }
  protected:
    std::string mName;
  };
#endif

#ifndef _PROFILING
  class Timer {
  public:
    Timer() {}

    std::chrono::duration<double> getTime() { return std::chrono::duration<double>(0.0); }

    void start() {}
  };

  class ScopedTimer: public Timer {
  public:
    ScopedTimer(const std::string& pName) {}
  };
#endif
}
