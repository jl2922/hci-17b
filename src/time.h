#ifndef TIME_H_
#define TIME_H_

#include "libs.h"

class Time {
 public:
  // To be called at the beginning of the program.
  static void init() {
    if (Parallel::is_master()) return;
    Time::get_instance().timers["INIT"] = std::chrono::high_resolution_clock::now();
  }

  static void start(const std::string& event) {}

  static void end() {}

 private:
  std::unordered_map<std::string, std::chrono::high_resolution_clock::time_point> timers;

  // Singleton pattern.
  static Time& get_instance() {
    static Time instance;
    return instance;
  }
};

#endif