#ifndef PARALLEL_H_
#define PARALLEL_H_

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
#include "libs.h"

#ifndef SERIAL
class Parallel {
 private:
  int id;
  int n;
  boost::mpi::environment* env;  // For MPI 1.1.
  boost::mpi::communicator world;

  Parallel() {
    id = this->world.rank();
    n = this->world.size();
  }

  // Singleton pattern boilerplate.
  static Parallel& get_instance() {
    static Parallel instance;
    return instance;
  }

 public:
  static void init(boost::mpi::environment& env) { Parallel::get_instance().env = &env; }

  static int get_id() { return Parallel::get_instance().id; }

  static bool is_master() { return Parallel::get_instance().id == 0; }

  static std::string get_host() { return Parallel::get_instance().env->processor_name(); }

  static void barrier() {
    fflush(stdout);
    Parallel::get_instance().world.barrier();
  }
};
#else
// Non-MPI stub for debugging and profiling.
class Parallel {
 public:
  static bool is_master() { return true; }

  static int get_id() { return 0; }

  static std::string get_host() { return "localhost"; }

  static void barrier() { return; }
};
#endif

#endif