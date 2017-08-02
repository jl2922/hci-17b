#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
#include "config.h"
#include "libs.h"
#include "parallel.h"

int main(int argc, char **argv) {
#ifndef SERIAL
  boost::mpi::environment env(argc, argv);  // For MPI 1.1.
  Parallel::init(env);
#endif

  std::setlocale(LC_NUMERIC, "");

  if (Parallel::is_master()) {
    printf("Heat-Bath Configuration Interaction\n");
    const time_t start_time = time(0);
    printf("%s\n", asctime(localtime(&start_time)));
  }

  Config::load("config.json");
  if (Parallel::is_master()) {
    Config::print();
  }

  return 0;
}