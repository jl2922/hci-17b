#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
#include "config.h"
#include "heg_solver/heg_solver.h"
#include "libs.h"
#include "parallel.h"
#include "time.h"

int main(int argc, char** argv) {
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

  Time::init();

  Config::load("config.json");
  if (Parallel::is_master()) {
    Config::print();
  }

  const std::string& type = Config::get<std::string>("type");
  if (type == "heg") {
    HEGSolver::run();
  } else {
    throw std::invalid_argument("System type not supported");
  }

  return 0;
}