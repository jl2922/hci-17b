#include "heg_solver.h"

#include "../config.h"
#include "../parallel.h"
#include "../time.h"

void HEGSolver::solve() {
  printf("Proc %d running on %s\n", Parallel::get_id(), Parallel::get_host().c_str());
  n_up = Config::get<std::size_t>("n_up");
  n_dn = Config::get<std::size_t>("n_dn");
  rcut_vars = Config::get_array<double>("rcut_vars");
  eps_vars = Config::get_array<double>("eps_vars");
  rcut_pts = Config::get_array<double>("rcut_pts");
  eps_pts = Config::get_array<double>("eps_pts");

  // Check configuration validity.
  assert(rcut_vars.size() == rcut_pts.size());
  assert(eps_vars.size() == eps_pts.size());
  for (std::size_t i = 1; i < rcut_vars.size(); i++) {
    assert(rcut_vars[i - 1] <= rcut_vars[i]);
  }
  for (std::size_t i = 1; i < eps_vars.size(); i++) {
    assert(eps_vars[i - 1] >= eps_vars[i]);
  }

  Time::start("variation");
  for (const double rcut_var : rcut_vars) {
    std::string rcut_var_event = str(boost::format("rcut_var: %#.4g") % rcut_var);
    Time::start(rcut_var_event);
    // this->rcut_var = rcut_var;
    // setup();
    for (const double eps_var : eps_vars) {
      std::string eps_var_event = str(boost::format("eps_var: %#.4g") % eps_var);
      Time::start(eps_var_event);
      // this->eps_var = eps_var;
      // if (!load_variation_result()) {
      //   variation();
      //   save_variation_result();
      // }
      Time::end();
    }
    Time::end();
  }
  Time::end();

  // Time::start("perturbation stage");
  // n_orbs_pts.clear();
  // for (const double rcut_pt : rcut_pts) {
  //   n_orbs_pts.push_back(KPointsUtil::get_n_k_points(rcut_pt) * 2);
  // }
  // // Start from the largest PT so that it fails earlier upon insufficient memory.
  // for (const double rcut_var : rcut_vars | boost::adaptors::reversed) {
  //   std::string rcut_var_event = str(boost::format("perturbation with rcut_var: %#.4g") %
  //   rcut_var);
  //   Time::start(rcut_var_event);
  //   this->rcut_var = rcut_var;
  //   for (const double eps_var : eps_vars | boost::adaptors::reversed) {
  //     std::string eps_var_event = str(boost::format("perturbation with eps_var: %#.4g") %
  //     eps_var);
  //     Time::start(eps_var_event);
  //     this->eps_var = eps_var;
  //     assert(load_variation_result());
  //     perturbation();
  //     Time::end(eps_var_event);
  //   }
  //   Time::end(rcut_var_event);
  // }
  // Time::end("perturbation stage");

  // Time::start("extrapolation");
  // extrapolate();
  // Time::end("extrapolation");
}