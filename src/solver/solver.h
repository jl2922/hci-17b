#ifndef SOLVER_H_
#define SOLVER_H_

#include <boost/functional/hash.hpp>
#include "../std.h"
#include "../wavefunction/wavefunction.h"
#include "excitation_store.h"

class Solver {
 protected:
  size_t n_up;
  size_t n_dn;
  double max_abs_H;
  Wavefunction wf;
  double energy_hf;
  double energy_var;
  double energy_pt;
  ExcitationStore ex;

  virtual void solve() {}

  virtual double hamiltonian(const Det&, const Det&) const = 0;

  void variation(const double, const double, const double);

  Det generate_hf_det();

  virtual std::list<Det> find_connected_dets(const Det&, const double eps) const = 0;

  double diagonalize(const std::size_t, const double, const double);

  std::vector<double> apply_hamiltonian(const std::vector<double>&);
};

#endif