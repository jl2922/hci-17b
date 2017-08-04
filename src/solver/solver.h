#ifndef SOLVER_H_
#define SOLVER_H_

#include "../std.h"
#include "../wavefunction/wavefunction.h"

class Solver {
 protected:
  size_t n_up;
  size_t n_dn;
  double max_abs_H;
  Wavefunction wf;
  double energy_hf;
  double energy_var;

  virtual void solve() {}

  virtual double hamiltonian(const Det&, const Det&) const = 0;

  void variation();

  Det generate_hf_det();
};

#endif