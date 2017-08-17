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
  bool end_variation;
  std::unordered_map<OrbitalsPair, size_t, boost::hash<OrbitalsPair>> var_dets_id_lut;
  std::unordered_map<OrbitalsPair, double, boost::hash<OrbitalsPair>> new_dets_coef_lut;
  ExcitationStore ex;

  virtual void solve() {}

  virtual double hamiltonian(const Det&, const Det&) const = 0;

  void variation(const double, const double, const double);

  Det generate_hf_det();

  virtual std::list<Det> find_connected_dets(const Det&, const double eps) const = 0;

  double diagonalize(const double, const double);

  std::vector<double> apply_hamiltonian(const std::vector<double>&, const double, const double);
};

#endif