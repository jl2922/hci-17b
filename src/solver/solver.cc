#include "solver.h"

#include <boost/functional/hash.hpp>
#include "../parallel.h"
#include "../std.h"
#include "../time.h"
#include "../wavefunction/wavefunction.h"
#include "davidson.h"

Det Solver::generate_hf_det() {
  Det det;
  for (size_t i = 0; i < n_up; i++) det.up.set_orb(i, true);
  for (size_t i = 0; i < n_dn; i++) det.dn.set_orb(i, true);
  return det;
}

void Solver::variation(
    const double eps_var, const double eps_var_ham_old, const double eps_var_ham_new) {
  const double THRESHOLD = 1.0e-6;

  // Setup HF or existing wf as initial wf and evaluate energy.
  if (wf.size() == 0) {
    const Det& det_hf = generate_hf_det();
    wf.append_term(det_hf, 1.0);
    energy_hf = energy_var = hamiltonian(det_hf, det_hf);
    if (Parallel::is_master()) printf("HF energy: %#.15g Ha\n", energy_hf);
  }

  double energy_var_new = 0.0;  // Ensures the first iteration will run.
  std::unordered_set<OrbitalsPair, boost::hash<OrbitalsPair>> var_dets_set;
  for (const auto& term : wf.get_terms()) var_dets_set.insert(term.det.encode());
  int iteration = 0;  // For print.
  while (fabs(energy_var - energy_var_new) > THRESHOLD) {
    Time::start("Variation Iteration: " + std::to_string(iteration));

    // Find connected determinants.
    // Mapping from new det to spawning det coef.
    std::unordered_map<OrbitalsPair, double, boost::hash<OrbitalsPair>> new_dets_set;
    for (const auto& term : wf.get_terms()) {
      const double abs_coef = fabs(term.coef);
      const auto& connected_dets = find_connected_dets(term.det, eps_var / abs_coef);
      for (const auto& new_det : connected_dets) {
        const auto& new_det_code = new_det.encode();
        if (var_dets_set.count(new_det_code) == 0 && new_dets_set.count(new_det_code) == 0) {
          new_dets_set[new_det_code] = abs_coef;
        }
      }
    }

    if (Parallel::get_id() == 0) {
      printf(
          "Number of new / total dets: %'llu / %'llu\n",
          static_cast<unsigned long long>(new_dets_set.size()),
          static_cast<unsigned long long>(new_dets_set.size() + var_dets_set.size()));
    }
    Time::checkpoint("found new dets");

    energy_var = energy_var_new;

    // for (const auto& new_det_info : new_dets_set) {
    //   const double abs_coef = fabs(new_det_info.second);
    //   Det new_det;
    //   new_det.decode(new_det_info.first);
    //   const auto& connected_dets = find_connected_dets(new_det, eps_var_ham_new / abs_coef);
    //   for (const auto& conn_det : connected_dets) {
    //     const auto& code = conn_det.encode();
    //     if (var_dets_set.count(code) == 0 && new_dets_set.count(code) == 0) continue;
    //     // update_helper_list(new_det, conn_det);
    //   }
    // }
    // Time::checkpoint("generated excitation store");

    for (const auto& new_det_info : new_dets_set) {
      const auto& code = new_det_info.first;
      var_dets_set.insert(code);
      Det det;
      det.decode(code);
      wf.append_term(det, 0.0);
    }

    energy_var_new =
        diagonalize(new_dets_set.size() > 0 ? 3 : 10, eps_var_ham_old, eps_var_ham_new);
    if (Parallel::get_id() == 0) printf("Variation energy: %#.15g Ha\n", energy_var_new);

    iteration++;

    Time::end();
  }

  energy_var = energy_var_new;
  if (Parallel::get_id() == 0) printf("Final variation energy: %#.15g Ha\n", energy_var);
}

double Solver::diagonalize(
    const std::size_t max_iterations, const double eps_var_ham_old, const double eps_var_ham_new) {
  std::vector<double> diagonal;
  std::vector<double> initial_vector;
  diagonal.reserve(wf.size());
  initial_vector.reserve(wf.size());
  for (const auto& term : wf.get_terms()) {
    const auto& det = term.det;
    diagonal.push_back(hamiltonian(det, det));
    initial_vector.push_back(term.coef);
  }

  Time::start("Diagonalization");
  std::function<std::vector<double>(std::vector<double>)> apply_hamiltonian_func =
      std::bind(&Solver::apply_hamiltonian, this, std::placeholders::_1);

  Davidson davidson(diagonal, apply_hamiltonian_func, wf.size());
  if (Parallel::get_id() == 0) davidson.set_verbose(true);
  davidson.diagonalize(initial_vector, max_iterations);
  Time::end();

  const double energy_var = davidson.get_lowest_eigenvalue();
  const auto& coefs_new = davidson.get_lowest_eigenvector();
  wf.set_coefs(coefs_new);
  wf.sort_by_coefs();

  return energy_var;
}

std::vector<double> Solver::apply_hamiltonian(const std::vector<double>& vec) {
  std::size_t n = vec.size();
  assert(n == wf.size());
  std::vector<double> res(n, 0.0);
  const auto& dets = wf.get_dets();

  for (size_t i = Parallel::get_id(); i < n; i += Parallel::get_n()) {
    const Det& det_i = dets[i];
    const double H_ii = hamiltonian(det_i, det_i);
  }
  Time::checkpoint("hamiltonian applied");

  Parallel::reduce_to_sum(res);
  return res;
}