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
  const double THRESHOLD = 1.0e-5;

  // Setup HF or existing wf as initial wf and evaluate energy.
  if (wf.size() == 0) {
    const Det& det_hf = generate_hf_det();
    wf.append_term(det_hf, 1.0);
    energy_hf = energy_var = hamiltonian(det_hf, det_hf);
    if (Parallel::is_master()) printf("HF energy: %#.15g Ha\n", energy_hf);
  }

  double energy_var_new = 0.0;  // Ensures the first iteration will run.

  int iteration = 0;  // For print.
  end_variation = false;
  while (fabs(energy_var - energy_var_new) > THRESHOLD && !end_variation) {
    Time::start("Variation Iteration: " + std::to_string(iteration));

    var_dets_id_lut.clear();
    size_t det_id = 0;
    for (const auto& term : wf.get_terms()) {
      var_dets_id_lut.insert({term.det.encode(), det_id++});
    }

    // Find connected determinants.
    // Mapping from new det to spawning det coef.
    new_dets_coef_lut.clear();
    for (const auto& term : wf.get_terms()) {
      const double abs_coef = fabs(term.coef);
      const auto& connected_dets = find_connected_dets(term.det, eps_var / abs_coef);
      for (const auto& new_det : connected_dets) {
        const auto& new_det_code = new_det.encode();
        if (var_dets_id_lut.count(new_det_code) == 0 &&
            new_dets_coef_lut.count(new_det_code) == 0) {
          new_dets_coef_lut[new_det_code] = abs_coef;
        }
      }
    }

    if (Parallel::get_id() == 0) {
      printf(
          "Number of new / total dets: %'llu / %'llu\n",
          static_cast<unsigned long long>(new_dets_coef_lut.size()),
          static_cast<unsigned long long>(new_dets_coef_lut.size() + var_dets_id_lut.size()));
    }
    Time::checkpoint("found new dets");

    energy_var = energy_var_new;

    for (const auto& new_det_info : new_dets_coef_lut) {
      const auto& code = new_det_info.first;
      var_dets_id_lut.insert({code, det_id++});
      Det det;
      det.decode(code);
      wf.append_term(det, 0.0);
    }

    energy_var_new = diagonalize(eps_var_ham_old, eps_var_ham_new);
    if (Parallel::get_id() == 0) printf("Variation energy: %#.15g Ha\n", energy_var_new);

    iteration++;

    Time::end();
  }

  energy_var = energy_var_new;
  if (Parallel::get_id() == 0) printf("Final variation energy: %#.15g Ha\n", energy_var);
}

double Solver::diagonalize(const double eps_var_ham_old, const double eps_var_ham_new) {
  const size_t max_iterations = new_dets_coef_lut.size() > 0 ? 5 : 10;
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
  std::function<std::vector<double>(std::vector<double>)> apply_hamiltonian_func = std::bind(
      &Solver::apply_hamiltonian, this, std::placeholders::_1, eps_var_ham_old, eps_var_ham_new);

  Davidson davidson(diagonal, apply_hamiltonian_func, wf.size());
  if (Parallel::get_id() == 0) davidson.set_verbose(true);
  const int n_iter = davidson.diagonalize(initial_vector, max_iterations);
  if (n_iter == 10) end_variation = true;
  Time::end();

  const double energy_var = davidson.get_lowest_eigenvalue();
  const auto& coefs_new = davidson.get_lowest_eigenvector();
  wf.set_coefs(coefs_new);
  wf.sort_by_coefs();

  return energy_var;
}

std::vector<double> Solver::apply_hamiltonian(
    const std::vector<double>& vec, const double eps_var_ham_old, const double eps_var_ham_new) {
  size_t n = vec.size();
  assert(n == wf.size());
  std::vector<double> res(n, 0.0);
  const auto& dets = wf.get_dets();
  const auto& coefs = wf.get_coefs();
  size_t n_old_dets = n - new_dets_coef_lut.size();

  for (size_t i = Parallel::get_id(); i < n; i += Parallel::get_n()) {
    const Det& det_i = dets[i];
    const auto& det_i_code = det_i.encode();
    const bool is_old_det = i < n_old_dets;
    const double eps_var_ham = is_old_det ? eps_var_ham_old : eps_var_ham_new;
    const double abs_coef = is_old_det ? coefs[i] : new_dets_coef_lut[det_i_code];
    const auto& connected_dets = find_connected_dets(det_i, eps_var_ham / abs_coef);
    for (const auto& det_j : connected_dets) {
      const auto& det_j_code = det_j.encode();
      if (var_dets_id_lut.count(det_j_code)) {
        const size_t j = var_dets_id_lut.at(det_j_code);
        if (j < i) continue;
        const double H_ij = hamiltonian(det_i, det_j);
        res[i] += H_ij * vec[j];
        if (j != i) {
          res[j] += H_ij * vec[i];
        }
      }
    }
  }
  Time::checkpoint("hamiltonian applied");
#ifdef __INTEL_COMPILER
  for (std::size_t i = 0; i < n; i++) Parallel::reduce_to_sum(res[i]);
#else
  Parallel::reduce_to_sum(res);
#endif
  return res;
}