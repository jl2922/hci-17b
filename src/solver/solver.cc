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
    std::unordered_map<OrbitalsPair, double, boost::hash<OrbitalsPair>> new_dets_set;
    for (const auto& term : wf.get_terms()) {
      const double abs_coef = fabs(term.coef);
      const auto& connected_dets = find_connected_dets(term.det, eps_var / abs_coef);
      for (const auto& new_det : connected_dets) {
        const auto& code = new_det.encode();
        if (var_dets_set.count(code) == 0 && new_dets_set.count(code) == 0) {
          new_dets_set[code] = abs_coef;
        }
      }
    }

    if (Parallel::get_id() == 0) {
      printf(
          "Number of new / total dets: %'llu / %'llu\n",
          static_cast<unsigned long long>(new_dets_set.size()),
          static_cast<unsigned long long>(new_dets_set.size() + var_dets_set.size()));
    }

    energy_var = energy_var_new;

    // update helper lists.
    Time::checkpoint("found connections");
    for (const auto& term : wf.get_terms()) {
      const double abs_coef = fabs(term.coef);
      const auto& connected_dets = find_connected_dets(term.det, eps_var_ham_old / abs_coef);
      for (const auto& conn_det : connected_dets) {
        const auto& code = conn_det.encode();
        if (var_dets_set.count(code) == 0 && new_dets_set.count(code) == 0) continue;
        update_helper_list(term.det, conn_det);
      }
    }
    for (const auto& new_det_info : new_dets_set) {
      const double abs_coef = fabs(new_det_info.second);
      Det new_det;
      new_det.decode(new_det_info.first);
      const auto& connected_dets = find_connected_dets(new_det, eps_var_ham_new / abs_coef);
      for (const auto& conn_det : connected_dets) {
        const auto& code = conn_det.encode();
        if (var_dets_set.count(code) == 0 && new_dets_set.count(code) == 0) continue;
        update_helper_list(new_det, conn_det);
      }
    }
    Time::checkpoint("generated excitation store");

    for (const auto& new_det_info : new_dets_set) {
      const auto& code = new_det_info.first;
      var_dets_set.insert(code);
      Det det;
      det.decode(code);
      wf.append_term(det, 0.0);
    }

    energy_var_new = diagonalize(new_dets_set.size() > 0 ? 5 : 10);
    if (Parallel::get_id() == 0) printf("Variation energy: %#.15g Ha\n", energy_var_new);

    iteration++;

    Time::end();
  }

  energy_var = energy_var_new;
  if (Parallel::get_id() == 0) printf("Final variation energy: %#.15g Ha\n", energy_var);
}

void Solver::update_helper_list(const Det& det_1, const Det& det_2) {
  if (det_1.up == det_2.up && det_1.dn != det_2.dn) {
    ex.add(det_1.dn.encode(), det_2.dn.encode(), true);
  } else if (det_1.up != det_2.up && det_1.dn == det_2.dn) {
    ex.add(det_1.up.encode(), det_2.up.encode(), true);
  } else if (det_1.up != det_2.up && det_1.dn != det_2.dn) {
    ex.add(det_1.up.encode(), det_2.up.encode(), false);
    ex.add(det_1.dn.encode(), det_2.dn.encode(), true);
  }
}

double Solver::diagonalize(std::size_t max_iterations) {
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
  double energy_var = davidson.get_lowest_eigenvalue();
  const auto& coefs_new = davidson.get_lowest_eigenvector();
  wf.set_coefs(coefs_new);
  wf.sort_by_coefs();

  return energy_var;
}

std::vector<double> Solver::apply_hamiltonian(const std::vector<double>& vec) {
  std::size_t n = vec.size();
  assert(n == wf.size());
  std::vector<double> res(n, 0.0);
  //   std::vector<long double> res_precise(n, 0.0);
  const auto& dets = wf.get_dets();
  std::unordered_map<
      std::pair<uint32_t, uint32_t>,
      size_t,
      boost::hash<std::pair<uint32_t, uint32_t>>>
      dets_code_set;
  for (size_t i = 0; i < dets.size(); i++) {
    const Det& det = dets[i];
    const uint32_t code_up = ex.lookup_id(det.up.encode());
    const uint32_t code_dn = ex.lookup_id(det.dn.encode());
    dets_code_set[std::make_pair(code_up, code_dn)] = i;
  }
  //   unsigned long long n_connections = 0;
  //   static unsigned long long n_connections_prev = 0;
  //   unsigned long long same_spin_count = 0, opposite_spin_count = 0;

  for (std::size_t i = 0; i < n; i++) {
    if (i % Parallel::get_n() != static_cast<std::size_t>(Parallel::get_id())) continue;
    const Det& det_i = dets[i];
    const uint32_t code_up = ex.lookup_id(det_i.up.encode());
    const uint32_t code_dn = ex.lookup_id(det_i.dn.encode());
    const double H_ii = hamiltonian(det_i, det_i);
    res[i] += H_ii * vec[i];
    for (const uint32_t code_up_conn : ex.find(code_up, true)) {
      if (dets_code_set.count(std::make_pair(code_up_conn, code_dn))) {
        const size_t j = dets_code_set[std::make_pair(code_up_conn, code_dn)];
        if (j < i) continue;
        const Det& det_j = dets[j];
        const double H_ij = hamiltonian(det_i, det_j);
        res[i] += H_ij * vec[j];
        res[j] += H_ij * vec[i];
      }
    }
    for (const uint32_t code_dn_conn : ex.find(code_dn, true)) {
      if (dets_code_set.count(std::make_pair(code_up, code_dn_conn))) {
        const size_t j = dets_code_set[std::make_pair(code_up, code_dn_conn)];
        if (j < i) continue;
        const Det& det_j = dets[j];
        const double H_ij = hamiltonian(det_i, det_j);
        res[i] += H_ij * vec[j];
        res[j] += H_ij * vec[i];
      }
    }
    for (const uint32_t code_up_conn : ex.find(code_up, false)) {
      for (const uint32_t code_dn_conn : ex.find(code_dn, false)) {
        if (dets_code_set.count(std::make_pair(code_up_conn, code_dn_conn))) {
          const size_t j = dets_code_set[std::make_pair(code_up_conn, code_dn_conn)];
          if (j < i) continue;
          const Det& det_j = dets[j];
          const double H_ij = hamiltonian(det_i, det_j);
          res[i] += H_ij * vec[j];
          res[j] += H_ij * vec[i];
        }
      }
    }
  }
  Time::checkpoint("hamiltonian applied");

  Parallel::reduce_to_sum(res);
  //   Parallel::reduce_to_sum(n_connections);
  //   Parallel::reduce_to_sum(same_spin_count);
  //   Parallel::reduce_to_sum(opposite_spin_count);
  //   if (Parallel::get_id() == 0 && n_connections != n_connections_prev) {
  //     printf("Same spin: %'llu, opposite spin: %'llu\n", same_spin_count, opposite_spin_count);
  //     printf("Number of connections: %'llu\n", n_connections);
  //     n_connections_prev = n_connections;
  //   }
  //   for (std::size_t i = 0; i < n; i++) res[i] = static_cast<double>(res_precise[i]);
  return res;
}