#include "solver.h"

#include "../parallel.h"
#include "../std.h"

Det Solver::generate_hf_det() {
  Det det;
  for (size_t i = 0; i < n_up; i++) det.up.set_orb(i, true);
  for (size_t i = 0; i < n_dn; i++) det.dn.set_orb(i, true);
  return det;
}

void Solver::variation() {
  const double THRESHOLD = 1.0e-6;

  // Setup HF or existing wf as initial wf and evaluate energy.
  if (wf.size() == 0) {
    const Det& det_hf = generate_hf_det();
    wf.append_term(det_hf, 1.0);
    energy_hf = energy_var = hamiltonian(det_hf, det_hf);
    if (Parallel::is_master()) printf("HF energy: %#.15g Ha\n", energy_hf);
  }
}