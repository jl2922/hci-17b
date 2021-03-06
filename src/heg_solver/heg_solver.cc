#include "heg_solver.h"

#include <boost/functional/hash.hpp>
#include "../array_math.h"
#include "../config.h"
#include "../parallel.h"
#include "../time.h"
#include "k_points_util.h"

void HEGSolver::solve() {
  printf("Proc %d running on %s\n", Parallel::get_id(), Parallel::get_host().c_str());
  n_up = Config::get<size_t>("n_up");
  n_dn = Config::get<size_t>("n_dn");
  rcut_vars = Config::get_array<double>("rcut_vars");
  eps_vars = Config::get_array<double>("eps_vars");

  // Check configuration validity.
  check_validity();

  Time::start("variation");
  for (size_t i = 0; i < rcut_vars.size(); i++) {
    if (i > 0 && rcut_vars[i] == rcut_vars[i - 1]) continue;
    const double rcut_var = rcut_vars[i];
    std::string rcut_var_event = str(boost::format("rcut_var: %#.4g") % rcut_var);
    Time::start(rcut_var_event);

    Time::start("setup");
    setup(rcut_var);
    Time::end();

    for (size_t j = 0; j < eps_vars.size(); j++) {
      if (j > 0 && eps_vars[j] == eps_vars[j - 1]) continue;
      const double eps_var = eps_vars[j];
      const double eps_var_ham_old = eps_var * Config::get<double>("eps_var_ham_old_ratio");
      const double eps_var_ham_new = eps_var * Config::get<double>("eps_var_ham_new_ratio");
      std::string eps_var_event = str(boost::format("eps_var: %#.4g") % eps_var);
      Time::start(eps_var_event);

      // if (!load_variation_result()) {
      variation(eps_var, eps_var_ham_old, eps_var_ham_new);
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

void HEGSolver::check_validity() {
  if (Parallel::is_master()) {
    for (size_t i = 1; i < rcut_vars.size(); i++) {
      assert(rcut_vars[i - 1] <= rcut_vars[i]);
    }
    for (size_t i = 1; i < eps_vars.size(); i++) {
      assert(eps_vars[i - 1] >= eps_vars[i]);
    }
  }
  Parallel::barrier();
}

void HEGSolver::setup(const double rcut_var) {
  const double r_s = Config::get<double>("r_s");
  const double density = 3.0 / (4.0 * M_PI * pow(r_s, 3));
  const double cell_length = pow((n_up + n_dn) / density, 1.0 / 3);
  k_unit = 2 * M_PI / cell_length;
  H_unit = 1.0 / (M_PI * cell_length);
  k_points = KPointsUtil::generate_k_points(rcut_var);
  for (size_t i = 0; i < k_points.size(); i++) k_lut[k_points[i]] = i;
  if (Parallel::is_master()) {
    printf("number of orbitals: %d\n", static_cast<int>(k_points.size() * 2));
  }
  Time::start("hci_queue");
  generate_hci_queue(rcut_var);
  Time::end();
}

void HEGSolver::generate_hci_queue(const double rcut) {
  same_spin_hci_queue.clear();
  opposite_spin_hci_queue.clear();
  max_abs_H = 0.0;

  // Common dependencies.
  const auto& k_diffs = KPointsUtil::get_k_diffs(k_points);

  // Same spin.
  for (const auto& diff_pq : k_diffs) {
    for (const auto& diff_pr : k_diffs) {
      const auto& diff_sr = diff_pr + diff_pr - diff_pq;  // Momentum conservation.
      if (diff_sr == 0 || norm(diff_sr) > rcut * 2) continue;
      const auto& diff_ps = diff_pr - diff_sr;
      if (diff_ps == 0) continue;
      if (squared_norm(diff_pr) == squared_norm(diff_ps)) continue;
      const double abs_H = fabs(1.0 / squared_norm(diff_pr) - 1.0 / squared_norm(diff_ps));
      if (abs_H < DBL_EPSILON) continue;
      const auto& item = std::make_pair(diff_pr, abs_H * H_unit);
      same_spin_hci_queue[diff_pq].push_back(item);
    }
  }
  for (auto& kv : same_spin_hci_queue) {
    auto& items = kv.second;
    std::stable_sort(
        items.begin(),
        items.end(),
        [](const std::pair<std::array<int8_t, 3>, double>& a,
           const std::pair<std::array<int8_t, 3>, double>& b) -> bool {
          return a.second > b.second;
        });
    max_abs_H = std::max(max_abs_H, items.front().second);
  }

  // Opposite spin.
  for (const auto& diff_pr : k_diffs) {
    const double abs_H = 1.0 / sum(square(diff_pr));
    if (abs_H < DBL_EPSILON) continue;
    const auto& item = std::make_pair(diff_pr, abs_H * H_unit);
    opposite_spin_hci_queue.push_back(item);
  }
  std::stable_sort(
      opposite_spin_hci_queue.begin(),
      opposite_spin_hci_queue.end(),
      [](const std::pair<std::array<int8_t, 3>, double>& a,
         const std::pair<std::array<int8_t, 3>, double>& b) -> bool {
        return a.second > b.second;
      });
  max_abs_H = std::max(max_abs_H, opposite_spin_hci_queue.front().second);
}

double HEGSolver::hamiltonian(const Det& det_pq, const Det& det_rs) const {
  double H = 0.0;

  if (det_pq == det_rs) {
    const auto& occ_pq_up = det_pq.up.get_elec_orbs();
    const auto& occ_pq_dn = det_pq.dn.get_elec_orbs();

    // One electron operator.
    for (const auto p : occ_pq_up) H += squared_norm(k_points[p] * k_unit) * 0.5;
    for (const auto p : occ_pq_dn) H += squared_norm(k_points[p] * k_unit) * 0.5;

    // Two electrons operator.
    for (size_t i = 0; i < n_up; i++) {
      const auto p = occ_pq_up[i];
      for (size_t j = i + 1; j < n_up; j++) {
        const auto q = occ_pq_up[j];
        H -= H_unit / squared_norm(k_points[p] - k_points[q]);
      }
    }
    for (size_t i = 0; i < n_dn; i++) {
      const auto p = occ_pq_dn[i];
      for (size_t j = i + 1; j < n_dn; j++) {
        const auto q = occ_pq_dn[j];
        H -= H_unit / squared_norm(k_points[p] - k_points[q]);
      }
    }
  } else {
    // Off-diagonal elements.
    Det det_eor;
    det_eor.from_eor(det_pq, det_rs);
    const size_t n_eor_up = det_eor.up.get_n_elecs();
    const size_t n_eor_dn = det_eor.dn.get_n_elecs();
    if (n_eor_up + n_eor_dn != 4) return 0.0;
    const auto& eor_up_set_bits = det_eor.up.get_elec_orbs();
    const auto& eor_dn_set_bits = det_eor.dn.get_elec_orbs();
    bool k_p_set = false, k_r_set = false;
    uint16_t orb_p = 0, orb_r = 0, orb_s = 0;

    // Obtain p, q, s.
    std::array<int8_t, 3> k_change;
    k_change.fill(0);
    for (const auto orb_i : eor_up_set_bits) {
      if (det_pq.up.get_orb(orb_i)) {
        k_change -= k_points[orb_i];
        if (!k_p_set) {
          orb_p = orb_i;
          k_p_set = true;
        }
      } else {
        k_change += k_points[orb_i];
        if (!k_r_set) {
          orb_r = orb_i;
          k_r_set = true;
        } else {
          orb_s = orb_i;
        }
      }
    }
    for (const auto orb_i : eor_dn_set_bits) {
      if (det_pq.dn.get_orb(orb_i)) {
        k_change -= k_points[orb_i];
        if (!k_p_set) {
          orb_p = orb_i;
          k_p_set = true;
        }
      } else {
        k_change += k_points[orb_i];
        if (!k_r_set) {
          orb_r = orb_i;
          k_r_set = true;
        } else {
          orb_s = orb_i;
        }
      }
    }

    // Check for momentum conservation.
    if (k_change != 0) return 0.0;

    H = H_unit / squared_norm(k_points[orb_p] - k_points[orb_r]);
    if (n_eor_up != 2) H -= H_unit / squared_norm(k_points[orb_p] - k_points[orb_s]);

    const int gamma_exp =
        get_gamma_exp(det_pq.up, eor_up_set_bits) + get_gamma_exp(det_pq.dn, eor_dn_set_bits) +
        get_gamma_exp(det_rs.up, eor_up_set_bits) + get_gamma_exp(det_rs.dn, eor_dn_set_bits);
    if ((gamma_exp & 1) == 1) H = -H;
  }
  return H;
}

int HEGSolver::get_gamma_exp(const SpinDet& spin_det, const std::vector<uint16_t>& eor) const {
  int gamma_exp = 0;
  int ptr = 0;
  const auto& occ = spin_det.get_elec_orbs();
  for (const uint16_t orb_id : eor) {
    if (!spin_det.get_orb(orb_id)) continue;
    ptr = std::lower_bound(occ.begin() + ptr, occ.end(), orb_id) - occ.begin();
    gamma_exp += ptr;
  }
  return gamma_exp;
}

std::list<OrbitalPair> HEGSolver::get_pq_pairs(const Det& det, const int dn_offset) const {
  const auto& occ_up = det.up.get_elec_orbs();
  const auto& occ_dn = det.dn.get_elec_orbs();
  const size_t n_up = det.up.get_n_elecs();
  const size_t n_dn = det.dn.get_n_elecs();

  std::list<OrbitalPair> pq_pairs;

  for (size_t i = 0; i < n_up; i++) {
    for (size_t j = i + 1; j < n_up; j++) {
      pq_pairs.push_back(std::make_pair(occ_up[i], occ_up[j]));
    }
  }
  for (size_t i = 0; i < n_dn; i++) {
    for (size_t j = i + 1; j < n_dn; j++) {
      pq_pairs.push_back(std::make_pair(occ_dn[i] + dn_offset, occ_dn[j] + dn_offset));
    }
  }
  for (size_t i = 0; i < n_up; i++) {
    for (size_t j = 0; j < n_dn; j++) {
      pq_pairs.push_back(std::make_pair(occ_up[i], occ_dn[j] + dn_offset));
    }
  }

  return pq_pairs;
}

std::list<Det> HEGSolver::find_connected_dets(const Det& det, const double eps) const {
  std::list<Det> connected_dets;
  connected_dets.push_back(det);

  if (max_abs_H < eps) return connected_dets;

  const Orbital dn_offset = static_cast<Orbital>(k_points.size());
  const auto& pq_pairs = get_pq_pairs(det, dn_offset);

  for (const auto& pq_pair : pq_pairs) {
    const Orbital p = pq_pair.first;
    const Orbital q = pq_pair.second;

    // Get rs pairs.
    Orbital pp = p, qq = q;
    if (p >= dn_offset && q >= dn_offset) {
      pp -= dn_offset;
      qq -= dn_offset;
    } else if (p < dn_offset && q >= dn_offset && p > q - dn_offset) {
      pp = q - dn_offset;
      qq = p + dn_offset;
    }
    bool same_spin = false;
    std::vector<std::pair<std::array<int8_t, 3>, double>> const* items_ptr;
    if (pp < dn_offset && qq < dn_offset) {
      same_spin = true;
      const auto& diff_pq = k_points[qq] - k_points[pp];
      items_ptr = &(same_spin_hci_queue.find(diff_pq)->second);
    } else {
      items_ptr = &(opposite_spin_hci_queue);
    }
    const auto& items = *items_ptr;
    Orbital qs_offset = 0;
    if (!same_spin) qs_offset = dn_offset;

    for (const auto& item : items) {
      if (item.second < eps) break;
      const auto& diff_pr = item.first;
      const auto& it_r = k_lut.find(diff_pr + k_points[pp]);
      if (it_r == k_lut.end()) continue;
      Orbital r = it_r->second;
      const auto& it_s = k_lut.find(k_points[pp] + k_points[qq - qs_offset] - k_points[r]);
      if (it_s == k_lut.end()) continue;
      Orbital s = it_s->second;
      if (same_spin && s < r) continue;
      s += qs_offset;
      if (p >= dn_offset && q >= dn_offset) {
        r += dn_offset;
        s += dn_offset;
      } else if (p < dn_offset && q >= dn_offset && p > q - dn_offset) {
        const Orbital tmp = s;
        s = r + dn_offset;
        r = tmp - dn_offset;
      }

      // Test whether pqrs is a valid excitation for det.
      if (det.get_orb(r, dn_offset) || det.get_orb(s, dn_offset)) continue;
      connected_dets.push_back(det);
      Det& new_det = connected_dets.back();
      new_det.set_orb(p, dn_offset, false);
      new_det.set_orb(q, dn_offset, false);
      new_det.set_orb(r, dn_offset, true);
      new_det.set_orb(s, dn_offset, true);
    }
  }

  return connected_dets;
};