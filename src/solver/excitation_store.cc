#include "excitation_store.h"

void ExcitationStore::add(const Orbitals& orbs_1, const Orbitals& orbs_2, const bool same_spin) {
  const uint32_t id_1 = get_id(orbs_1);
  const uint32_t id_2 = get_id(orbs_2);
  if (same_spin) {
    same_spin_excitations[id_1].insert(id_2);
    same_spin_excitations[id_2].insert(id_1);
  } else {
    opposite_spin_excitations[id_1].insert(id_2);
    opposite_spin_excitations[id_2].insert(id_1);
  }
}

std::vector<uint32_t> ExcitationStore::find(const Orbitals& orbs, const bool same_spin) const {
  std::vector<uint32_t> res;
  if (lut.count(orbs) == 0) return res;
  const uint32_t id = lut.at(orbs);
  res = find(id, same_spin);
  return res;
}

std::vector<uint32_t> ExcitationStore::find(const uint32_t orbs_id, const bool same_spin) const {
  std::vector<uint32_t> res;
  if (same_spin) {
    if (same_spin_excitations.count(orbs_id) == 0) return res;
    for (const uint32_t ex_id : same_spin_excitations.at(orbs_id)) {
      res.push_back(ex_id);
    }
  } else {
    if (opposite_spin_excitations.count(orbs_id) == 0) return res;
    for (const uint32_t ex_id : opposite_spin_excitations.at(orbs_id)) {
      res.push_back(ex_id);
    }
  }
  return res;
}

uint32_t ExcitationStore::get_id(const Orbitals& orbs) {
  if (lut.count(orbs) == 0) {
    const uint32_t id = unique_orbs.size();
    unique_orbs.push_back(orbs);
    lut[orbs] = id;
    return id;
  } else {
    return lut[orbs];
  }
}