#ifndef EXCITATION_H_
#define EXCITATION_H_

#include <boost/functional/hash.hpp>
#include "../std.h"
#include "../wavefunction/spin_det.h"

class ExcitationStore {
 public:
  void add(const Orbitals&, const Orbitals&, const bool same_spin);

  std::vector<uint32_t> find(const Orbitals&, const bool same_spin) const;

  std::vector<uint32_t> find(const uint32_t, const bool same_spin) const;

  uint32_t lookup_id(const Orbitals& orbs) { return lut[orbs]; }

 private:
  std::vector<Orbitals> unique_orbs;
  std::unordered_map<Orbitals, uint32_t, boost::hash<Orbitals>> lut;
  std::unordered_map<uint32_t, std::unordered_set<uint32_t>> same_spin_excitations;
  std::unordered_map<uint32_t, std::unordered_set<uint32_t>> opposite_spin_excitations;

  uint32_t get_id(const Orbitals&);
};

#endif