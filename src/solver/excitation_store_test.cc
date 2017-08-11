#include "excitation_store.h"
#include "gtest/gtest.h"

TEST(ExcitationStoreTest, AddAndFind) {
  ExcitationStore ex;
  Orbitals orbs_1({1, 2, 3});
  Orbitals orbs_2({1, 2, 4});
  Orbitals orbs_3({1, 5, 6});
  ex.add(orbs_1, orbs_2, false);
  EXPECT_EQ(ex.find(orbs_1, true).size(), 0);
  EXPECT_EQ(ex.find(orbs_1, false).size(), 1);
  EXPECT_EQ(ex.find(orbs_1, false)[0], ex.lookup_id(orbs_2));

  ex.add(orbs_1, orbs_3, false);
  EXPECT_EQ(ex.find(orbs_1, true).size(), 0);
  EXPECT_EQ(ex.find(orbs_1, false).size(), 2);

  ex.add(orbs_1, orbs_3, true);
  EXPECT_EQ(ex.find(orbs_1, true).size(), 1);
  EXPECT_EQ(ex.find(orbs_1, true)[0], ex.lookup_id(orbs_3));
  EXPECT_EQ(ex.find(orbs_1, false).size(), 2);
}