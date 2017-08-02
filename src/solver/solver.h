#ifndef SOLVER_H_
#define SOLVER_H_

class Solver {
 protected:
  std::size_t n_up;
  std::size_t n_dn;

  virtual void solve() {}
};

#endif