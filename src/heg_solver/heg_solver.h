#ifndef HEG_SOLVER_H_
#define HEG_SOLVER_H_

#include "../libs.h"
#include "../solver/solver.h"

class HEGSolver : public Solver {
 public:
  static void run() { HEGSolver::get_instance().solve(); }

  void solve() override;

  void setup();

 private:
  std::vector<double> rcut_vars;
  std::vector<double> eps_vars;
  std::vector<double> rcut_pts;
  std::vector<double> eps_pts;

  static HEGSolver get_instance() {
    static HEGSolver heg_solver;
    return heg_solver;
  }
};

#endif