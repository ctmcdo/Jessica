#include <immintrin.h>
#include <inttypes.h>
#include <nmmintrin.h>
#include <stdint.h>
#include <stdlib.h>

#include <algorithm>

#include "filter_pawn.h"
#include "ortools/base/logging.h"
#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/util/sorted_interval_list.h"

#define BOARD_SIDE_LENGTH 8
#define NUM_CAPTURABLE_CHESSMEN_PSIDE 15
#define NUM_SQUARES (BOARD_SIDE_LENGTH * BOARD_SIDE_LENGTH)

namespace operations_research {
namespace sat {

const Domain ZERO_DOMAIN(0, 0);
const Domain COST_DOMAIN(-BOARD_SIDE_LENGTH, BOARD_SIDE_LENGTH);

extern "C" {
int get_row_num(int n) { return (n / BOARD_SIDE_LENGTH); }
int get_col_num(int n) { return (n % BOARD_SIDE_LENGTH); }
}

extern "C" bool FilterPawn(uint64_t _pawns[NUM_SIDES], uint64_t enpassant,
                           int num_capturable_chessmen[NUM_SIDES],
                           int min_promotions[NUM_SIDES]) {
  CpModelBuilder cp_model;

  std::vector<std::vector<IntVar>> pawnVars(NUM_SIDES);
  std::vector<std::vector<IntVar>> pawnCosts(NUM_SIDES);
  uint64_t pawns[] = {_pawns[0], _pawns[1]};
  for (int i = 0; i < NUM_SIDES; i++) {
    int num_pawns = _mm_popcnt_u64(pawns[i]);
    pawnVars[i].reserve(num_pawns);
    pawnCosts[i].reserve(num_pawns);
    for (int p = 0; p < num_pawns; p++) {
      int tz = __tzcnt_u64(pawns[i]);
      // std::cout << "tz: " << tz << "\n";
      pawns[i] = __blsr_u64(pawns[i]);

      int col = get_col_num(tz);
      if (i == 0 && enpassant == (1UL << tz)) {
        const Domain single_col_domain(col, col);
        pawnVars[i].push_back(cp_model.NewIntVar(single_col_domain));
        pawnCosts[i].push_back(cp_model.NewIntVar(ZERO_DOMAIN));
      } else {
        int row = get_row_num(tz);

        int start = std::max(col - row, 0);
        int end = std::min(col + row, BOARD_SIDE_LENGTH - 1);
        const Domain d(start, end);
        pawnVars[i].push_back(cp_model.NewIntVar(d));
        pawnCosts[i].push_back(cp_model.NewIntVar(COST_DOMAIN));
        cp_model.AddAbsEquality(pawnCosts[i][p], pawnVars[i][p] - p);
        // if (is_kth_bit_set(mask, tz)) {
        //  if (tz >= (2 * BOARD_SIDE_LENGTH)) {
        //    vertices[v][col] = 2;
        //  } else {
        //    vertices[v][col] = INFIN;
        //  }
        //}
        //}
      }
    }

    cp_model.AddAllDifferent(pawnVars[i]);
    cp_model.AddLessOrEqual(LinearExpr::Sum(pawnCosts[i]),
                            NUM_CAPTURABLE_CHESSMEN_PSIDE -
                                num_capturable_chessmen[i]);
  }

  const CpSolverResponse response = Solve(cp_model.Build());
  if (response.status() == CpSolverStatus::OPTIMAL ||
      response.status() == CpSolverStatus::FEASIBLE) {
    return true;
  }
  return false;
}

} // namespace sat
} // namespace operations_research
