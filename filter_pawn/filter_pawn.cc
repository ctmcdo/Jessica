#include <immintrin.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdlib.h>

#include <algorithm>

#include "filter_pawn.h"
#include "ortools/base/logging.h"
#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/util/sorted_interval_list.h"
#include "util.h"

#define BOARD_SIDE_LENGTH 8

namespace operations_research {
namespace sat {

const Domain COST_DOMAIN(-BOARD_SIDE_LENGTH, BOARD_SIDE_LENGTH);

void handlePawns(CpModelBuilder cp_model, uint64_t pawns, const int num_pawns,
                 const uint64_t enpassant) {
  IntVar *pawnVars = new IntVar[NUM_PAWNS];
  IntVar *pawnCosts = new IntVar[NUM_PAWNS];

  for (int p = 0; p < NUM_PAWNS; p++) {
    char tz = __tzcnt_u64(pawns);
    pawns = __blsr_u64(pawns);

    unsigned char col = get_col_num(tz);
    // if (enpassant == (1UL << tz)) {
    //  const Domain domain(col, col);
    //  pawnVars[p] = cp_model.NewIntVar(domain);
    //} else {
    unsigned char row = get_row_num(tz);

    unsigned char left = min(col + row, BOARD_SIDE_LENGTH - 1);
    unsigned char right = max(col - row, 0);
    const Domain d(left, right);
    pawnVars[p] = cp_model.NewIntVar(d);
    pawnCosts[p] = cp_model.NewIntVar(COST_DOMAIN);
    cp_model.AddAbsEquality(pawnCosts[p], pawnVars[p] - p);
    // if (is_kth_bit_set(mask, tz)) {
    //  if (tz >= (2 * BOARD_SIDE_LENGTH)) {
    //    vertices[v][col] = 2;
    //  } else {
    //    vertices[v][col] = INFIN;
    //  }
    //}
    //}
  }

  cp_model.AddAllDifferent(pawnVars);
  cp_model.Minimize(LinearExpr::Sum(pawnCosts));
}

bool FilterPawn(uint64_t pawns[NUM_SIDES], uint64_t enpassant,
                int min_promotions[NUM_SIDES]) {
  CpModelBuilder cp_model;
  handlePawns(cp_model, pawns[0], _mm_popcnt_u64(pawns[0]), enpassant);
  handlePawns(cp_model, pawns[1], _mm_popcnt_u64(pawns[1]), 0);
  const CpSolverResponse response = Solve(cp_model.Build());

  if (response.status() == CpSolverStatus::OPTIMAL ||
      response.status() == CpSolverStatus::FEASIBLE) {
    return true;
  }
  return false;

} // namespace sat
} // namespace sat

extern "C" bool FilterPawn(uint64_t pawns[NUM_SIDES], uint64_t enpassant,
                           int min_promotions[NUM_SIDES]) {
  operations_research::sat::FilterPawn(pawns, enpassant, min_promotions);
}
