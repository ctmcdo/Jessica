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
                           int min_promotions_to_acc_for[NUM_SIDES]) {
  CpModelBuilder cp_model;

  BoolVar matchedStartingSquares[NUM_SIDES][BOARD_SIDE_LENGTH];
  BoolVar pawnWhichStartedOn_isConsumedByPawnOnTheBoard[NUM_SIDES]
                                                       [BOARD_SIDE_LENGTH];
  for (int i = 0; i < NUM_SIDES; i++) {
    for (int j = 0; j < BOARD_SIDE_LENGTH; j++) {
      matchedStartingSquares[i][j] = cp_model.NewBoolVar();
      pawnWhichStartedOn_isConsumedByPawnOnTheBoard[i][j] =
          cp_model.NewBoolVar();
    }
  }

  std::vector<std::vector<IntVar>> pawnVars(NUM_SIDES);
  std::vector<std::vector<IntVar>> pawnCosts(NUM_SIDES);
  uint64_t pawns[] = {_pawns[0], _pawns[1]};
  for (int i = 0; i < NUM_SIDES; i++) {
    int opp = (i + 1) % NUM_SIDES;
    int num_pawns = _mm_popcnt_u64(pawns[i]);
    pawnVars[i].reserve(num_pawns);
    pawnCosts[i].reserve(num_pawns);
    for (int p = 0; p < num_pawns; p++) {
      int tz = __tzcnt_u64(pawns[i]);
      pawns[i] = __blsr_u64(pawns[i]);

      int col = get_col_num(tz);
      if (i == 0 && enpassant == (1UL << tz)) {
        const Domain single_col_domain(col, col);
        pawnVars[i].push_back(cp_model.NewIntVar(single_col_domain));
        pawnCosts[i].push_back(cp_model.NewIntVar(ZERO_DOMAIN));
        cp_model.AddEquality(matchedStartingSquares[i][col], true);

      } else {
        int row = get_row_num(tz);
        int start = std::max(col - row, 0);
        int end = std::min(col + row, BOARD_SIDE_LENGTH - 1);
        const Domain d(start, end);
        pawnVars[i].push_back(cp_model.NewIntVar(d));
        for (int k = start; k < end + 1; k++) {
          BoolVar p_is_matched_to_k = cp_model.NewBoolVar();
          cp_model.AddEquality(pawnVars[i][p], k)
              .onlyEnforceIf(p_is_matched_to_k);
          cp_model.AddImplication(p_is_matched_to_k,
                                  matchedStartingSquares[i][k]);

          BoolVar k_is_covered_by_p = cp_model.NewBoolVar();
          if (k < col) {
            cp_model.AddLessOrEqual(pawnVars[i][p], k)
                .onlyEnfoceIf(k_is_covered_by_p);
          } else {
            cp_model.AddGreaterOrEqual(pawnVars[i][p], k)
                .onlyEnfoceIf(k_is_covered_by_p);
          }
          BoolVar covered_but_not_matched_literals[2] = {
              k_is_covered_by_p, p_is_matched_to_k.Not()};
          cp_model.AddBoolAnd(covered_but_not_matched_literals)
              .onlyEnforceIf(
                  pawnWhichStartedOn_isConsumedByPawnOnTheBoard[opp][k]);
        }

        pawnCosts[i].push_back(cp_model.NewIntVar(COST_DOMAIN));
        cp_model.AddAbsEquality(pawnCosts[i][p], pawnVars[i][p] - p);
      }
    }

    std::vector<std::vector<std::vector<BoolVar>>> pawnWhichStartedOn_isPaired(
        NUM_SIDES);
    int d = 0;
    int other_direction = 1;
    for (int i = -1; i <= 1; i += 2) {
      for (int j = 1; j < BOARD_SIDE_LENGTH; j++) {
        cp_model.AddEquality(pawnWhichStartedOn_isPaired[d][0][j],
                             pawnWhichStartedOn_isPaired[d][1][j + i]);
        BoolVar literals[] = {
            matchedStartingSquare[0][j].Not(),
            matchedStartingSquare[1][j + i].Not(),
            pawnWhichStartedOn_isConsumedByPawnOnTheBoard[0][j].Not(),
            pawnWhichStartedOn_isConsumedByPawnOnTheBoard[1][j + i].Not(),
            pawnWhichStartedOn_isPaired[other_direction][0][j].Not(),
            pawnWhichStartedOn_isPaired[other_direction][1][j + i].Not()};
        cp.model.AddBoolAnd(literals).OnlyEnforceIf(
            pawnWhichStartedOn_isPaired[d][0][j]);
      }
      d = 1;
      other_direction = 0;
    }

    std::vector<std::vector<BoolVar>> freePromotions(NUM_SIDES);
    for (int i = 0; i < NUM_SIDES; i++) {
      int opp = (i + 1) % NUM_SIDES;
      for (int j = 0; j < BOARD_SIDE_LENGTH; j++) {
        BoolVar opp_never_left_file = cp_model.NewBoolVar();
        cp_model.Add(pawnVar[opp][j] == j).OnlyEnforceIf(opp_never_left_file);
        BoolVar literals[] = {
            matchedStartingSquare[i][j].Not(),
            pawnWhichStartedOn_isConsumedByPawnOnTheBoard[i][j].Not(),
            pawnWhichStartedOn_isPaired[0][i][j].Not(),
            pawnWhichStartedOn_isPaired[1][i][j].Not(),
            opp_never_left_file.Not()};
        cp_model.AddBoolAnd(literals).OnlyEnforceIf(freePromotions[i][j]);
      }
    }

    cp_model.AddImplication();
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
