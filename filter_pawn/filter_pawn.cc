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
#define MAX_PAWN_COST_SUM 25
#define NUM_CAPTURABLE_CHESSMEN_PSIDE 15
#define NUM_SQUARES (BOARD_SIDE_LENGTH * BOARD_SIDE_LENGTH)

namespace operations_research {
namespace sat {

const Domain BOARD_SIDE_LENGTH_DOMAIN(0, BOARD_SIDE_LENGTH);
const Domain PAWN_COST_SUM_DOMAIN(0, MAX_PAWN_COST_SUM);
const Domain SINGLE_PAWN_COST_DOMAIN(-BOARD_SIDE_LENGTH, BOARD_SIDE_LENGTH);
const Domain ZERO_DOMAIN(0, 0);

extern "C" {
int get_row_num(int n) { return (n / BOARD_SIDE_LENGTH); }
int get_col_num(int n) { return (n % BOARD_SIDE_LENGTH); }
}

extern "C" bool FilterPawn(uint64_t _pawns[NUM_SIDES], uint64_t enpassant,
                           int captured_base_pieces[NUM_SIDES],
                           int min_promotions_to_acc_for[NUM_SIDES]) {
  CpModelBuilder cp_model;

  BoolVar matchedStartingSquares[NUM_SIDES][BOARD_SIDE_LENGTH];
  BoolVar pawnWhichStartedOn_isConsumedByPawnOnTheBoard[NUM_SIDES]
                                                       [BOARD_SIDE_LENGTH];
  BoolVar matchedToSameFileAsPawn[NUM_SIDES][BOARD_SIDE_LENGTH];
  BoolVar sameFilePromotions[NUM_SIDES][BOARD_SIDE_LENGTH];
  BoolVar capturedBasePiecePromotions[NUM_SIDES][BOARD_SIDE_LENGTH];
  for (int i = 0; i < NUM_SIDES; i++) {
    for (int j = 0; j < BOARD_SIDE_LENGTH; j++) {
      matchedStartingSquares[i][j] = cp_model.NewBoolVar();
      pawnWhichStartedOn_isConsumedByPawnOnTheBoard[i][j] =
          cp_model.NewBoolVar();
      matchedToSameFileAsPawn[i][j] = cp_model.NewBoolVar();
      sameFilePromotions[i][j] = cp_model.NewBoolVar();
      capturedBasePiecePromotions[i][j] = cp_model.NewBoolVar();
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
        for (int k = start; k < col; k++) {
          BoolVar p_is_matched_to_k = cp_model.NewBoolVar();
          cp_model.AddEquality(pawnVars[i][p], k)
              .OnlyEnforceIf(p_is_matched_to_k);
          cp_model.AddImplication(p_is_matched_to_k,
                                  matchedStartingSquares[i][k]);

          BoolVar k_is_covered_by_p = cp_model.NewBoolVar();
          cp_model.AddLessOrEqual(pawnVars[i][p], k)
              .OnlyEnforceIf(k_is_covered_by_p);

          BoolVar covered_but_not_matched_literals[] = {
              k_is_covered_by_p, p_is_matched_to_k.Not()};
          cp_model.AddBoolAnd(covered_but_not_matched_literals)
              .OnlyEnforceIf(
                  pawnWhichStartedOn_isConsumedByPawnOnTheBoard[opp][k]);
        }

        BoolVar p_is_matched_to_col = cp_model.NewBoolVar();
        cp_model.AddEquality(pawnVars[i][p], col)
            .OnlyEnforceIf(p_is_matched_to_col);
        cp_model.AddImplication(p_is_matched_to_col,
                                matchedStartingSquares[i][col]);
        cp_model.AddImplication(p_is_matched_to_col,
                                matchedToSameFileAsPawn[i][col]);

        for (int k = col + 1; k < end + 1; k++) {
          BoolVar p_is_matched_to_k = cp_model.NewBoolVar();
          cp_model.AddEquality(pawnVars[i][p], k)
              .OnlyEnforceIf(p_is_matched_to_k);
          cp_model.AddImplication(p_is_matched_to_k,
                                  matchedStartingSquares[i][k]);

          BoolVar k_is_covered_by_p = cp_model.NewBoolVar();
          cp_model.AddGreaterOrEqual(pawnVars[i][p], k)
              .OnlyEnforceIf(k_is_covered_by_p);

          BoolVar covered_but_not_matched_literals[2] = {
              k_is_covered_by_p, p_is_matched_to_k.Not()};
          cp_model.AddBoolAnd(covered_but_not_matched_literals)
              .OnlyEnforceIf(
                  pawnWhichStartedOn_isConsumedByPawnOnTheBoard[opp][k]);
        }

        pawnCosts[i].push_back(cp_model.NewIntVar(SINGLE_PAWN_COST_DOMAIN));
        cp_model.AddAbsEquality(pawnCosts[i][p], pawnVars[i][p] - p);
      }
    }
    cp_model.AddAllDifferent(pawnVars[i]);
  }

  BoolVar pawnWhichStartedOn_isPaired[2][NUM_SIDES][BOARD_SIDE_LENGTH];
  for (int i = 0; i < 2; i++) { // left or right
    for (int j = 0; j < NUM_SIDES; j++) {
      for (int k = 0; k < BOARD_SIDE_LENGTH; k++) {
        pawnWhichStartedOn_isPaired[i][j][k] = cp_model.NewBoolVar();
      }
    }
  }
  int d = 0;
  int other_direction = 1;
  for (int i = -1; i <= 1; i += 2) {
    for (int j = 1; j < BOARD_SIDE_LENGTH; j++) {
      cp_model.AddEquality(pawnWhichStartedOn_isPaired[d][0][j],
                           pawnWhichStartedOn_isPaired[d][1][j + i]);
      BoolVar literals[] = {
          matchedStartingSquares[0][j].Not(),
          matchedStartingSquares[1][j + i].Not(),
          pawnWhichStartedOn_isConsumedByPawnOnTheBoard[0][j].Not(),
          pawnWhichStartedOn_isConsumedByPawnOnTheBoard[1][j + i].Not(),
          pawnWhichStartedOn_isPaired[other_direction][0][j].Not(),
          pawnWhichStartedOn_isPaired[other_direction][1][j + i].Not()};
      cp_model.AddBoolAnd(literals).OnlyEnforceIf(
          pawnWhichStartedOn_isPaired[d][0][j]);
    }
    d = 1;
    other_direction = 0;
  }
  IntVar v[] = {cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN),
                cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN)};
  cp_model.AddEquality(v[0],
                       LinearExpr::Sum(pawnWhichStartedOn_isPaired[0][0]));
  cp_model.AddEquality(v[1],
                       LinearExpr::Sum(pawnWhichStartedOn_isPaired[1][0]));
  IntVar vsum = cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN);
  cp_model.AddEquality(vsum, LinearExpr::Sum(v));
  IntVar pairsGoTo[] = {cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN),
                        cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN)};
  cp_model.AddLessOrEqual(pairsGoTo[0], vsum);
  cp_model.AddEquality(pairsGoTo[1], BOARD_SIDE_LENGTH - pairsGoTo[0]);

  for (int i = 0; i < NUM_SIDES; i++) {
    int opp = (i + 1) % NUM_SIDES;
    for (int j = 0; j < BOARD_SIDE_LENGTH; j++) {
      BoolVar literals[] = {
          matchedStartingSquares[i][j].Not(),
          pawnWhichStartedOn_isConsumedByPawnOnTheBoard[i][j].Not(),
          pawnWhichStartedOn_isPaired[0][i][j].Not(),
          pawnWhichStartedOn_isPaired[1][i][j].Not(),
          matchedToSameFileAsPawn[opp][j].Not(),
          sameFilePromotions[opp][j].Not(),
      };
      cp_model.AddBoolAnd(literals).OnlyEnforceIf(sameFilePromotions[i][j]);
    }
  }

  for (int i = 0; i < NUM_SIDES; i++) {
    for (int j = 0; j < BOARD_SIDE_LENGTH; j++) {
      BoolVar literals[] = {
          matchedStartingSquares[i][j].Not(),
          pawnWhichStartedOn_isConsumedByPawnOnTheBoard[i][j].Not(),
          pawnWhichStartedOn_isPaired[0][i][j].Not(),
          pawnWhichStartedOn_isPaired[1][i][j].Not(),
          sameFilePromotions[i][j].Not(),
      };
      cp_model.AddBoolAnd(literals).OnlyEnforceIf(
          capturedBasePiecePromotions[i][j]);
    }
  }

  for (int i = 0; i < NUM_SIDES; i++) {
    IntVar pcSum = cp_model.NewIntVar(PAWN_COST_SUM_DOMAIN);
    cp_model.AddEquality(pcSum, LinearExpr::Sum(pawnCosts[i]));

    IntVar capturedBasePromSum = cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN);
    cp_model.AddEquality(capturedBasePromSum,
                         LinearExpr::Sum(capturedBasePiecePromotions[i]));

    IntVar coveredSum = cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN);
    int opp = (i + 1) % NUM_SIDES;
    cp_model.AddEquality(
        coveredSum,
        LinearExpr::Sum(pawnWhichStartedOn_isConsumedByPawnOnTheBoard[opp]));

    cp_model.AddLessOrEqual(pcSum + capturedBasePromSum - coveredSum,
                            captured_base_pieces[opp]);

    IntVar sameFilePromSum = cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN);
    cp_model.AddEquality(sameFilePromSum,
                         LinearExpr::Sum(sameFilePromotions[i]));

    cp_model.AddGreaterOrEqual(sameFilePromSum + capturedBasePromSum,
                               min_promotions_to_acc_for[i]);
  }

  const CpSolverResponse response = Solve(cp_model.Build());
  if (response.status() == CpSolverStatus::OPTIMAL ||
      response.status() == CpSolverStatus::FEASIBLE) {
    return true;
  }
  return false;
} // namespace sat

} // namespace sat
} // namespace operations_research
