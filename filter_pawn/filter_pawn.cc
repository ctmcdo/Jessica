#include <immintrin.h>
#include <inttypes.h>
#include <nmmintrin.h>
#include <stdint.h>
#include <stdlib.h>

#include <algorithm>
#include <vector>

#include "filter_pawn.h"
#include "ortools/base/logging.h"
#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/util/sorted_interval_list.h"

using namespace std;

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

  vector<IntVar> pawnVars[NUM_SIDES];
  vector<IntVar> pawnCosts[NUM_SIDES];

  BoolVar matchedStartingSquares[NUM_SIDES][BOARD_SIDE_LENGTH];
  BoolVar startingSquareIsCovered[NUM_SIDES][BOARD_SIDE_LENGTH];
  BoolVar matchedToSameFileAsPawn[NUM_SIDES][BOARD_SIDE_LENGTH];
  BoolVar sameFilePromotions[NUM_SIDES][BOARD_SIDE_LENGTH];
  for (int i = 0; i < NUM_SIDES; i++) {
    for (int j = 0; j < BOARD_SIDE_LENGTH; j++) {
      matchedStartingSquares[i][j] = cp_model.NewBoolVar();
      startingSquareIsCovered[i][j] = cp_model.NewBoolVar();
      matchedToSameFileAsPawn[i][j] = cp_model.NewBoolVar();
      sameFilePromotions[i][j] = cp_model.NewBoolVar();
    }
  }

  vector<BoolVar> fileMatchOptionsNegated[NUM_SIDES][BOARD_SIDE_LENGTH];
  vector<BoolVar> fileCoverOptionsNegated[NUM_SIDES][BOARD_SIDE_LENGTH];

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
        cp_model.FixVariable(matchedStartingSquares[i][col], true);
        cp_model.FixVariable(matchedToSameFileAsPawn[i][col], true);

      } else {
        int row = get_row_num(tz);
        int start = std::max(col - row, 0);
        int end = std::min(col + row, BOARD_SIDE_LENGTH - 1);
        const Domain d(start, end);
        pawnVars[i].push_back(cp_model.NewIntVar(d));
        for (int k = start; k < col; k++) {
          BoolVar pIsMatchedToK = cp_model.NewBoolVar();
          cp_model.AddEquality(pawnVars[i][p], k).OnlyEnforceIf(pIsMatchedToK);
          cp_model.AddNotEqual(pawnVars[i][p], k).OnlyEnforceIf(pIsMatchedToK);
          cp_model.AddImplication(pIsMatchedToK, matchedStartingSquares[i][k]);
          fileMatchOptionsNegated[i][k].push_back(pIsMatchedToK.Not());

          BoolVar pCoversK = cp_model.NewBoolVar();
          cp_model.AddLessOrEqual(pawnVars[i][p], k).OnlyEnforceIf(pCoversK);
          cp_model.AddGreaterThan(pawnVars[i][p], k)
              .OnlyEnforceIf(pCoversK.Not());
          cp_model.AddImplication(pCoversK, startingSquareIsCovered[opp][k]);
          fileCoverOptionsNegated[i][k].push_back(pCoversK.Not());
        }

        BoolVar pIsMatchedToFile = cp_model.NewBoolVar();
        cp_model.AddEquality(pawnVars[i][p], col)
            .OnlyEnforceIf(pIsMatchedToFile);
        cp_model.AddNotEqual(pawnVars[i][p], col)
            .OnlyEnforceIf(pIsMatchedToFile.Not());
        cp_model.AddImplication(pIsMatchedToFile,
                                matchedStartingSquares[i][col]);
        cp_model.AddImplication(pIsMatchedToFile,
                                matchedToSameFileAsPawn[i][col]);
        fileMatchOptionsNegated[i][col].push_back(pIsMatchedToFile);

        for (int k = col + 1; k <= end; k++) {
          BoolVar pIsMatchedToK = cp_model.NewBoolVar();
          cp_model.AddEquality(pawnVars[i][p], k).OnlyEnforceIf(pIsMatchedToK);
          cp_model.AddNotEqual(pawnVars[i][p], k).OnlyEnforceIf(pIsMatchedToK);
          cp_model.AddImplication(pIsMatchedToK, matchedStartingSquares[i][k]);
          fileMatchOptionsNegated[i][k].push_back(pIsMatchedToK.Not());

          BoolVar pCoversK = cp_model.NewBoolVar();
          cp_model.AddGreaterOrEqual(pawnVars[i][p], k).OnlyEnforceIf(pCoversK);
          cp_model.AddLessThan(pawnVars[i][p], k).OnlyEnforceIf(pCoversK.Not());
          cp_model.AddImplication(pCoversK, startingSquareIsCovered[opp][k]);
          fileCoverOptionsNegated[i][k].push_back(pCoversK.Not());
        }

        pawnCosts[i].push_back(cp_model.NewIntVar(SINGLE_PAWN_COST_DOMAIN));
        cp_model.AddAbsEquality(pawnCosts[i][p], pawnVars[i][p] - p);
      }
    }
    cp_model.AddAllDifferent(pawnVars[i]);
  }

  for (int i = 0; i < NUM_SIDES; i++) {
    for (int j = 0; j < BOARD_SIDE_LENGTH; j++) {
      cp_model.AddBoolAnd(fileMatchOptionsNegated[i][j])
          .OnlyEnforceIf(matchedStartingSquares[i][j].Not());

      cp_model.AddBoolAnd(fileCoverOptionsNegated[i][j])
          .OnlyEnforceIf(startingSquareIsCovered[i][j].Not());
    }
  }

  BoolVar pawnWhichStartedOn_isConsumedByPawnOnBoard[NUM_SIDES]
                                                    [BOARD_SIDE_LENGTH];
  for (int i = 0; i < NUM_SIDES; i++) {
    int opp = (i + 1) % NUM_SIDES;
    for (int j = 0; j < BOARD_SIDE_LENGTH; j++) {
      BoolVar andLiterals[] = {matchedStartingSquares[i][j].Not(),
                               startingSquareIsCovered[i][j]};
      pawnWhichStartedOn_isConsumedByPawnOnBoard[i][j] = cp_model.NewBoolVar();
      cp_model.AddBoolAnd(andLiterals)
          .OnlyEnforceIf(pawnWhichStartedOn_isConsumedByPawnOnBoard[i][j]);
    }
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
  int col_lim = BOARD_SIDE_LENGTH;
  for (int i = -1; i <= 1; i += 2) {
    for (int j = 1; j < col_lim; j++) {
      cp_model.AddEquality(pawnWhichStartedOn_isPaired[d][0][j],
                           pawnWhichStartedOn_isPaired[d][1][j + i]);
      BoolVar andLiterals[] = {
          matchedStartingSquares[0][j].Not(),
          matchedStartingSquares[1][j + i].Not(),
          pawnWhichStartedOn_isConsumedByPawnOnBoard[0][j].Not(),
          pawnWhichStartedOn_isConsumedByPawnOnBoard[1][j + i].Not(),
          pawnWhichStartedOn_isPaired[other_direction][0][j].Not(),
          pawnWhichStartedOn_isPaired[other_direction][1][j + i].Not()};
      cp_model.AddBoolAnd(andLiterals)
          .OnlyEnforceIf(pawnWhichStartedOn_isPaired[d][0][j]);
    }
    d = 1;
    other_direction = 0;
    col_lim = BOARD_SIDE_LENGTH - 1;
  }
  IntVar num_pairs_dir0 = cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN);
  cp_model.AddEquality(num_pairs_dir0,
                       LinearExpr::Sum(pawnWhichStartedOn_isPaired[0][0]));
  IntVar num_pairs_dir1 = cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN);
  cp_model.AddEquality(num_pairs_dir1,
                       LinearExpr::Sum(pawnWhichStartedOn_isPaired[1][0]));
  IntVar num_pairs = cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN);
  cp_model.AddEquality(num_pairs, num_pairs_dir0 + num_pairs_dir1);
  IntVar pairsGoTo[] = {cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN),
                        cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN)};
  cp_model.AddLessOrEqual(pairsGoTo[0], num_pairs);
  cp_model.AddEquality(pairsGoTo[1], num_pairs - pairsGoTo[0]);

  for (int i = 0; i < NUM_SIDES; i++) {
    int opp = (i + 1) % NUM_SIDES;
    for (int j = 0; j < BOARD_SIDE_LENGTH; j++) {
      BoolVar andLiterals[] = {
          matchedStartingSquares[i][j].Not(),
          pawnWhichStartedOn_isConsumedByPawnOnBoard[i][j].Not(),
          pawnWhichStartedOn_isPaired[0][i][j].Not(),
          pawnWhichStartedOn_isPaired[1][i][j].Not(),
          matchedToSameFileAsPawn[opp][j].Not(),
          sameFilePromotions[opp][j].Not(),
      };
      cp_model.AddBoolAnd(andLiterals).OnlyEnforceIf(sameFilePromotions[i][j]);
    }
  }

  BoolVar capturedBasePiecePromotions[NUM_SIDES][BOARD_SIDE_LENGTH];
  for (int i = 0; i < NUM_SIDES; i++) {
    for (int j = 0; j < BOARD_SIDE_LENGTH; j++) {
      BoolVar andLiterals[] = {
          matchedStartingSquares[i][j].Not(),
          pawnWhichStartedOn_isConsumedByPawnOnBoard[i][j].Not(),
          pawnWhichStartedOn_isPaired[0][i][j].Not(),
          pawnWhichStartedOn_isPaired[1][i][j].Not(),
          sameFilePromotions[i][j].Not(),
      };
      capturedBasePiecePromotions[i][j] = cp_model.NewBoolVar();
      cp_model.AddBoolAnd(andLiterals)
          .OnlyEnforceIf(capturedBasePiecePromotions[i][j]);
    }
  }

  for (int i = 0; i < NUM_SIDES; i++) {
    int opp = (i + 1) % NUM_SIDES;

    IntVar pcSum = cp_model.NewIntVar(PAWN_COST_SUM_DOMAIN);
    cp_model.AddEquality(pcSum, LinearExpr::Sum(pawnCosts[i]));

    IntVar capturedBasePromSum = cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN);
    cp_model.AddEquality(capturedBasePromSum,
                         LinearExpr::Sum(capturedBasePiecePromotions[i]));

    IntVar consumedSum = cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN);
    cp_model.AddEquality(
        consumedSum,
        LinearExpr::Sum(pawnWhichStartedOn_isConsumedByPawnOnBoard[opp]));
    cp_model.AddLessOrEqual(pcSum + capturedBasePromSum - consumedSum,
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
