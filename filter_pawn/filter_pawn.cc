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
#define SQUARES_ON_PAWN_FILE (BOARD_SIDE_LENGTH - 2)
#define NUM_CAPTURABLE_CHESSMEN_PSIDE 15
#define NUM_SQUARES (BOARD_SIDE_LENGTH * BOARD_SIDE_LENGTH)

namespace operations_research {
namespace sat {

const Domain BOARD_SIDE_LENGTH_DOMAIN(0, BOARD_SIDE_LENGTH);
const Domain PAWN_COST_SUM_DOMAIN(0, MAX_PAWN_COST_SUM);
const Domain SINGLE_PAWN_DISPLACEMENT_DOMAIN(-SQUARES_ON_PAWN_FILE + 1,
                                             SQUARES_ON_PAWN_FILE -
                                                 1); // less starting file
const Domain ZERO_DOMAIN(0, 0);
const Domain NUM_CAPT_DOMAIN(0, NUM_CAPTURABLE_CHESSMEN_PSIDE);
const Domain PROMOTIONS_DOMAIN(-BOARD_SIDE_LENGTH, BOARD_SIDE_LENGTH);

extern "C" {
int get_row_num(int n) { return (n / BOARD_SIDE_LENGTH); }
int get_col_num(int n) { return (n % BOARD_SIDE_LENGTH); }
}

extern "C" bool FilterPawn(uint64_t _pawns[NUM_SIDES], uint64_t enpassant,
                           int capturedBasePieces[NUM_SIDES],
                           int minPromotionsToAccountFor[NUM_SIDES]) {
  CpModelBuilder cp_model;
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

  vector<BoolVar> fileMatchOptions[NUM_SIDES][BOARD_SIDE_LENGTH];
  vector<BoolVar> fileMatchOptionsNegated[NUM_SIDES][BOARD_SIDE_LENGTH];
  vector<BoolVar> fileCoverOptions[NUM_SIDES][BOARD_SIDE_LENGTH];
  vector<BoolVar> fileCoverOptionsNegated[NUM_SIDES][BOARD_SIDE_LENGTH];

  vector<IntVar> pawnVars[NUM_SIDES];
  vector<IntVar> pawnCosts[NUM_SIDES];
  IntVar pawnCostSum[NUM_SIDES];
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
          cp_model.AddNotEqual(pawnVars[i][p], k)
              .OnlyEnforceIf(pIsMatchedToK.Not());
          fileMatchOptions[i][k].push_back(pIsMatchedToK);
          fileMatchOptionsNegated[i][k].push_back(pIsMatchedToK.Not());

          BoolVar pCoversK = cp_model.NewBoolVar();
          cp_model.AddLessOrEqual(pawnVars[i][p], k).OnlyEnforceIf(pCoversK);
          cp_model.AddGreaterThan(pawnVars[i][p], k)
              .OnlyEnforceIf(pCoversK.Not());
          fileCoverOptions[i][k].push_back(pCoversK);
          fileCoverOptionsNegated[i][k].push_back(pCoversK.Not());
        }

        BoolVar pIsMatchedToFile = cp_model.NewBoolVar();
        cp_model.AddEquality(pawnVars[i][p], col)
            .OnlyEnforceIf(pIsMatchedToFile);
        cp_model.AddNotEqual(pawnVars[i][p], col)
            .OnlyEnforceIf(pIsMatchedToFile.Not());
        cp_model.AddEquality(pIsMatchedToFile, matchedToSameFileAsPawn[i][col]);
        fileMatchOptions[i][col].push_back(pIsMatchedToFile);
        fileMatchOptionsNegated[i][col].push_back(pIsMatchedToFile.Not());

        for (int k = col + 1; k <= end; k++) {
          BoolVar pIsMatchedToK = cp_model.NewBoolVar();
          cp_model.AddEquality(pawnVars[i][p], k).OnlyEnforceIf(pIsMatchedToK);
          cp_model.AddNotEqual(pawnVars[i][p], k)
              .OnlyEnforceIf(pIsMatchedToK.Not());
          fileMatchOptions[i][k].push_back(pIsMatchedToK);
          fileMatchOptionsNegated[i][k].push_back(pIsMatchedToK.Not());

          BoolVar pCoversK = cp_model.NewBoolVar();
          cp_model.AddGreaterOrEqual(pawnVars[i][p], k).OnlyEnforceIf(pCoversK);
          cp_model.AddLessThan(pawnVars[i][p], k).OnlyEnforceIf(pCoversK.Not());
          fileCoverOptions[i][k].push_back(pCoversK);
          fileCoverOptionsNegated[i][k].push_back(pCoversK.Not());
        }

        pawnCosts[i].push_back(
            cp_model.NewIntVar(SINGLE_PAWN_DISPLACEMENT_DOMAIN));
        cp_model.AddAbsEquality(pawnCosts[i][p], pawnVars[i][p] - p);
      }
    }
    cp_model.AddAllDifferent(pawnVars[i]);

    pawnCostSum[i] = cp_model.NewIntVar(PAWN_COST_SUM_DOMAIN);
    cp_model.AddEquality(pawnCostSum[i], LinearExpr::Sum(pawnCosts[i]));
  }

  for (int i = 0; i < NUM_SIDES; i++) {
    for (int j = 0; j < BOARD_SIDE_LENGTH; j++) {
      cp_model.AddBoolOr(fileMatchOptions[i][j])
          .OnlyEnforceIf(matchedStartingSquares[i][j]);
      cp_model.AddBoolAnd(fileMatchOptionsNegated[i][j])
          .OnlyEnforceIf(matchedStartingSquares[i][j].Not());

      cp_model.AddBoolOr(fileCoverOptions[i][j])
          .OnlyEnforceIf(startingSquareIsCovered[i][j]);
      cp_model.AddBoolAnd(fileCoverOptionsNegated[i][j])
          .OnlyEnforceIf(startingSquareIsCovered[i][j].Not());
    }
  }

  BoolVar pawnWhichStartedOn_isConsumedByPawnOnBoard[NUM_SIDES]
                                                    [BOARD_SIDE_LENGTH];
  IntVar numConsumed[NUM_SIDES];
  for (int i = 0; i < NUM_SIDES; i++) {
    int opp = (i + 1) % NUM_SIDES;
    for (int j = 0; j < BOARD_SIDE_LENGTH; j++) {
      BoolVar literals[] = {matchedStartingSquares[i][j].Not(),
                            startingSquareIsCovered[i][j]};
      pawnWhichStartedOn_isConsumedByPawnOnBoard[i][j] = cp_model.NewBoolVar();
      cp_model.AddBoolAnd(literals).OnlyEnforceIf(
          pawnWhichStartedOn_isConsumedByPawnOnBoard[i][j]);
    }
    numConsumed[i] = cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN);
    cp_model.AddEquality(
        numConsumed[i],
        LinearExpr::Sum(pawnWhichStartedOn_isConsumedByPawnOnBoard[i]));
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
      BoolVar literals[] = {
          matchedStartingSquares[0][j].Not(),
          matchedStartingSquares[1][j + i].Not(),
          pawnWhichStartedOn_isConsumedByPawnOnBoard[0][j].Not(),
          pawnWhichStartedOn_isConsumedByPawnOnBoard[1][j + i].Not(),
          pawnWhichStartedOn_isPaired[other_direction][0][j].Not(),
          pawnWhichStartedOn_isPaired[other_direction][1][j + i].Not()};
      cp_model.AddBoolAnd(literals).OnlyEnforceIf(
          pawnWhichStartedOn_isPaired[d][0][j]);
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

  IntVar numSameFilePromotions[NUM_SIDES];
  for (int i = 0; i < NUM_SIDES; i++) {
    int opp = (i + 1) % NUM_SIDES;
    for (int j = 0; j < BOARD_SIDE_LENGTH; j++) {
      BoolVar literals[] = {
          matchedStartingSquares[i][j].Not(),
          pawnWhichStartedOn_isConsumedByPawnOnBoard[i][j].Not(),
          pawnWhichStartedOn_isPaired[0][i][j].Not(),
          pawnWhichStartedOn_isPaired[1][i][j].Not(),
          matchedToSameFileAsPawn[opp][j].Not(),
          sameFilePromotions[opp][j].Not(),
      };
      cp_model.AddBoolAnd(literals).OnlyEnforceIf(sameFilePromotions[i][j]);
    }
    numSameFilePromotions[i] = cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN);
    cp_model.AddEquality(numSameFilePromotions[i],
                         LinearExpr::Sum(sameFilePromotions[i]));
  }

  BoolVar capturedPiecePromotions[NUM_SIDES][BOARD_SIDE_LENGTH];
  IntVar numCapturedPiecePromotions[NUM_SIDES];
  IntVar numCapturedBasePiecePromotions[NUM_SIDES];
  IntVar numCapturedPromotionPromotions[NUM_SIDES];
  BoolVar sideBorrowed[NUM_SIDES] = {cp_model.NewBoolVar(),
                                     cp_model.NewBoolVar()};
  for (int i = 0; i < NUM_SIDES; i++) {
    int opp = (i + 1) % NUM_SIDES;
    for (int j = 0; j < BOARD_SIDE_LENGTH; j++) {
      BoolVar literals[] = {
          matchedStartingSquares[i][j].Not(),
          pawnWhichStartedOn_isConsumedByPawnOnBoard[i][j].Not(),
          pawnWhichStartedOn_isPaired[0][i][j].Not(),
          pawnWhichStartedOn_isPaired[1][i][j].Not(),
          sameFilePromotions[i][j].Not(),
      };
      capturedPiecePromotions[i][j] = cp_model.NewBoolVar();
      cp_model.AddBoolAnd(literals).OnlyEnforceIf(
          capturedPiecePromotions[i][j]);
    }

    numCapturedPiecePromotions[i] =
        cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN);
    cp_model.AddEquality(numCapturedPiecePromotions[i],
                         LinearExpr::Sum(capturedPiecePromotions[i]));

    numCapturedBasePiecePromotions[i] =
        cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN);
    cp_model.AddLessOrEqual(numCapturedBasePiecePromotions[i],
                            capturedBasePieces[opp]);
    numCapturedPromotionPromotions[i] =
        cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN);
    cp_model.AddEquality(numCapturedPiecePromotions[i],
                         numCapturedBasePiecePromotions[i] +
                             numCapturedPromotionPromotions[i]);

    cp_model.AddGreaterThan(numCapturedPromotionPromotions[i], 0)
        .OnlyEnforceIf(sideBorrowed[i]);
    cp_model.AddEquality(numCapturedPromotionPromotions[i], 0)
        .OnlyEnforceIf(sideBorrowed[i].Not());

    cp_model.AddLessOrEqual(numCapturedPromotionPromotions[i],
                            numSameFilePromotions[opp] +
                                numCapturedBasePiecePromotions[opp]);
  }
  cp_model.AddNotEqual(sideBorrowed[0], sideBorrowed[1]);

  IntVar promotionSurplus[NUM_SIDES];
  for (int i = 0; i < NUM_SIDES; i++) {
    int opp = (i + 1) % NUM_SIDES;
    promotionSurplus[i] = cp_model.NewIntVar(PROMOTIONS_DOMAIN);
    cp_model.AddEquality(promotionSurplus[i],
                         pairsGoTo[i] + numCapturedPiecePromotions[i] +
                             numSameFilePromotions[i] -
                             numCapturedPromotionPromotions[opp] -
                             minPromotionsToAccountFor[i]);
    cp_model.AddGreaterOrEqual(promotionSurplus[i], 0);

    cp_model.AddLessOrEqual(numCapturedPromotionPromotions[i],
                            promotionSurplus[i]);
  }

  IntVar numPiecesPotentiallyCapturedForPawnDiagram[NUM_SIDES];
  for (int i = 0; i < NUM_SIDES; i++) {
    int opp = (i + 1) % NUM_SIDES;
    numPiecesPotentiallyCapturedForPawnDiagram[i] =
        cp_model.NewIntVar(NUM_CAPT_DOMAIN);
    cp_model.AddEquality(
        numPiecesPotentiallyCapturedForPawnDiagram[i],
        (capturedBasePieces[opp] - numCapturedBasePiecePromotions[opp]) +
            (promotionSurplus[i] - numCapturedPromotionPromotions[i]));
    cp_model.AddGreaterOrEqual(
        numConsumed[i] + numPiecesPotentiallyCapturedForPawnDiagram[i],
        pawnCostSum[i]);
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
