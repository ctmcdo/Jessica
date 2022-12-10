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
#define MAX_COST_PSIDE 25
#define NUM_CAPTURABLE_CHESSMEN_PSIDE 15
#define NUM_SQUARES (BOARD_SIDE_LENGTH * BOARD_SIDE_LENGTH)
#define NUM_PAWN_OCCUPIABLE_SQUARES_PER_FILE (BOARD_SIDE_LENGTH - 2)

#define LEFT 1
#define RIGHT 0

namespace operations_research {
namespace sat {

const Domain BOARD_SIDE_LENGTH_DOMAIN(0, BOARD_SIDE_LENGTH);
const Domain
    NUM_CAPTURABLE_CHESSMEN_PSIDE_DOMAIN(0, NUM_CAPTURABLE_CHESSMEN_PSIDE);
const Domain COST_PSIDE_DOMAIN(0, MAX_COST_PSIDE);
const Domain
    PAWN_RELATIVE_RANGE_DOMAIN(-NUM_PAWN_OCCUPIABLE_SQUARES_PER_FILE + 1,
                               NUM_PAWN_OCCUPIABLE_SQUARES_PER_FILE -
                                   1); // their magnitude less 1 to account
                                       // for the file they're currently on
const Domain PROMOTIONS_DOMAIN(-BOARD_SIDE_LENGTH, BOARD_SIDE_LENGTH);
const Domain SINGLE_PAWN_COST_DOMAIN(0,
                                     NUM_PAWN_OCCUPIABLE_SQUARES_PER_FILE - 1);
const Domain ZERO_DOMAIN(0, 0);

extern "C" {
int get_row_num(int n) { return (n / BOARD_SIDE_LENGTH); }
int get_col_num(int n) { return (n % BOARD_SIDE_LENGTH); }
}

extern "C" bool FilterPawn(uint64_t _pawns[NUM_SIDES], uint64_t enpassant,
                           int numCapturedBasePieces[NUM_SIDES],
                           int minPromotions[NUM_SIDES]) {
  CpModelBuilder cp_model;
  BoolVar matchedStartingSquares[NUM_SIDES][BOARD_SIDE_LENGTH];
  BoolVar coveredStartingSquares[NUM_SIDES][BOARD_SIDE_LENGTH];
  BoolVar fileSamePawnFileMatches[NUM_SIDES][BOARD_SIDE_LENGTH];
  BoolVar zeroCostPromotions[NUM_SIDES][BOARD_SIDE_LENGTH];
  for (int i = 0; i < NUM_SIDES; i++) {
    for (int j = 0; j < BOARD_SIDE_LENGTH; j++) {
      matchedStartingSquares[i][j] = cp_model.NewBoolVar();
      coveredStartingSquares[i][j] = cp_model.NewBoolVar();
      fileSamePawnFileMatches[i][j] = cp_model.NewBoolVar();
      zeroCostPromotions[i][j] = cp_model.NewBoolVar();
    }
  }

  vector<BoolVar> pawnVarMatchesFile[NUM_SIDES][BOARD_SIDE_LENGTH];
  vector<BoolVar> pawnVarDoesNotMatchFile[NUM_SIDES][BOARD_SIDE_LENGTH];
  vector<BoolVar> pawnVarCoversFile[NUM_SIDES][BOARD_SIDE_LENGTH];
  vector<BoolVar> pawnVarDoesNotCoverFile[NUM_SIDES][BOARD_SIDE_LENGTH];

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
        pawnVarMatchesFile[i][col].push_back(cp_model.TrueVar());
        cp_model.FixVariable(matchedStartingSquares[i][col], true);
        cp_model.FixVariable(fileSamePawnFileMatches[i][col], true);
        pawnCosts[i].push_back(cp_model.NewIntVar(ZERO_DOMAIN));

      } else {
        int row = get_row_num(tz);
        int start = std::max(col - row, 0);
        int end = std::min(col + row, BOARD_SIDE_LENGTH - 1);
        const Domain d(start, end);
        pawnVars[i].push_back(cp_model.NewIntVar(d));
        for (int k = start; k < col; k++) {
          // TODO: try make this into a function
          BoolVar pIsMatchedToK = cp_model.NewBoolVar();
          cp_model.AddEquality(pawnVars[i][p], k).OnlyEnforceIf(pIsMatchedToK);
          cp_model.AddNotEqual(pawnVars[i][p], k)
              .OnlyEnforceIf(pIsMatchedToK.Not());
          pawnVarMatchesFile[i][k].push_back(pIsMatchedToK);
          pawnVarDoesNotMatchFile[i][k].push_back(pIsMatchedToK.Not());

          BoolVar pCoversK = cp_model.NewBoolVar();
          cp_model.AddLessOrEqual(pawnVars[i][p], k).OnlyEnforceIf(pCoversK);
          cp_model.AddGreaterThan(pawnVars[i][p], k)
              .OnlyEnforceIf(pCoversK.Not());
          pawnVarCoversFile[i][k].push_back(pCoversK);
          pawnVarDoesNotCoverFile[i][k].push_back(pCoversK.Not());
        }

        BoolVar pIsMatchedToFile = cp_model.NewBoolVar();
        cp_model.AddEquality(pawnVars[i][p], col)
            .OnlyEnforceIf(pIsMatchedToFile);
        cp_model.AddNotEqual(pawnVars[i][p], col)
            .OnlyEnforceIf(pIsMatchedToFile.Not());
        cp_model.AddEquality(pIsMatchedToFile, fileSamePawnFileMatches[i][col]);
        pawnVarMatchesFile[i][col].push_back(pIsMatchedToFile);
        pawnVarDoesNotMatchFile[i][col].push_back(pIsMatchedToFile.Not());

        for (int k = col + 1; k <= end; k++) {
          BoolVar pIsMatchedToK = cp_model.NewBoolVar();
          cp_model.AddEquality(pawnVars[i][p], k).OnlyEnforceIf(pIsMatchedToK);
          cp_model.AddNotEqual(pawnVars[i][p], k)
              .OnlyEnforceIf(pIsMatchedToK.Not());
          pawnVarMatchesFile[i][k].push_back(pIsMatchedToK);
          pawnVarDoesNotMatchFile[i][k].push_back(pIsMatchedToK.Not());

          BoolVar pCoversK = cp_model.NewBoolVar();
          cp_model.AddGreaterOrEqual(pawnVars[i][p], k).OnlyEnforceIf(pCoversK);
          cp_model.AddLessThan(pawnVars[i][p], k).OnlyEnforceIf(pCoversK.Not());
          pawnVarCoversFile[i][k].push_back(pCoversK);
          pawnVarDoesNotCoverFile[i][k].push_back(pCoversK.Not());
        }

        pawnCosts[i].push_back(cp_model.NewIntVar(SINGLE_PAWN_COST_DOMAIN));
        IntVar offset = cp_model.NewIntVar(PAWN_RELATIVE_RANGE_DOMAIN);
        cp_model.AddEquality(offset, pawnVars[i][p] - col);
        cp_model.AddAbsEquality(pawnCosts[i][p], offset);
      }
    }
    cp_model.AddAllDifferent(pawnVars[i]);

    pawnCostSum[i] = cp_model.NewIntVar(COST_PSIDE_DOMAIN);
    cp_model.AddEquality(pawnCostSum[i], LinearExpr::Sum(pawnCosts[i]));
  }

  for (int i = 0; i < NUM_SIDES; i++) {
    for (int j = 0; j < BOARD_SIDE_LENGTH; j++) {
      // I filed a bug related to this check in ortools:
      // https://github.com/google/or-tools/issues/3591
      if (pawnVarMatchesFile[i][j].size() > 0) {
        cp_model.AddBoolOr(pawnVarMatchesFile[i][j])
            .OnlyEnforceIf(matchedStartingSquares[i][j]);
        cp_model.AddBoolAnd(pawnVarDoesNotMatchFile[i][j])
            .OnlyEnforceIf(matchedStartingSquares[i][j].Not());
      } else {
        cp_model.FixVariable(matchedStartingSquares[i][j], false);
      }

      if (pawnVarCoversFile[i][j].size()) {
        cp_model.AddBoolOr(pawnVarCoversFile[i][j])
            .OnlyEnforceIf(coveredStartingSquares[i][j]);
        cp_model.AddBoolAnd(pawnVarDoesNotCoverFile[i][j])
            .OnlyEnforceIf(coveredStartingSquares[i][j].Not());
      } else {
        cp_model.FixVariable(coveredStartingSquares[i][j], false);
      }
    }
  }

  BoolVar pawnConsumedPawnWhichStartedOn[NUM_SIDES][BOARD_SIDE_LENGTH];
  IntVar numPawnsConsumedByPawns[] = {
      cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN),
      cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN)};
  for (int i = 0; i < NUM_SIDES; i++) {
    int opp = (i + 1) % NUM_SIDES;
    for (int j = 0; j < BOARD_SIDE_LENGTH; j++) {
      BoolVar notMatchedButCovered[] = {matchedStartingSquares[i][j].Not(),
                                        coveredStartingSquares[i][j]};
      pawnConsumedPawnWhichStartedOn[i][j] = cp_model.NewBoolVar();
      cp_model.AddBoolAnd(notMatchedButCovered)
          .OnlyEnforceIf(pawnConsumedPawnWhichStartedOn[i][j]);
    }
    cp_model.AddEquality(numPawnsConsumedByPawns[i],
                         LinearExpr::Sum(pawnConsumedPawnWhichStartedOn[i]));
  }

  BoolVar pawnWhichStartedOn_isPaired[2][NUM_SIDES][BOARD_SIDE_LENGTH];
  for (int d = 0; d < 2; d++) { // d for direction which is either left or right
    for (int j = 0; j < NUM_SIDES; j++) {
      for (int k = 0; k < BOARD_SIDE_LENGTH; k++) {
        pawnWhichStartedOn_isPaired[d][j][k] = cp_model.NewBoolVar();
      }
    }
  }
  cp_model.FixVariable(pawnWhichStartedOn_isPaired[LEFT][0][0], false);
  cp_model.FixVariable(
      pawnWhichStartedOn_isPaired[LEFT][1][BOARD_SIDE_LENGTH - 1], false);

  cp_model.FixVariable(
      pawnWhichStartedOn_isPaired[RIGHT][0][BOARD_SIDE_LENGTH - 1], false);
  cp_model.FixVariable(pawnWhichStartedOn_isPaired[RIGHT][1][0], false);

  // TODO: set j to 1 and d and see different positions produced
  int d = 1;
  int file_lim = BOARD_SIDE_LENGTH;
  int other_direction = 0;
  for (int i = -1; i <= 1; i += 2) {
    for (int j = d; j < file_lim; j++) {
      cp_model.AddEquality(pawnWhichStartedOn_isPaired[d][0][j],
                           pawnWhichStartedOn_isPaired[d][1][j + i]);

      BoolVar notMatchedNorConsumedNorPairedInOtherDirection[] = {
          matchedStartingSquares[0][j].Not(),
          matchedStartingSquares[1][j + i].Not(),
          pawnConsumedPawnWhichStartedOn[0][j].Not(),
          pawnConsumedPawnWhichStartedOn[1][j + i].Not(),
          pawnWhichStartedOn_isPaired[other_direction][0][j].Not(),
          pawnWhichStartedOn_isPaired[other_direction][1][j + i].Not()};
      cp_model.AddBoolAnd(notMatchedNorConsumedNorPairedInOtherDirection)
          .OnlyEnforceIf(pawnWhichStartedOn_isPaired[d][0][j]);
    }
    d = 0;
    file_lim = BOARD_SIDE_LENGTH - 1;
    other_direction = 1;
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

  IntVar numSameFilePromotions[] = {
      cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN),
      cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN)};
  for (int i = 0; i < NUM_SIDES; i++) {
    int opp = (i + 1) % NUM_SIDES;
    for (int j = 0; j < BOARD_SIDE_LENGTH; j++) {
      BoolVar literals[] = {
          matchedStartingSquares[i][j].Not(),
          pawnConsumedPawnWhichStartedOn[i][j].Not(),
          pawnWhichStartedOn_isPaired[0][i][j].Not(),
          pawnWhichStartedOn_isPaired[1][i][j].Not(),
          fileSamePawnFileMatches[opp][j].Not(),
          zeroCostPromotions[opp][j].Not(),
      };
      cp_model.AddBoolAnd(literals).OnlyEnforceIf(zeroCostPromotions[i][j]);
    }
    cp_model.AddEquality(numSameFilePromotions[i],
                         LinearExpr::Sum(zeroCostPromotions[i]));
  }

  BoolVar capturedPiecePromotions[NUM_SIDES][BOARD_SIDE_LENGTH];
  IntVar numCapturedPiecePromotions[] = {
      cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN),
      cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN)};
  IntVar numCapturedBasePiecePromotions[] = {
      cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN),
      cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN)};
  IntVar numCapturedPromotionPromotions[] = {
      cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN),
      cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN)};
  BoolVar sideBorrowed[] = {cp_model.NewBoolVar(), cp_model.NewBoolVar()};
  for (int i = 0; i < NUM_SIDES; i++) {
    int opp = (i + 1) % NUM_SIDES;
    for (int j = 0; j < BOARD_SIDE_LENGTH; j++) {
      BoolVar literals[] = {
          matchedStartingSquares[i][j].Not(),
          pawnConsumedPawnWhichStartedOn[i][j].Not(),
          pawnWhichStartedOn_isPaired[0][i][j].Not(),
          pawnWhichStartedOn_isPaired[1][i][j].Not(),
          zeroCostPromotions[i][j].Not(),
      };
      capturedPiecePromotions[i][j] = cp_model.NewBoolVar();
      cp_model.AddBoolAnd(literals).OnlyEnforceIf(
          capturedPiecePromotions[i][j]);
    }
    cp_model.AddEquality(numCapturedPiecePromotions[i],
                         LinearExpr::Sum(capturedPiecePromotions[i]));
    cp_model.AddEquality(numCapturedPiecePromotions[i],
                         numCapturedBasePiecePromotions[i] +
                             numCapturedPromotionPromotions[i]);

    cp_model.AddLessOrEqual(numCapturedBasePiecePromotions[i],
                            numCapturedBasePieces[opp]);

    cp_model.AddLessOrEqual(numCapturedPromotionPromotions[i],
                            numSameFilePromotions[opp] +
                                numCapturedBasePiecePromotions[opp]);
    cp_model.AddGreaterThan(numCapturedPromotionPromotions[i], 0)
        .OnlyEnforceIf(sideBorrowed[i]);
    cp_model.AddEquality(numCapturedPromotionPromotions[i], 0)
        .OnlyEnforceIf(sideBorrowed[i].Not());
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
                             minPromotions[i]);
    cp_model.AddGreaterOrEqual(promotionSurplus[i], 0);
  }

  IntVar numPiecesPotentiallyCapturedForPawnDiagram[NUM_SIDES];
  for (int i = 0; i < NUM_SIDES; i++) {
    int opp = (i + 1) % NUM_SIDES;
    numPiecesPotentiallyCapturedForPawnDiagram[i] =
        cp_model.NewIntVar(NUM_CAPTURABLE_CHESSMEN_PSIDE_DOMAIN);
    cp_model.AddEquality(
        numPiecesPotentiallyCapturedForPawnDiagram[i],
        (numCapturedBasePieces[opp] - numCapturedBasePiecePromotions[opp]) +
            (promotionSurplus[i] - numCapturedPromotionPromotions[i]));
    cp_model.AddGreaterOrEqual(
        numPawnsConsumedByPawns[i] +
            numPiecesPotentiallyCapturedForPawnDiagram[i],
        pawnCostSum[i]);
  }

  Model model;
  model.Add(NewFeasibleSolutionObserver([&](const CpSolverResponse &r) {
    /*
    cout << "matchedStartingSquares:\n";
    for (int i = 0; i < NUM_SIDES; i++) {
      for (int j = BOARD_SIDE_LENGTH - 1; j >= 0; j--) {
        cout << SolutionIntegerValue(r, matchedStartingSquares[i][j]) << " ";
      }
      cout << "\n";
    }

    cout << "zeroCostPromotions:\n";
    for (int i = 0; i < NUM_SIDES; i++) {
      for (int j = BOARD_SIDE_LENGTH - 1; j >= 0; j--) {
        cout << SolutionIntegerValue(r, zeroCostPromotions[i][j]) << " ";
      }
      cout << "\n";
    }

    cout << "pawnVars\n";
    for (int i = 0; i < NUM_SIDES; i++) {
      for (int j = 0; j < pawnVars[i].size(); j++) {
        cout << SolutionIntegerValue(r, pawnVars[i][j]) << " ";
      }
      cout << "\n";
    }

    cout << "pawnCosts\n";
    for (int i = 0; i < NUM_SIDES; i++) {
      for (int j = 0; j < pawnCosts[i].size(); j++) {
        cout << SolutionIntegerValue(r, pawnCosts[i][j]) << " ";
      }
      cout << "\n";
    }
    */
  }));

  const CpSolverResponse response = SolveCpModel(cp_model.Build(), &model);
  if (response.status() == CpSolverStatus::OPTIMAL ||
      response.status() == CpSolverStatus::FEASIBLE) {
    return true;
  }
  return false;
} // namespace sat

} // namespace sat
} // namespace operations_research
