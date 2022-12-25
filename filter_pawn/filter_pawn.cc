#include "filter_pawn.h"

#include <immintrin.h>
#include <inttypes.h>
#include <nmmintrin.h>
#include <stdint.h>
#include <stdlib.h>

#include "ortools/base/logging.h"
#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/util/sorted_interval_list.h"
#include <algorithm>
#include <vector>

using namespace std;

#define BOARD_SIDE_LENGTH                    8
#define MAX_COST_PSIDE                       25
#define NUM_CAPTURABLE_CHESSMEN_PSIDE        15
#define NUM_SQUARES                          (BOARD_SIDE_LENGTH * BOARD_SIDE_LENGTH)
#define NUM_PAWN_OCCUPIABLE_SQUARES_PER_FILE (BOARD_SIDE_LENGTH - 2)

#define FORWARD_SLASH 1
#define BACKSLASH     0

namespace operations_research
{
namespace sat
{

const Domain BOARD_SIDE_LENGTH_DOMAIN(0, BOARD_SIDE_LENGTH);
const Domain NUM_CAPTURABLE_CHESSMEN_PSIDE_DOMAIN(0, NUM_CAPTURABLE_CHESSMEN_PSIDE);
const Domain COST_PSIDE_DOMAIN(0, MAX_COST_PSIDE);
const Domain PAWN_RELATIVE_RANGE_DOMAIN(-NUM_PAWN_OCCUPIABLE_SQUARES_PER_FILE + 1,
                                        NUM_PAWN_OCCUPIABLE_SQUARES_PER_FILE
                                            - 1); // their magnitude less 1 to account
                                                  // for the file they're currently on
const Domain PROMOTIONS_DOMAIN(-BOARD_SIDE_LENGTH, BOARD_SIDE_LENGTH);
const Domain SINGLE_PAWN_COST_DOMAIN(0, NUM_PAWN_OCCUPIABLE_SQUARES_PER_FILE - 1);
const Domain ZERO_DOMAIN(0, 0);

extern "C"
{
    int get_row_num(int n)
    {
        return (n / BOARD_SIDE_LENGTH);
    }
    int get_col_num(int n)
    {
        return (n % BOARD_SIDE_LENGTH);
    }
}

extern "C" bool FilterPawn(uint64_t _pawns[NUM_SIDES],
                           uint64_t enpassant,
                           int      numCapturedBasePieces[NUM_SIDES],
                           int      minPromotions[NUM_SIDES])
{
    CpModelBuilder cp_model;
    BoolVar        matchedStartingSquares[NUM_SIDES][BOARD_SIDE_LENGTH];
    BoolVar        coveredStartingSquares[NUM_SIDES][BOARD_SIDE_LENGTH];
    BoolVar        fileSamePawnFileMatches[NUM_SIDES][BOARD_SIDE_LENGTH];
    BoolVar        zeroCostPromotions[NUM_SIDES][BOARD_SIDE_LENGTH];
    for (int i = 0; i < NUM_SIDES; i++)
    {
        for (int j = 0; j < BOARD_SIDE_LENGTH; j++)
        {
            matchedStartingSquares[i][j]  = cp_model.NewBoolVar();
            coveredStartingSquares[i][j]  = cp_model.NewBoolVar();
            fileSamePawnFileMatches[i][j] = cp_model.NewBoolVar();
            zeroCostPromotions[i][j]      = cp_model.NewBoolVar();
        }
    }

    vector<BoolVar> startingSquarePotentialMatches[NUM_SIDES][BOARD_SIDE_LENGTH];
    vector<BoolVar> startingSquarePotentialMatchesNegated[NUM_SIDES][BOARD_SIDE_LENGTH];
    vector<BoolVar> startingSquarePotentialSameFileMatches[NUM_SIDES][BOARD_SIDE_LENGTH];
    vector<BoolVar> startingSquarePotentialSameFileMatchesNegated[NUM_SIDES][BOARD_SIDE_LENGTH];
    vector<BoolVar> startingSquareCoverLiterals[NUM_SIDES][BOARD_SIDE_LENGTH];
    vector<BoolVar> startingSquareCoverLiteralsNegated[NUM_SIDES][BOARD_SIDE_LENGTH];

    vector<IntVar> pawnVars[NUM_SIDES];
    vector<IntVar> pawnCosts[NUM_SIDES];
    IntVar         pawnCostSum[NUM_SIDES];
    uint64_t       pawns[] = { _pawns[0], _pawns[1] };
    for (int i = 0; i < NUM_SIDES; i++)
    {
        if (pawns[i] == 0)
        {
            pawnCostSum[i] = cp_model.NewIntVar(ZERO_DOMAIN);
            continue;
        }

        int opp       = (i + 1) % NUM_SIDES;
        int num_pawns = _mm_popcnt_u64(pawns[i]);
        pawnVars[i].reserve(num_pawns);
        pawnCosts[i].reserve(num_pawns);

        for (int p = 0; p < num_pawns; p++)
        {
            int tz   = __tzcnt_u64(pawns[i]);
            pawns[i] = __blsr_u64(pawns[i]);

            int col = get_col_num(tz);
            if (enpassant == (1UL << tz))
            {
                const Domain single_col_domain(col, col);
                pawnVars[0].push_back(cp_model.NewIntVar(single_col_domain));
                startingSquarePotentialMatches[0][col].push_back(cp_model.TrueVar());
                cp_model.FixVariable(matchedStartingSquares[0][col], true);
                cp_model.FixVariable(fileSamePawnFileMatches[0][col], true);
                pawnCosts[0].push_back(cp_model.NewIntVar(ZERO_DOMAIN));
            }
            else
            {
                int row = get_row_num(tz);
                int start;
                int end;
                if (i == 0)
                {
                    start = std::max(col - (row - 1), 0);
                    end   = std::min(col + (row - 1), BOARD_SIDE_LENGTH - 1);
                }
                else
                {
                    start = std::max(col - ((BOARD_SIDE_LENGTH - 2) - row), 0);
                    end   = std::min(col + ((BOARD_SIDE_LENGTH - 2) - row), BOARD_SIDE_LENGTH - 1);
                }
                const Domain d(start, end);
                pawnVars[i].push_back(cp_model.NewIntVar(d));
                for (int k = start; k < col; k++)
                {
                    //  TODO: try make this into a function
                    BoolVar pIsMatchedToK = cp_model.NewBoolVar();
                    cp_model.AddEquality(pawnVars[i][p], k).OnlyEnforceIf(pIsMatchedToK);
                    cp_model.AddNotEqual(pawnVars[i][p], k).OnlyEnforceIf(pIsMatchedToK.Not());
                    startingSquarePotentialMatches[i][k].push_back(pIsMatchedToK);
                    startingSquarePotentialMatchesNegated[i][k].push_back(pIsMatchedToK.Not());

                    BoolVar pCoversK = cp_model.NewBoolVar();
                    cp_model.AddLessThan(pawnVars[i][p], k).OnlyEnforceIf(pCoversK);
                    cp_model.AddGreaterOrEqual(pawnVars[i][p], k).OnlyEnforceIf(pCoversK.Not());
                    startingSquareCoverLiterals[opp][k].push_back(pCoversK);
                    startingSquareCoverLiteralsNegated[opp][k].push_back(pCoversK.Not());
                }

                BoolVar pIsMatchedToFile = cp_model.NewBoolVar();
                cp_model.AddEquality(pawnVars[i][p], col).OnlyEnforceIf(pIsMatchedToFile);
                cp_model.AddNotEqual(pawnVars[i][p], col).OnlyEnforceIf(pIsMatchedToFile.Not());
                startingSquarePotentialMatches[i][col].push_back(pIsMatchedToFile);
                startingSquarePotentialMatchesNegated[i][col].push_back(pIsMatchedToFile.Not());
                startingSquarePotentialSameFileMatches[i][col].push_back(pIsMatchedToFile);
                startingSquarePotentialSameFileMatchesNegated[i][col].push_back(
                    pIsMatchedToFile.Not());

                startingSquareCoverLiterals[opp][col].push_back(pIsMatchedToFile.Not());
                startingSquareCoverLiteralsNegated[opp][col].push_back(pIsMatchedToFile);

                for (int k = col + 1; k <= end; k++)
                {
                    BoolVar pIsMatchedToK = cp_model.NewBoolVar();
                    cp_model.AddEquality(pawnVars[i][p], k).OnlyEnforceIf(pIsMatchedToK);
                    cp_model.AddNotEqual(pawnVars[i][p], k).OnlyEnforceIf(pIsMatchedToK.Not());
                    startingSquarePotentialMatches[i][k].push_back(pIsMatchedToK);
                    startingSquarePotentialMatchesNegated[i][k].push_back(pIsMatchedToK.Not());

                    BoolVar pCoversK = cp_model.NewBoolVar();
                    cp_model.AddGreaterThan(pawnVars[i][p], k).OnlyEnforceIf(pCoversK);
                    cp_model.AddLessOrEqual(pawnVars[i][p], k).OnlyEnforceIf(pCoversK.Not());
                    startingSquareCoverLiterals[opp][k].push_back(pCoversK);
                    startingSquareCoverLiteralsNegated[opp][k].push_back(pCoversK.Not());
                }

                pawnCosts[i].push_back(cp_model.NewIntVar(SINGLE_PAWN_COST_DOMAIN));
                cp_model.AddAbsEquality(pawnCosts[i][p], pawnVars[i][p] - col);
            }
        }
        cp_model.AddAllDifferent(pawnVars[i]);

        pawnCostSum[i] = cp_model.NewIntVar(COST_PSIDE_DOMAIN);
        cp_model.AddEquality(pawnCostSum[i], LinearExpr::Sum(pawnCosts[i]));
    }

    for (int i = 0; i < NUM_SIDES; i++)
    {
        for (int j = 0; j < BOARD_SIDE_LENGTH; j++)
        {
            // I filed a bug related to this check in ortools:
            // https://github.com/google/or-tools/issues/3591
            if (startingSquarePotentialMatches[i][j].size() > 0)
            {
                cp_model.AddBoolOr(startingSquarePotentialMatches[i][j])
                    .OnlyEnforceIf(matchedStartingSquares[i][j]);
                cp_model.AddBoolAnd(startingSquarePotentialMatchesNegated[i][j])
                    .OnlyEnforceIf(matchedStartingSquares[i][j].Not());
            }
            else
            {
                cp_model.FixVariable(matchedStartingSquares[i][j], false);
            }

            if (startingSquarePotentialSameFileMatches[i][j].size() > 0)
            {
                cp_model.AddBoolOr(startingSquarePotentialSameFileMatches[i][j])
                    .OnlyEnforceIf(fileSamePawnFileMatches[i][j]);
                cp_model.AddBoolAnd(startingSquarePotentialSameFileMatchesNegated[i][j])
                    .OnlyEnforceIf(fileSamePawnFileMatches[i][j].Not());
            }
            else
            {
                cp_model.FixVariable(fileSamePawnFileMatches[i][j], false);
            }

            if (startingSquareCoverLiterals[i][j].size() > 0)
            {
                cp_model.AddBoolOr(startingSquareCoverLiterals[i][j])
                    .OnlyEnforceIf(coveredStartingSquares[i][j]);
                cp_model.AddBoolAnd(startingSquareCoverLiteralsNegated[i][j])
                    .OnlyEnforceIf(coveredStartingSquares[i][j].Not());
            }
            else
            {
                cp_model.FixVariable(coveredStartingSquares[i][j], false);
            }
        }
    }

    BoolVar pawnConsumedPawnWhichStartedOn[NUM_SIDES][BOARD_SIDE_LENGTH];
    IntVar  numPawnsConsumedByPawnsOnBoard[] = { cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN),
                                                 cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN) };
    for (int i = 0; i < NUM_SIDES; i++)
    {
        for (int j = 0; j < BOARD_SIDE_LENGTH; j++)
        {
            BoolVar notMatchedButCovered[]
                = { matchedStartingSquares[i][j].Not(), coveredStartingSquares[i][j] };
            pawnConsumedPawnWhichStartedOn[i][j] = cp_model.NewBoolVar();
            cp_model.AddBoolAnd(notMatchedButCovered)
                .OnlyEnforceIf(pawnConsumedPawnWhichStartedOn[i][j]);
        }
        cp_model.AddEquality(numPawnsConsumedByPawnsOnBoard[i],
                             LinearExpr::Sum(pawnConsumedPawnWhichStartedOn[i]));
    }

    BoolVar pawnWhichStartedOn_isPaired[2][NUM_SIDES][BOARD_SIDE_LENGTH];
    for (int d = 0; d < 2; d++)
    { // d for direction which is either left or right
        for (int j = 0; j < NUM_SIDES; j++)
        {
            for (int k = 0; k < BOARD_SIDE_LENGTH; k++)
            {
                pawnWhichStartedOn_isPaired[d][j][k] = cp_model.NewBoolVar();
            }
        }
    }
    cp_model.FixVariable(pawnWhichStartedOn_isPaired[BACKSLASH][0][BOARD_SIDE_LENGTH - 1], false);
    cp_model.FixVariable(pawnWhichStartedOn_isPaired[BACKSLASH][1][0], false);

    cp_model.FixVariable(pawnWhichStartedOn_isPaired[FORWARD_SLASH][0][0], false);
    cp_model.FixVariable(pawnWhichStartedOn_isPaired[FORWARD_SLASH][1][BOARD_SIDE_LENGTH - 1],
                         false);

    int d               = 1; // FORWARD_SLASH
    int file_lim        = BOARD_SIDE_LENGTH;
    int other_direction = BACKSLASH;
    for (int i = -1; i <= 1; i += 2)
    {
        for (int j = d; j < file_lim; j++)
        {
            cp_model.AddEquality(pawnWhichStartedOn_isPaired[d][0][j],
                                 pawnWhichStartedOn_isPaired[d][1][j + i]);

            BoolVar notMatchedNorConsumedNorPairedInOtherDirection[]
                = { matchedStartingSquares[0][j].Not(),
                    matchedStartingSquares[1][j + i].Not(),
                    pawnConsumedPawnWhichStartedOn[0][j].Not(),
                    pawnConsumedPawnWhichStartedOn[1][j + i].Not(),
                    pawnWhichStartedOn_isPaired[other_direction][0][j].Not(),
                    pawnWhichStartedOn_isPaired[other_direction][1][j + i].Not() };
            cp_model.AddBoolAnd(notMatchedNorConsumedNorPairedInOtherDirection)
                .OnlyEnforceIf(pawnWhichStartedOn_isPaired[d][0][j]);
        }
        d               = 0; // BACKSLASH
        file_lim        = BOARD_SIDE_LENGTH - 1;
        other_direction = FORWARD_SLASH;
    }
    IntVar num_backslash_pairs = cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN);
    cp_model.AddEquality(num_backslash_pairs,
                         LinearExpr::Sum(pawnWhichStartedOn_isPaired[BACKSLASH][0]));
    IntVar num_forward_slash_pairs = cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN);
    cp_model.AddEquality(num_forward_slash_pairs,
                         LinearExpr::Sum(pawnWhichStartedOn_isPaired[FORWARD_SLASH][0]));
    IntVar num_pairs = cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN);
    cp_model.AddEquality(num_pairs, num_backslash_pairs + num_forward_slash_pairs);
    IntVar pairsGoTo[] = { cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN),
                           cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN) };
    cp_model.AddLessOrEqual(pairsGoTo[0], num_pairs);
    cp_model.AddEquality(pairsGoTo[1], num_pairs - pairsGoTo[0]);

    IntVar numZeroCostPromotions[] = { cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN),
                                       cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN) };
    for (int i = 0; i < NUM_SIDES; i++)
    {
        int opp = (i + 1) % NUM_SIDES;
        for (int j = 0; j < BOARD_SIDE_LENGTH; j++)
        {
            BoolVar zeroCostPromotionLiterals[] = {
                matchedStartingSquares[i][j].Not(),
                pawnConsumedPawnWhichStartedOn[i][j].Not(),
                pawnWhichStartedOn_isPaired[BACKSLASH][i][j].Not(),
                pawnWhichStartedOn_isPaired[FORWARD_SLASH][i][j].Not(),
                fileSamePawnFileMatches[opp][j].Not(),
                zeroCostPromotions[opp][j].Not(),
            };
            cp_model.AddBoolAnd(zeroCostPromotionLiterals).OnlyEnforceIf(zeroCostPromotions[i][j]);
        }
        cp_model.AddEquality(numZeroCostPromotions[i], LinearExpr::Sum(zeroCostPromotions[i]));
    }

    BoolVar promsViaCapturedPiece[NUM_SIDES][BOARD_SIDE_LENGTH];
    IntVar  numPromsViaCapturedPiece[]     = { cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN),
                                               cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN) };
    IntVar  numPromsViaCapturedBasePiece[] = { cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN),
                                               cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN) };
    IntVar  numPromsViaCapturedPromotion[] = { cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN),
                                               cp_model.NewIntVar(BOARD_SIDE_LENGTH_DOMAIN) };
    BoolVar atLeastOnePromCapturedByPawn[] = { cp_model.NewBoolVar(), cp_model.NewBoolVar() };
    for (int i = 0; i < NUM_SIDES; i++)
    {
        int opp = (i + 1) % NUM_SIDES;
        for (int j = 0; j < BOARD_SIDE_LENGTH; j++)
        {
            BoolVar promViaCapturedPieceLiterals[] = {
                matchedStartingSquares[i][j].Not(),
                pawnConsumedPawnWhichStartedOn[i][j].Not(),
                pawnWhichStartedOn_isPaired[BACKSLASH][i][j].Not(),
                pawnWhichStartedOn_isPaired[FORWARD_SLASH][i][j].Not(),
                zeroCostPromotions[i][j].Not(),
            };
            promsViaCapturedPiece[i][j] = cp_model.NewBoolVar();
            cp_model.AddBoolAnd(promViaCapturedPieceLiterals)
                .OnlyEnforceIf(promsViaCapturedPiece[i][j]);
        }
        cp_model.AddEquality(numPromsViaCapturedPiece[i],
                             LinearExpr::Sum(promsViaCapturedPiece[i]));

        cp_model.AddLessOrEqual(numPromsViaCapturedBasePiece[i], numCapturedBasePieces[opp]);
        cp_model.AddLessOrEqual(numPromsViaCapturedPromotion[i],
                                numZeroCostPromotions[opp] + numPromsViaCapturedBasePiece[opp]);
        cp_model.AddEquality(numPromsViaCapturedPiece[i],
                             numPromsViaCapturedBasePiece[i] + numPromsViaCapturedPromotion[i]);

        cp_model.AddGreaterThan(numPromsViaCapturedPromotion[i], 0)
            .OnlyEnforceIf(atLeastOnePromCapturedByPawn[i]);
        cp_model.AddEquality(numPromsViaCapturedPromotion[i], 0)
            .OnlyEnforceIf(atLeastOnePromCapturedByPawn[i].Not());
    }
    cp_model.AddAtMostOne(atLeastOnePromCapturedByPawn);

    IntVar promsLeftForOppDiagram[NUM_SIDES];
    for (int i = 0; i < NUM_SIDES; i++)
    {
        int opp                   = (i + 1) % NUM_SIDES;
        promsLeftForOppDiagram[i] = cp_model.NewIntVar(PROMOTIONS_DOMAIN);
        cp_model.AddEquality(promsLeftForOppDiagram[i],
                             pairsGoTo[i] + numZeroCostPromotions[i] + numPromsViaCapturedPiece[i]
                                 - numPromsViaCapturedPromotion[opp] - minPromotions[i]);
        cp_model.Maximize(promsLeftForOppDiagram[i]);
        cp_model.AddGreaterOrEqual(promsLeftForOppDiagram[i], 0);
    }

    for (int i = 0; i < NUM_SIDES; i++)
    {
        int    opp = (i + 1) % NUM_SIDES;
        IntVar pawnDiagramPieceCaptureBudget
            = cp_model.NewIntVar(NUM_CAPTURABLE_CHESSMEN_PSIDE_DOMAIN);
        cp_model.AddEquality(pawnDiagramPieceCaptureBudget,
                             (numCapturedBasePieces[opp] - numPromsViaCapturedBasePiece[i])
                                 + promsLeftForOppDiagram[opp]);

        cp_model.AddGreaterOrEqual(
            numPawnsConsumedByPawnsOnBoard[opp] + pawnDiagramPieceCaptureBudget, pawnCostSum[i]);
    }

    Model model;
    model.Add(NewFeasibleSolutionObserver([&](const CpSolverResponse& r) {
        /*
    cout << "pawnVars\n";
    for (int i = 1; i >= 0; i--) {
      if (pawnVars[i].size() == 0) {
        printf("No pawns for side %d\n", i);
        continue;
      }
      for (int j = 0; j < pawnVars[i].size(); j++) {
        cout << SolutionIntegerValue(r, pawnVars[i][j]) << " ";
      }
      cout << "\n";
    }
    cout << "\n";
    cout << "pawnCosts\n";
    for (int i = 1; i >= 0; i--) {
      if (pawnVars[i].size() == 0) {
        printf("No pawns for side %d\n", i);
        continue;
      }
      for (int j = 0; j < pawnCosts[i].size(); j++) {
        cout << SolutionIntegerValue(r, pawnCosts[i][j]) << " ";
      }
      cout << "\n";
    }

    /*
    cout << "pawnConsumedPawnWhichStartedOn:\n";
    for (int i = 1; i >= 0; i--) {
    for (int j = BOARD_SIDE_LENGTH - 1; j >= 0; j--) {
    cout << SolutionIntegerValue(r, pawnConsumedPawnWhichStartedOn[i][j])
         << " ";
    }
    cout << "\n";
    }
    cout << "\n";

    cout << "pawnWhichStartedOn_isPaired:\n";
    for (int d = 0; d < 2; d++) {
    printf("direction %d:\n", d);
    for (int i = 1; i >= 0; i--) {
    for (int j = BOARD_SIDE_LENGTH - 1; j >= 0; j--) {
      cout << SolutionIntegerValue(r, pawnWhichStartedOn_isPaired[d][i][j])
           << " ";
    }
    cout << "\n";
    }
    }
    cout << "\n";

    cout << "zeroCostPromotions:\n";
    for (int i = 1; i >= 0; i--) {
    for (int j = BOARD_SIDE_LENGTH - 1; j >= 0; j--) {
    cout << SolutionIntegerValue(r, zeroCostPromotions[i][j]) << " ";
    }
    cout << "\n";
    }
    cout << "\n";

    /*
    cout << "promsViaCapturedPiece:\n";
    for (int i = 1; i >= 0; i--) {
    for (int j = BOARD_SIDE_LENGTH - 1; j >= 0; j--) {
    cout << SolutionIntegerValue(r, promsViaCapturedPiece[i][j]) << " ";
    }
    cout << "\n";
    }
    cout << "\n";

    cout << "numPromsViaCapturedBasePiece:\n";
    for (int i = 1; i >= 0; i--) {
    cout << SolutionIntegerValue(r, numPromsViaCapturedBasePiece[i]) << endl;
    }
    cout << "\n";

    cout << "numPromsViaCapturedPromotion:\n";
    for (int i = 1; i >= 0; i--) {
    cout << SolutionIntegerValue(r, numPromsViaCapturedPromotion[i]) << endl;
    }
    cout << "\n";

    cout << "promsLeftForOppDiagram:\n";
    for (int i = 1; i >= 0; i--) {
    cout << SolutionIntegerValue(r, promsLeftForOppDiagram[i]) << endl;
    }
    cout << "\n";
    /*
    for (int i = 1; i >= 0; i--) {
      int opp = (i + 1) % NUM_SIDES;
      printf("minPromotions[%d]: ", i);
      cout << SolutionIntegerValue(r, minPromotions[i]) << endl;
      printf("pawnCostSum[%d]: ", i);
      cout << SolutionIntegerValue(r, pawnCostSum[i]) << endl;
      cout << endl;

      printf("promsLeftForOppDiagram[%d]: ", i);
      cout << SolutionIntegerValue(r, promsLeftForOppDiagram[i]) << endl;
      printf("pairsGoTo[%d]: ", i);
      cout << SolutionIntegerValue(r, pairsGoTo[i]) << endl;
      printf("numZeroCostPromotions[%d]: ", i);
      cout << SolutionIntegerValue(r, numZeroCostPromotions[i]) << endl;
      printf("numPromsViaCapturedPiece[%d]: ", i);
      cout << SolutionIntegerValue(r, numPromsViaCapturedPiece[i]) << endl;
      cout << endl;

      printf("numCapturedBasePieces[%d]: ", opp);
      cout << SolutionIntegerValue(r, numCapturedBasePieces[opp]) << endl;
      printf("numPromsViaCapturedBasePiece[%d]:, ", i);
      cout << SolutionIntegerValue(r, numPromsViaCapturedBasePiece[i]) <<
    endl; printf("promsLeftForOppDiagram[%d]: ", opp); cout <<
    SolutionIntegerValue(r, promsLeftForOppDiagram[opp]) << endl;
      printf("numPromsViaCapturedPromotion[%d]: ", i);
      cout << SolutionIntegerValue(r, numPromsViaCapturedPromotion[i]) <<
    endl; cout << endl;
    }

    cout << "numPromsViaCapturedBasePiece:\n";
    for (int i = 1; i >= 0; i--) {
      cout << SolutionIntegerValue(r, numPromsViaCapturedBasePiece[i]) <<
    endl;
    }
    cout << "promsLeftForOppDiagram:\n";
    for (int i = 1; i >= 0; i--) {
      cout << SolutionIntegerValue(r, promsLeftForOppDiagram[i]) << endl;
    }
    cout << "numPawnsConsumedByPawnsOnBoard:\n";
    for (int i = 1; i >= 0; i--) {
      cout << SolutionIntegerValue(r, numPawnsConsumedByPawnsOnBoard[i])
           << endl;
    }
    cout << "pawnCostSum:\n";
    for (int i = 1; i >= 0; i--) {
      cout << SolutionIntegerValue(r, pawnCostSum[i]) << endl;

    }
    */
    }));

    const CpSolverResponse response = SolveCpModel(cp_model.Build(), &model);
    if (response.status() == CpSolverStatus::OPTIMAL
        || response.status() == CpSolverStatus::FEASIBLE)
    {
        return true;
    }
    else if (response.status() == CpSolverStatus::MODEL_INVALID)
    {
        cout << "Unexpected invalid model\n";
        exit(1);
    }
    return false;
} // namespace sat

} // namespace sat
} // namespace operations_research
