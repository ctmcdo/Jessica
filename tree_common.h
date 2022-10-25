#pragma once
#include <gmp.h>
#include <inttypes.h>

#include "chess_constants.h"

#define MAX_OF_ANY_BASE_PIECE 2
#define NUM_COVERED_SETS 14
#define MAX_UNIQUE_COSTS 4

#define FOUR_FACTORIAL (4 * 3 * 2)

#define ENPASSANT_EDGE_AND_RIGHT_VARIATIONS ((BOARD_SIDE_LENGTH - 1) + 1)
#define ENPASSANT_LEFT_LESS_EDGE_VARIATIONS (BOARD_SIDE_LENGTH - 2)

#define ONE_FIXED_ROOK_VARIATIONS 2

extern int
    fr_coveredSet_indices[NUM_FIXED_ROOK_SCENARIOS][MAX_OF_ANY_BASE_PIECE + 1]
                         [MAX_OF_ANY_BASE_PIECE + 1][MAX_OF_ANY_BASE_PIECE + 1]
                         [MAX_OF_ANY_BASE_PIECE + 1];
extern int fr_coveredSet_perms[NUM_FIXED_ROOK_SCENARIOS][NUM_COVERED_SETS]
                              [FOUR_FACTORIAL][NUM_PIECE_TYPES_LESS_KING];
extern int fr_coveredSet_perm_cost_boundaries[NUM_FIXED_ROOK_SCENARIOS]
                                             [NUM_COVERED_SETS]
                                             [MAX_UNIQUE_COSTS];

typedef struct position_node {
  mpz_t num_positions;
  struct position_node **children;
  char num_children;
} position_node;

enum castling_rights_cases {
  NO_CASTLING_RIGHTS,
  CASTLING_RIGHTS_ONE_SIDE,
  CASTLING_RIGHTS_BOTH_SIDES
};

enum enpassant_cases {
  NO_ENPASSANT,
  ENPASSANT_EDGE_AND_RIGHT,
  ENPASSANT_LEFT_LESS_EDGE
};
