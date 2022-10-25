#include <assert.h>
#include <immintrin.h>
#include <inttypes.h>
#include <nmmintrin.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "chess_constants.h"
#include "hungarian.h"
#include "position.h"
#include "util.h"

#define max(a, b) ((a) > (b) ? (a) : (b))

#define INFIN (uint64_t)100

enum PAWN_COST_ERROR_CODES { SIDE0_TOO_EXPENSIVE = 1, SIDE1_TOO_EXPENSIVE };

typedef struct CostTestStruct {
  bool enpassant;
  int left;
  int right;
  int tz_sum;
  int v;
  double **vertices;
} CostTestStruct;

typedef struct CostStruct {
  int side0;
  int side1;
} CostStruct;

const CostTestStruct DefaultCostTestStruct = {.enpassant = false,
                                              .left = 0,
                                              .right = BOARD_SIDE_LENGTH,
                                              .v = 0,
                                              .vertices = NULL,
                                              .tz_sum = 0};

void assert_mask(int n, uint64_t column) {
  assert((0 <= n) && (n < 56));
  if (column != 0) {
    assert(column < rcb(7, 0));
  }
  assert(((1UL << n) & column) == 0);
}

bool is_kth_bit_set(uint64_t n, int k) {
  assert((k >= 0) && (k < NUM_SQUARES));
  return (1UL << k) & n;
}

uint64_t get_mask_for_side0(int bp, uint64_t side0_column) {
  assert_mask(bp, side0_column);
  return (~__blsmsk_u64(1UL << bp)) & side0_column;
}

uint64_t get_mask_for_side1(int wp, uint64_t side1_column) {
  assert_mask(wp, side1_column);
  return __blsmsk_u64(1UL << wp) & side1_column;
}

void free_cost_test_struct(CostTestStruct *test_struct) {
  for (unsigned int i = 0; i < BOARD_SIDE_LENGTH; i++) {
    free((*test_struct).vertices[i]);
  }
  free((*test_struct).vertices);
  free(test_struct);
}

uint64_t get_cost(uint64_t pawns, uint64_t mask, uint64_t enpassant,
                  CostTestStruct *test_struct_) {
  // mallocing instead of creating local variables for test_struct
  // fields for readability.
  CostTestStruct *test_struct;
  if (test_struct_ == NULL) {
    test_struct = (CostTestStruct *)malloc(sizeof(CostTestStruct));
  } else {
    test_struct = test_struct_;
  }
  (*test_struct) = DefaultCostTestStruct;

  double **vertices;
  vertices = (double **)malloc(BOARD_SIDE_LENGTH * sizeof(double *));
  for (unsigned int i = 0; i < BOARD_SIDE_LENGTH; i++) {
    vertices[i] = calloc(BOARD_SIDE_LENGTH, sizeof(double));
  }

  uint64_t curr_position = pawns;
  int tz = __tzcnt_u64(curr_position);
  unsigned int v;
  for (v = 0; tz != NUM_SQUARES; v++) {
    (*test_struct).tz_sum += tz;

    for (unsigned int i = 0; i < BOARD_SIDE_LENGTH; i++) {
      vertices[v][i] = INFIN;
    }

    unsigned int col = get_col_num(tz);
    if (enpassant == (1UL << tz)) {
      assert(!is_kth_bit_set(mask, tz));

      vertices[v][col] = 0;
      (*test_struct).enpassant = true;
    } else {
      unsigned int row = get_row_num(tz);

      unsigned int left = min(col + row, 7);
      if (left > (*test_struct).left) {
        (*test_struct).left = left;
      }
      unsigned int right = max(col - row, 0);
      if (right < (*test_struct).right) {
        (*test_struct).right = right;
      }
      for (unsigned int i = right; i <= left; i++) {
        vertices[v][i] = abs(col - i);
      }
      if (is_kth_bit_set(mask, tz)) {
        if (tz >= (2 * BOARD_SIDE_LENGTH)) {
          vertices[v][col] = 2;
        } else {
          vertices[v][col] = INFIN;
        }
      }
    }

    curr_position = __blsr_u64(curr_position);
    tz = __tzcnt_u64(curr_position);
  }
  assert(v <= BOARD_SIDE_LENGTH);

  (*test_struct).v = v;
  (*test_struct).vertices = vertices;

  hungarian_problem_t h;
  hungarian_init(&h, vertices, BOARD_SIDE_LENGTH, BOARD_SIDE_LENGTH,
                 HUNGARIAN_MODE_MINIMIZE_COST);
  uint64_t cost = hungarian_solve(&h);

  hungarian_free(&h);

  if (test_struct_ == NULL) {
    free_cost_test_struct(test_struct);
  }

  return cost;
}

int get_column_mask(uint64_t side0, uint64_t side1,
                    uint64_t masks[BOARD_SIDE_LENGTH][BOARD_SIDE_LENGTH],
                    int col) {
  uint64_t side0_column = side0 & (COLUMN_1 << col);
  uint64_t side1_column = side1 & (COLUMN_1 << col);

  // := (63 - 64) == -1 if no side0 pawns on column
  int highest_side0 = (NUM_SQUARES - 1) - __builtin_clzll(side0_column);
  int lowest_side1 = __tzcnt_u64(side1_column);

  int num_masks = 0;
  if (highest_side0 > lowest_side1) {
    uint64_t wcc = side0_column;
    int hwc = highest_side0;
    while (hwc > lowest_side1) {
      masks[col][num_masks] = get_mask_for_side1(hwc, side1_column);
      wcc &= ~(1UL << hwc);
      hwc = (NUM_SQUARES - 1) - __builtin_clzll(wcc);
      ++num_masks;
    }

    while (lowest_side1 < highest_side0) {
      masks[col][num_masks] = get_mask_for_side0(lowest_side1, side0_column);
      side1_column &= ~(1UL << lowest_side1);
      lowest_side1 = __tzcnt_u64(side1_column);
      ++num_masks;
    }
  }

  return num_masks;
}

int num_mask_combinations(int num_masks[BOARD_SIDE_LENGTH]) {
  int combinations = num_masks[BOARD_SIDE_LENGTH - 1];
  if (combinations == 0) {
    combinations = 1;
  }
  for (int i = BOARD_SIDE_LENGTH - 2; i >= 0; i--) {
    assert(num_masks[i] >= 0);

    if (num_masks[i] > 0) {
      combinations *= num_masks[i];
      combinations += num_masks[i];
    } else {
      combinations++;
    }
  }
  return combinations;
}

int find_affordable_mask(uint64_t side0, uint64_t side1, uint64_t enpassant,
                         uint64_t masks[BOARD_SIDE_LENGTH][BOARD_SIDE_LENGTH],
                         int num_masks[BOARD_SIDE_LENGTH], CostStruct max_costs,
                         uint64_t mask, unsigned int col, int *num_calls) {
  assert(col <= BOARD_SIDE_LENGTH && col >= 0);

  if (col == BOARD_SIDE_LENGTH) {
    int side0_cost =
        get_cost(side0 >> BOARD_SIDE_LENGTH, mask >> BOARD_SIDE_LENGTH,
                 enpassant >> BOARD_SIDE_LENGTH, NULL);
    if (side0_cost > max_costs.side0) {
      return SIDE0_TOO_EXPENSIVE;
    }

    int side1_cost = get_cost(
        rotate_bitboard_across_central_rows(side1) >> BOARD_SIDE_LENGTH,
        rotate_bitboard_across_central_rows(mask) >> BOARD_SIDE_LENGTH,
        rotate_bitboard_across_central_rows(enpassant) >> BOARD_SIDE_LENGTH,
        NULL);
    if (side1_cost > max_costs.side1) {
      return SIDE1_TOO_EXPENSIVE;
    }

    return 0;
  }

  bool couldSatisfySide0 = false;
  if (num_masks[col] > 0) {
    for (unsigned int i = 0; i < num_masks[col]; i++) {
      ++(*num_calls);
      int code = find_affordable_mask(side0, side1, enpassant, masks, num_masks,
                                      max_costs, mask + masks[col][i], col + 1,
                                      num_calls);
      if (code == SIDE1_TOO_EXPENSIVE) {
        couldSatisfySide0 = true;
      } else if (code == 0) {
        return 0;
      }
    }
  } else {
    ++(*num_calls);
    int code = find_affordable_mask(side0, side1, enpassant, masks, num_masks,
                                    max_costs, mask, col + 1, num_calls);
    if (code == SIDE1_TOO_EXPENSIVE) {
      couldSatisfySide0 = true;
    } else if (code == 0) {
      return 0;
    }
  }

  if (col == 0) {
    assert(*num_calls == num_mask_combinations(num_masks));
  }

  if (couldSatisfySide0) {
    return SIDE1_TOO_EXPENSIVE;
  }

  return SIDE0_TOO_EXPENSIVE;
}

int validate_cost_helper(position p, CostStruct max_costs) {
  uint64_t masks[BOARD_SIDE_LENGTH][BOARD_SIDE_LENGTH] = {0};
  int num_masks[BOARD_SIDE_LENGTH] = {0};
  for (unsigned int i = 0; i < BOARD_SIDE_LENGTH; i++) {
    num_masks[i] =
        get_column_mask(p.sides[0].pawns, p.sides[1].pawns, masks, i);
  }

  int num_calls = 0;
  return find_affordable_mask(p.sides[0].pawns, p.sides[1].pawns, p.enpassant,
                              masks, num_masks, max_costs, 0, 0, &num_calls);
}

int validate_pawn_cost(position p, int max_cost0, int max_cost1) {
  int num_capturable_chessmen[NUM_SIDES] = {0};
  for (unsigned int i = 0; i < NUM_SIDES; i++) {
    num_capturable_chessmen[i] += _mm_popcnt_u64(p.sides[i].pawns);
    for (unsigned int j = 0; j < NUM_PIECE_TYPES_LESS_KING; j++) {
      num_capturable_chessmen[i] += _mm_popcnt_u64(p.sides[i].pieces[j]);
    }
  }

  CostStruct max_costs;
  max_costs.side0 = max_cost0;
  max_costs.side1 = max_cost1;

  return validate_cost_helper(p, max_costs);
}
