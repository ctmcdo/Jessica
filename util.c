#include <assert.h>
#include <inttypes.h>

#include "chess_constants.h"

extern uint64_t binomials[NUM_SQUARES + 1][MAX_BISHOPS_PSIDE + 1];

// (row, column) to bit (char)
char rcc(char row, char col) { return (row * BOARD_SIDE_LENGTH) + col; }

// (row, column) to bitboard
uint64_t rcb(char row, char col) { return 1UL << rcc(row, col); }

char get_row_num(char n) {
  assert((n >= 0) && (n < NUM_SQUARES));
  return (n / BOARD_SIDE_LENGTH);
}

char get_col_num(char n) {
  assert((n >= 0) && (n < NUM_SQUARES));
  return (n % BOARD_SIDE_LENGTH);
}

char get_index_of_1st_set_bit(uint64_t n) {
  for (int i = 0; i < NUM_SQUARES; i++) {
    if (n & 1)
      return i;
    n = n >> 1;
  }
  return -1;
}

// taken from chessprogrammingwiki->flipping,
// mirroring and rotating.
uint64_t rotate_bitboard_across_central_rows(uint64_t x) {
  return ((x << 56)) | ((x << 40) & (uint64_t)(0x00ff000000000000)) |
         ((x << 24) & (uint64_t)(0x0000ff0000000000)) |
         ((x << 8) & (uint64_t)(0x000000ff00000000)) |
         ((x >> 8) & (uint64_t)(0x00000000ff000000)) |
         ((x >> 24) & (uint64_t)(0x0000000000ff0000)) |
         ((x >> 40) & (uint64_t)(0x000000000000ff00)) | ((x >> 56));
}

void compute_binomials() {
  for (int n = 0; n <= NUM_SQUARES; n++) {
    binomials[n][0] = 1;
  }

  for (int k = 1; k <= MAX_BISHOPS_PSIDE; k++) {
    for (int n = k; n <= NUM_SQUARES; n++) {
      binomials[n][k] = binomials[n - 1][k - 1] + binomials[n - 1][k];
    }
  }
}
