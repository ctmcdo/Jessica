#include "chess.h"

// (row, column) to bitboard
uint64_t rcb(int row, int col) {
  return 1UL << ((row * BOARD_SIDE_LENGTH) + col);
}
