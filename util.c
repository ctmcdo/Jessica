#include "chess.h"

char get_row_num(char n) {
  assert((n >= 0) && (n < NUM_SQUARES));
  return (n / BOARD_SIDE_LENGTH);
}

char get_col_num(char n) {
  assert((n >= 0) && (n < NUM_SQUARES));
  return (n % BOARD_SIDE_LENGTH);
}
