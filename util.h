#include <inttypes.h>

#include "chess_constants.h"

#define min(a, b) ((a) < (b) ? (a) : (b))
#define min3(a, b, c) min(min(a, b), c)

#define ROW_1 ((1UL << BOARD_SIDE_LENGTH) - 1)
#define COLUMN_1 0x0101010101010101

uint64_t binomials[NUM_SQUARES + 1][MAX_BISHOPS_PSIDE + 1];

// (row, column) to bit (char)
char rcc(char row, char col);

// (row, column) to bitboard
uint64_t rcb(char row, char col);

char get_row_num(char n);

char get_col_num(char n);

char get_index_of_1st_set_bit(uint64_t n);

// taken from chessprogrammingwiki->flipping,
// mirroring and rotating.
uint64_t rotate_bitboard_across_central_rows(uint64_t x);

void compute_binomials();
