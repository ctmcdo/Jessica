#include "chess.h"

#define PROCESSED_FEN_MAX_LENGTH (NUM_SQUARES + 1 + 4 + 2 + 3)

uint64_t rcb(int row, int col);

char get_row_num(char n);
char get_col_num(char n);

uint64_t rotate_bitboard_across_central_rows(uint64_t x);
void     rotate_position_across_central_rows(position* p);
