#include "chess.h"

#define ROOK0_BIT 0
#define ROOK1_BIT (BOARD_SIDE_LENGTH - 1)
#define ROOK2_BIT (BOARD_SIDE_LENGTH * (BOARD_SIDE_LENGTH - 1))
#define ROOK3_BIT (ROOK2_BIT + (BOARD_SIDE_LENGTH - 1))
#define a_ASCII_decimal 97

void fen_helper(uint64_t chessmen, char symbol, int squares[NUM_SQUARES]) {
  int tz = _tzcnt_u64(chessmen);
  while (tz != NUM_SQUARES) {
    squares[tz] = symbol;
    chessmen ^= 1UL << tz;
    tz = _tzcnt_u64(chessmen);
  }
}

char get_index_of_1st_set_bit(uint64_t n) {
  for (int i = 0; i < NUM_SQUARES; i++) {
    if (n & 1)
      return i;
    n = n >> 1;
  }
  return -1;
}

void print_fen(position p) {
  int squares[NUM_SQUARES];
  for (int i = 0; i < NUM_SQUARES; i++) {
    squares[i] = '1';
  }

  fen_helper(p.sides[p.side0isBlack].pawns, 'P', squares);
  fen_helper(p.sides[p.side0isBlack].pieces[BISHOP], 'B', squares);
  fen_helper(p.sides[p.side0isBlack].pieces[KNIGHT], 'N', squares);
  fen_helper(p.sides[p.side0isBlack].pieces[ROOK], 'R', squares);
  fen_helper(p.sides[p.side0isBlack].pieces[QUEEN], 'Q', squares);
  fen_helper(p.sides[p.side0isBlack].pieces[KING], 'K', squares);

  fen_helper(p.sides[!p.side0isBlack].pawns, 'p', squares);
  fen_helper(p.sides[!p.side0isBlack].pieces[BISHOP], 'b', squares);
  fen_helper(p.sides[!p.side0isBlack].pieces[KNIGHT], 'n', squares);
  fen_helper(p.sides[!p.side0isBlack].pieces[ROOK], 'r', squares);
  fen_helper(p.sides[!p.side0isBlack].pieces[QUEEN], 'q', squares);
  fen_helper(p.sides[!p.side0isBlack].pieces[KING], 'k', squares);

  for (int i = BOARD_SIDE_LENGTH - 1; i >= 0; i--) {
    for (int j = BOARD_SIDE_LENGTH - 1; j >= 0; j--) {
      printf("%c", squares[i * BOARD_SIDE_LENGTH + j]);
    }
    if (i > 0) {
      printf("/");
    }
  }
  printf(" ");

  char colour = 'w';
  if (p.side0isBlack) {
    colour = 'b';
  }
  printf("%c ", colour);

  bool cr0 = p.fixed_rooks & (1UL << ROOK0_BIT);
  bool cr1 = p.fixed_rooks & (1UL << ROOK1_BIT);
  bool cr2 = p.fixed_rooks & (1UL << ROOK2_BIT);
  bool cr3 = p.fixed_rooks & (1UL << ROOK3_BIT);
  if (!cr0 && !cr1 && !cr2 && !cr3) {
    printf("-");
  } else {
    if (cr0) {
      printf("K");
    }
    if (cr1) {
      printf("Q");
    }
    if (cr2) {
      printf("k");
    }
    if (cr3) {
      printf("q");
    }
  }
  printf(" ");

  if (p.enpassant) {
    char ep = get_index_of_1st_set_bit(p.enpassant);

    char column_letter =
        a_ASCII_decimal + (BOARD_SIDE_LENGTH - 1) - (ep % BOARD_SIDE_LENGTH);
    assert(a_ASCII_decimal <= column_letter && column_letter <= 104);
    printf("%c", column_letter);

    char row_num = 1 + ((int)(ep / BOARD_SIDE_LENGTH));
    assert(row_num == 4 || row_num == 5);
    if (row_num == 4) {
      // the enpassant square is the square behind the pawn that just moved 2
      printf("3");
    } else {
      printf("6");
    }

  } else {
    printf("-");
  }

  printf("\n");
}
