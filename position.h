#pragma once
#include <inttypes.h>
#include <stdbool.h>

#include "chess_constants.h"

typedef struct piece {
  int type;
  uint64_t bitboard;
} piece;

typedef struct side {
  uint64_t pieces[NUM_PIECE_TYPES];
  uint64_t pawns;
  uint64_t _fr;
} side;

typedef struct position {
  side sides[NUM_SIDES];
  uint64_t enpassant;
  uint64_t fixed_rooks;
  bool side0isBlack;
  bool side0toMove;
} position;
