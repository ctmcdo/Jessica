#pragma once
#include <inttypes.h>
#include <stdbool.h>

#define NUM_SIDES 2

#define BOARD_SIDE_LENGTH 8
#define NUM_SQUARES       (BOARD_SIDE_LENGTH * BOARD_SIDE_LENGTH)

#define ENPASSANT_ROW_0INDEX    3
#define KING_HOME_COLUMN_0INDEX 3

#define NUM_PIECE_TYPES                  5
#define BISHOP                           0
#define NUM_BASE_BISHOPS_PSIDE           2
#define KNIGHT                           1
#define NUM_BASE_KNIGHTS_PSIDE           2
#define ROOK                             2
#define NUM_BASE_ROOKS_PSIDE             2
#define NUM_BASE_ROOKS                   (NUM_SIDES * NUM_BASE_ROOKS_PSIDE)
#define QUEEN                            3
#define NUM_BASE_QUEENS_PSIDE            1
#define KING                             4
#define NUM_KINGS_PSIDE                  1
#define NUM_PIECE_TYPES_LESS_KING        (NUM_PIECE_TYPES - NUM_KINGS_PSIDE)
#define NUM_BASE_CAPTURABLE_PIECES_PSIDE (BOARD_SIDE_LENGTH - NUM_KINGS_PSIDE)

#define NUM_PAWNS_PSIDE      BOARD_SIDE_LENGTH
#define MAX_PROMOTIONS_PSIDE NUM_PAWNS_PSIDE

#define NUM_PIECES_PSIDE           BOARD_SIDE_LENGTH
#define NUM_PIECES_PSIDE_LESS_KING (NUM_PIECES_PSIDE - NUM_KINGS_PSIDE)

#define NUM_CHESSMEN_PSIDE            (NUM_PAWNS_PSIDE + NUM_PIECES_PSIDE)
#define NUM_CAPTURABLE_CHESSMEN_PSIDE (NUM_CHESSMEN_PSIDE - NUM_KINGS_PSIDE)

#define MAX_BISHOPS_PSIDE (NUM_BASE_BISHOPS_PSIDE + MAX_PROMOTIONS_PSIDE)

#define NUM_FIXED_ROOK_SCENARIOS 3

#define WHITE_KINGSIDE_ROOK_BIT  0
#define WHITE_QUEENSIDE_ROOK_BIT (BOARD_SIDE_LENGTH - 1)
#define BLACK_KINGSIDE_ROOK_BIT  (BOARD_SIDE_LENGTH * (BOARD_SIDE_LENGTH - 1))
#define BLACK_QUEENSIDE_ROOK_BIT (BLACK_KINGSIDE_ROOK_BIT + (BOARD_SIDE_LENGTH - 1))

extern int BASE_PIECES[NUM_FIXED_ROOK_SCENARIOS][NUM_PIECE_TYPES_LESS_KING];

typedef struct side
{
    uint64_t pieces[NUM_PIECE_TYPES];
    uint64_t pawns;
    uint64_t _fr;
} side;

typedef struct position
{
    side     sides[NUM_SIDES];
    uint64_t enpassant;
    uint64_t fixed_rooks;
    bool     side0isBlack;
    bool     side0toMove;
} position;
