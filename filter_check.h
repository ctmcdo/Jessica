#include <inttypes.h>

#include "chess_constants.h"
#include "position.h"
#include "util.h"

#define KNIGHT_ON_18_ATTACKS 0x0A1100110A

#define X_ROT_45_DEG 0x0102040810204080
#define X_ROT_NEG_45_DEG 0x8040201008040201

#define NUM_SLIDING_TYPES 3

#define SOUTH_EAST_ORIGIN 1
#define SOUTH_ORIGIN (1UL << 1)
#define SOUTH_WEST_ORIGIN (1UL << 2)
#define ANY_ORIGIN (1UL << 3) - 1

#define MOORE_NEIGHBOURHOOD_SIZE 8
#define SLIDING_QUEEN_INDEX 1
#define WEST 0
#define NORTH_WEST 1
#define EAST 4
#define SOUTH_EAST 5
#define SOUTH 6
#define SOUTH_WEST 7

#define ENPASSANT_RANK 3

enum CHECK_ERROR_CODES {
  TOUCHING_KINGS = 1,
  SIDE0_CHECKED_BY_KNIGHT,                                     // 2
  SIDE0_CHECKED_BY_SLIDING_PIECE,                              // 3
  SIDE0_CHECKED_BY_PAWN,                                       // 4
  MORE_THAN_1_CHECKING_KNIGHT,                                 // 5
  MORE_THAN_1_CHECKING_PAWN,                                   // 6
  SIDE1_IN_CHECK_WHILE_ENPASSANT,                              // 7
  MORE_THAN_2_CHECKING_CHESSMEN,                               // 8
  MORE_THAN_1_CHECKING_BISHOP,                                 // 9
  CHECKED_BY_ROOK_WITH_CASTLING_RIGHTS,                        // 10
  LAYERED_SOUTH_EAST_SLIDING_ATTACK_JUST_BEHIND_CHECKING_PAWN, // 11
  LAYERED_SOUTH_WEST_SLIDING_ATTACK_JUST_BEHIND_CHECKING_PAWN, // 12
  NON_SOUTHERN_CHECKING_CHESSMAN_WITH_PAWN_ATTACK,             // 13
  CHECKED_BY_STARTING_PAWN,                                    // 14
  CHECKED_BY_PAWN_NO_PREVIOUS_SQUARE,                          // 15
  CHECKED_BY_KNIGHT_NOT_ON_PROMOTION_RANK_NO_PREVIOUS_SQUARE,  // 16
  CHECKED_BY_KNIGHT_ON_PROMOTION_RANK_NO_PREVIOUS_SQUARE,      // 17
  LAYERED_SLIDING_ATTACK_WITH_ONE_CHECKING_PIECE,
  KNIGHT_DOUBLE_CHECK_NO_INTERSECTION,
  OBTUSE_ASPECTS,
  ACUTE_SLIDING_ATTACKS_BLOCKED,
  ACUTE_SLIDING_ATTACKS_NO_INTERSECTION,
  SLIDING_ATTACK_BEHIND_CHECKING_PIECE_WHICH_JUST_MOVED
};

typedef struct checking_info {
  char code;
  char forced_capture_sq;
  char forced_promotion_sq;
  char previous_square_options;
} checking_info;

typedef struct checking_piece {
  bool type_is_queen;
  char bit;
  uint64_t ray;
  char aspect;
} checking_piece;

uint64_t get_knight_moves(char knight_bit) {
  char col = get_col_num(knight_bit);
  char max = col + 2;
  if (max > (BOARD_SIDE_LENGTH - 1)) {
    max = BOARD_SIDE_LENGTH - 1;
  }
  char min = col - 2;
  if (min < 0) {
    min = 0;
  }
  uint64_t mask = 0;
  for (int i = min; i <= max; i++) {
    mask += COLUMN_1 << i;
  }

  uint64_t knight = KNIGHT_ON_18_ATTACKS;
  if (knight_bit < 18) {
    knight = knight >> (18 - knight_bit);
  } else {
    knight = knight << (knight_bit - 18);
  }

  return knight & mask;
}

uint64_t get_occupied_squares_for_side(side s) {
  uint64_t checkers = 0;
  checkers += s.pawns;
  for (int i = 0; i < NUM_PIECE_TYPES; i++) {
    checkers += s.pieces[i];
  }

  return checkers;
}

uint64_t get_occupied_squares(position p) {
  return get_occupied_squares_for_side(p.sides[0]) +
         get_occupied_squares_for_side(p.sides[1]);
}

// taken from chessprogrammingwiki->flipping,
// mirroring and rotating.
uint64_t rotate_bitboard_across_board_center(uint64_t n) {
  const uint64_t h1 = 0x5555555555555555;
  const uint64_t h2 = 0x3333333333333333;
  const uint64_t h4 = 0x0F0F0F0F0F0F0F0F;
  const uint64_t v1 = 0x00FF00FF00FF00FF;
  const uint64_t v2 = 0x0000FFFF0000FFFF;
  n = ((n >> 1) & h1) | ((n & h1) << 1);
  n = ((n >> 2) & h2) | ((n & h2) << 2);
  n = ((n >> 4) & h4) | ((n & h4) << 4);
  n = ((n >> 8) & v1) | ((n & v1) << 8);
  n = ((n >> 16) & v2) | ((n & v2) << 16);
  n = (n >> 32) | (n << 32);
  return n;
}

uint64_t get_ray_exclusive(char exclusive_bit, char aspect) {
  if (MOORE_NEIGHBOURHOOD_SIZE / 2 <= aspect) {
    exclusive_bit = (NUM_SQUARES - 1) - exclusive_bit;
  }

  char row = get_row_num(exclusive_bit);
  char col = get_col_num(exclusive_bit);

  char org_aspect = aspect;
  aspect %= MOORE_NEIGHBOURHOOD_SIZE / 2;

  uint64_t ray;
  switch (aspect) {
  case 0:
    ray = (ROW_1 << exclusive_bit) & (ROW_1 << (row * BOARD_SIDE_LENGTH));
    break;

  case 1:
    ray = X_ROT_NEG_45_DEG << exclusive_bit;
    if (row < col) {
      char cutoff_row = row + ((BOARD_SIDE_LENGTH - 1) - col);
      ray &= (1UL << ((cutoff_row + 1) * BOARD_SIDE_LENGTH)) - 1;
    }
    break;

  case 2:
    ray = (COLUMN_1 << exclusive_bit) & -(1UL << exclusive_bit);
    break;

  case 3:
    ray = X_ROT_45_DEG;
    if (exclusive_bit < (BOARD_SIDE_LENGTH - 1)) {
      ray >>= (BOARD_SIDE_LENGTH - 1) - exclusive_bit;
    } else {
      ray <<= exclusive_bit - (BOARD_SIDE_LENGTH - 1);
    }

    if (row + col <= BOARD_SIDE_LENGTH - 1) {
      char cutoff_row = row + col;
      ray &= (1UL << (cutoff_row * BOARD_SIDE_LENGTH + 1)) - 1;
    }
    break;

  default:
    assert(false);
  }

  ray &= ~(1UL << exclusive_bit);

  if (MOORE_NEIGHBOURHOOD_SIZE / 2 <= org_aspect) {
    ray = rotate_bitboard_across_board_center(ray);
  }

  return ray;
}

uint64_t get_trimmed_ray(uint64_t ray, char bit, char aspect) {
  if (aspect < MOORE_NEIGHBOURHOOD_SIZE / 2) {
    return ray & ((1UL << bit) - 1);
  }

  return ray & ~((1UL << (bit + 1)) - 1);
}

checking_piece get_sliding_attack_from_ahead(uint64_t pot_checking_queens,
                                             uint64_t type2_bb,
                                             uint64_t blocking_chessmen,
                                             uint64_t checking_ray) {
  char blocking_tz = _tzcnt_u64(checking_ray & blocking_chessmen);
  char tz = blocking_tz;

  checking_piece cp = {0};
  cp.bit = NUM_SQUARES;
  char queen_tz = _tzcnt_u64(checking_ray & pot_checking_queens);
  if (queen_tz < tz) {
    tz = queen_tz;
    cp.type_is_queen = true;
  }

  char type2_tz = _tzcnt_u64(checking_ray & type2_bb);
  if (type2_tz < tz) {
    tz = type2_tz;
  }

  if (tz < blocking_tz) {
    cp.bit = tz;
    cp.ray = get_trimmed_ray(checking_ray, tz, 0);
  }

  return cp;
}

checking_piece get_sliding_attack_from_behind(uint64_t pot_checking_queens,
                                              uint64_t type2_bb,
                                              uint64_t blocking_chessmen,
                                              uint64_t checking_ray) {
  char blocking_lz = _lzcnt_u64(blocking_chessmen & checking_ray);
  char lz = blocking_lz;

  checking_piece cp;
  cp.bit = NUM_SQUARES;
  char queen_lz = _lzcnt_u64(pot_checking_queens & checking_ray);
  if (queen_lz < lz) {
    lz = queen_lz;
    cp.type_is_queen = true;
  }

  char type2_lz = _lzcnt_u64(type2_bb & checking_ray);
  if (type2_lz < lz) {
    lz = type2_lz;
  }

  if (lz < blocking_lz) {
    cp.bit = (NUM_SQUARES - 1) - lz;
    cp.ray = get_trimmed_ray(checking_ray, (NUM_SQUARES - 1) - lz,
                             MOORE_NEIGHBOURHOOD_SIZE / 2);
  }

  return cp;
}

checking_piece get_sliding_attack(uint64_t pot_checking_queens,
                                  uint64_t type2_bb, uint64_t blocking_chessmen,
                                  char target_bit, char aspect) {
  checking_piece cp;
  if (aspect < MOORE_NEIGHBOURHOOD_SIZE / 2) {
    cp = get_sliding_attack_from_ahead(pot_checking_queens, type2_bb,
                                       blocking_chessmen,
                                       get_ray_exclusive(target_bit, aspect));
  } else if (aspect < MOORE_NEIGHBOURHOOD_SIZE) {
    cp = get_sliding_attack_from_behind(pot_checking_queens, type2_bb,
                                        blocking_chessmen,
                                        get_ray_exclusive(target_bit, aspect));
  } else {
    assert(false);
  }
  cp.aspect = aspect;

  return cp;
}

uint64_t get_checking_pawns(uint64_t pawns, uint64_t king) {
  uint64_t checking_pawns = 0;
  char pawn_index = _tzcnt_u64(pawns);
  while (pawn_index < NUM_SQUARES) {
    char row = get_row_num(pawn_index);
    char col = get_col_num(pawn_index);
    if (row > 0) {
      if (col > 0) {
        if ((1UL << (pawn_index + BOARD_SIDE_LENGTH - 1)) & king) {
          checking_pawns += 1UL << pawn_index;
        }
      }

      if (col < BOARD_SIDE_LENGTH - 1) {
        if ((1UL << (pawn_index + BOARD_SIDE_LENGTH + 1)) & king) {
          checking_pawns += 1UL << pawn_index;
        }
      }
    }

    pawns -= 1UL << pawn_index;
    pawn_index = _tzcnt_u64(pawns);
  }

  return checking_pawns;
}

// side 0 is never to move and hence can never be in check
char validate_checks_side0(position p) {
  char king_bit = get_index_of_1st_set_bit(p.sides[0].pieces[KING]);
  if (_mm_popcnt_u64(get_knight_moves(king_bit) & p.sides[1].pieces[KNIGHT])) {
    return SIDE0_CHECKED_BY_KNIGHT;
  }

  uint64_t blocking_chessmen =
      get_occupied_squares(p) - p.sides[1].pieces[QUEEN];
  for (int i = 0; i < MOORE_NEIGHBOURHOOD_SIZE; i++) {
    char other_piece_type = 2 * ((i + 1) % 2);
    uint64_t other_piece_bb = p.sides[1].pieces[other_piece_type];
    if (get_sliding_attack(p.sides[1].pieces[QUEEN], other_piece_bb,
                           blocking_chessmen - other_piece_bb, king_bit, i)
            .bit != NUM_SQUARES) {
      return SIDE0_CHECKED_BY_SLIDING_PIECE;
    }
  }

  if (get_checking_pawns(rotate_bitboard_across_board_center(p.sides[1].pawns),
                         1UL << ((NUM_SQUARES - 1) - king_bit))) {
    return SIDE0_CHECKED_BY_PAWN;
  }

  return 0;
}

checking_info validate_checks_side1(position p) {
  checking_info ci = {0};

  char king_bit = get_index_of_1st_set_bit(p.sides[1].pieces[KING]);
  // pretend knight on opposition's king's square. The moves such a knight can
  // make are also the squares on which knights can check the king.
  uint64_t checking_knights =
      get_knight_moves(king_bit) & p.sides[0].pieces[KNIGHT];
  char num_checking_knights = _mm_popcnt_u64(checking_knights);
  // no discovery possible with 1 < knights
  if (num_checking_knights > 1) {
    ci.code = MORE_THAN_1_CHECKING_KNIGHT;
    return ci;
  }

  char num_sliding_attacks[NUM_SLIDING_TYPES] = {0};
  // sliding pieces which are attacking opposition king.
  checking_piece sliding_attacks_by_piece_type[NUM_SLIDING_TYPES] = {0};
  for (int i = 0; i < NUM_SLIDING_TYPES; i++) {
    sliding_attacks_by_piece_type[i].bit = NUM_SQUARES;
    sliding_attacks_by_piece_type[i].aspect = MOORE_NEIGHBOURHOOD_SIZE;
  }

  char sliding_attacks_by_direction[MOORE_NEIGHBOURHOOD_SIZE];
  for (int i = 0; i < MOORE_NEIGHBOURHOOD_SIZE; i++) {
    sliding_attacks_by_direction[i] = NUM_SQUARES;
  }

  uint64_t occupied_squares = get_occupied_squares(p);
  uint64_t blocking_chessmen = occupied_squares - p.sides[0].pieces[QUEEN];
  // for i in {NORTH-WEST, NORTH, NORTH-EAST, ..., WEST}
  for (int i = 0; i < MOORE_NEIGHBOURHOOD_SIZE; i++) {

    // A queen can check in any line, but depending on whether the direction (i)
    // is orthogonal or diagonal the other potentially checking piece type is
    // rook or bishop, respectively.
    char other_piece_type = 2 * ((i + 1) % 2);
    uint64_t other_piece_bb = p.sides[0].pieces[other_piece_type];
    checking_piece cp = get_sliding_attack(
        p.sides[0].pieces[QUEEN], p.sides[0].pieces[other_piece_type],
        blocking_chessmen - other_piece_bb, king_bit, i);

    // If there is a checking piece on this line
    if (cp.bit != NUM_SQUARES) {
      sliding_attacks_by_direction[i] = cp.bit;
      if (cp.type_is_queen) {
        sliding_attacks_by_piece_type[SLIDING_QUEEN_INDEX] = cp;
        num_sliding_attacks[SLIDING_QUEEN_INDEX]++;
      } else {
        sliding_attacks_by_piece_type[other_piece_type] = cp;
        num_sliding_attacks[other_piece_type]++;
      }
    }
  }

  uint64_t checking_pawns =
      get_checking_pawns(p.sides[0].pawns, p.sides[1].pieces[KING]);
  char num_checking_pawns = _mm_popcnt_u64(checking_pawns);
  if (num_checking_pawns > 1) {
    ci.code = MORE_THAN_1_CHECKING_PAWN;
    return ci;
  }

  char num_checking_chessmen = num_checking_knights +
                               num_sliding_attacks[BISHOP] +
                               num_sliding_attacks[SLIDING_QUEEN_INDEX] +
                               num_sliding_attacks[ROOK] + num_checking_pawns;
  assert(num_checking_chessmen >= 0 && num_checking_chessmen <= 10);

  if (num_checking_chessmen > 2) {
    ci.code = MORE_THAN_2_CHECKING_CHESSMEN;
    return ci;
  }

  // Every double check has to be through discovery, but two bishops can never
  // cover each other. Two rooks or two queens can though because of discovery
  // through a just-promoted pawn
  if (num_sliding_attacks[BISHOP] > 1) {
    ci.code = MORE_THAN_1_CHECKING_BISHOP;
    return ci;
  }

  if (num_checking_chessmen == 0) {
    return ci;

    // NOTE: this is missing test
  } else if (p.enpassant) {
    // Side 1 can't perform enpassant because it's in check so I don't consider
    // this position different from its non-enpassant equivalent.
    ci.code = SIDE1_IN_CHECK_WHILE_ENPASSANT;
    return ci;
  }

  if (p.sides[0].fixed_rooks) {
    // Side 0 rook hasn't moved, and side 1 king can't move into check
    if ((1UL << sliding_attacks_by_piece_type[ROOK].bit) &
        p.sides[0].fixed_rooks) {
      ci.code = CHECKED_BY_ROOK_WITH_CASTLING_RIGHTS;
      return ci;
    }
  }

  // Starting pawn hasn't moved and opposition king can't move into check
  if (checking_pawns && (king_bit < ENPASSANT_RANK * BOARD_SIDE_LENGTH)) {
    ci.code = CHECKED_BY_STARTING_PAWN;
    return ci;
  }

  if (num_checking_chessmen == 1) {
    return ci;
  }

  // else we have 2 checking chessmen
  if (checking_knights) {
    return ci;
  }

  /*
  char aspect_diff = 0;
  char aspect_total = 0;
  for (int i = 0; i < NUM_SLIDING_TYPES; i++) {
    if (sliding_attacks[i].aspect != MOORE_NEIGHBOURHOOD_SIZE) {
      if (aspect_diff < 0) {
        aspect_diff += sliding_attacks_by_piece_type[i].aspect;
      } else {
        aspect_diff -= sliding_attacks_by_piece_type[i].aspect;
      }
      aspect_total += sliding_attacks_by_piece_type[i].aspect;
    }
  }
  aspect_diff = abs(aspect_diff);
  if (aspect_diff != 1 &&
      !(aspect_diff == SOUTH_WEST && aspect_total == SOUTH_WEST)) {
    ci.code = OBTUSE_ASPECTS;
    return ci;
  }
  */
  return ci;
}

checking_info validate_checks(position p) {
  checking_info ci = {0};
  char k0 = get_index_of_1st_set_bit(p.sides[0].pieces[KING]);
  char k1 = get_index_of_1st_set_bit(p.sides[1].pieces[KING]);
  if (abs(get_row_num(k0) - get_row_num(k1)) < 2 &&
      abs(get_col_num(k0) - get_col_num(k1)) < 2) {
    ci.code = TOUCHING_KINGS;
    return ci;
  }

  char code = validate_checks_side0(p);
  if (code) {
    ci.code = code;
    return ci;
  }

  return validate_checks_side1(p);
}
