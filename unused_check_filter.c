#include <assert.h>
#include <immintrin.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdlib.h>

#include "check.h"
#include "chess.h"
#include "filter_common.h"
#include "position.h"
#include "util.h"

#define MOORE_NEIGHBOURHOOD_SIZE 8
#define SLIDING_QUEEN_INDEX 1
#define WEST 0
#define NORTH_WEST 1
#define EAST 4
#define SOUTH_EAST 5
#define SOUTH 6
#define SOUTH_WEST 7

#define ENPASSANT_RANK 3

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

char get_previous_pawn_square_possibilities(char pbit, char kbit,
                                            uint64_t occupied_squares) {
  assert(BOARD_SIDE_LENGTH <= pbit);
  assert(pbit < BOARD_SIDE_LENGTH * BOARD_SIDE_LENGTH);
  assert(kbit < BOARD_SIDE_LENGTH * BOARD_SIDE_LENGTH);

  char pcol = get_col_num(pbit);
  char options = 0;
  uint64_t south = 1UL << (pbit - BOARD_SIDE_LENGTH);
  if (((~occupied_squares) & south) &&
      (get_checking_pawns(south, 1UL << kbit) == 0)) {
    options = FILE_ORIGIN;
  }

  uint64_t east = 1UL << (pbit - (BOARD_SIDE_LENGTH + 1));
  if ((pcol > 0) && (get_checking_pawns(east, 1UL << kbit) == 0)) {
    if ((~occupied_squares) & east) {
      options += SOUTH_EAST_ORIGIN;
    }
  }

  uint64_t west = 1UL << (pbit - (BOARD_SIDE_LENGTH - 1));
  if ((pcol < BOARD_SIDE_LENGTH - 1) &&
      (get_checking_pawns(west, 1UL << kbit) == 0)) {
    if ((~occupied_squares) & west) {
      options += SOUTH_WEST_ORIGIN;
    }
  }

  return options;
}

checking_info validate_checks_side1(position p) {
  checking_info ci = {0};

  char king_bit = get_index_of_1st_set_bit(p.sides[1].pieces[KING]);
  // pretend knight on opposition's king's square. The moves such a knight can
  // make are also the squares on which knights can check the king.
  uint64_t checking_knights =
      get_knight_moves(king_bit) & p.sides[0].pieces[KNIGHT];
  char num_checking_knights = _mm_popcnt_u64(checking_knights);
  // no discovery possible with 2>= knights
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

    // TODO: add test
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

  if (checking_pawns) {
    // Starting pawn hasn't moved and opposition king can't move into check
    if (king_bit < ENPASSANT_RANK * BOARD_SIDE_LENGTH) {
      ci.code = CHECKED_BY_STARTING_PAWN;
      return ci;
    }

    char pawn_bit = get_index_of_1st_set_bit(checking_pawns);
    char pawn_col = get_col_num(pawn_bit);
    char king_col = get_col_num(king_bit);
    if (pawn_col < king_col) {
      checking_piece cp = get_sliding_attack(
          p.sides[0].pieces[QUEEN], p.sides[0].pieces[BISHOP],
          blocking_chessmen - p.sides[0].pieces[BISHOP], pawn_bit, SOUTH_EAST);

      // bishop/queen on diagonal just before pawn so no previous pawn square
      if (cp.bit == (pawn_bit - (BOARD_SIDE_LENGTH + 1))) {
        ci.code = LAYERED_SOUTH_EAST_SLIDING_ATTACK_JUST_BEHIND_CHECKING_PAWN;
        return ci;

        // otherwise pawn could have come from south-east
      } else if (cp.bit != NUM_SQUARES) {
        ci.forced_capture = pawn_bit;
        ci.checking_pawn_origin_options = SOUTH_EAST_ORIGIN;
      }

    } else {
      checking_piece cp = get_sliding_attack(
          p.sides[0].pieces[QUEEN], p.sides[0].pieces[BISHOP],
          blocking_chessmen - p.sides[0].pieces[BISHOP], pawn_bit, SOUTH_WEST);

      // bishop/queen on diagonal just before pawn so no previous pawn square
      if (cp.bit == (pawn_bit - (BOARD_SIDE_LENGTH - 1))) {
        ci.code = LAYERED_SOUTH_WEST_SLIDING_ATTACK_JUST_BEHIND_CHECKING_PAWN;
        return ci;

        // otherwise pawn could have come from south-west
      } else if (cp.bit != NUM_SQUARES) {
        ci.forced_capture = pawn_bit;
        ci.checking_pawn_origin_options = SOUTH_WEST_ORIGIN;
      }
    }

    char pawn_rank = get_row_num(pawn_bit);
    // if no sliding attack on diagonals behind pawn
    if (ci.checking_pawn_origin_options == 0) {
      ci.checking_pawn_origin_options = get_previous_pawn_square_possibilities(
          pawn_bit, king_bit, occupied_squares);

      // if pawn has no possible previous square
      if (ci.checking_pawn_origin_options == 0) {
        ci.code = CHECKED_BY_PAWN_NO_PREVIOUS_SQUARE;
        return ci;
      } else if (!(ci.checking_pawn_origin_options & FILE_ORIGIN)) {
        ci.forced_capture = pawn_bit;
      }
    }

    // Only discovered double attack with pawn is if rook or queen on same file
    // as king with lower rank.
    if (num_checking_chessmen > 1) {
      char pawn_previous_square =
          (pawn_rank - 1) * BOARD_SIDE_LENGTH + king_col;
      if (sliding_attacks_by_direction[SOUTH] == NUM_SQUARES ||
          sliding_attacks_by_direction[SOUTH] >= pawn_previous_square) {
        char ci.code = NON_SOUTHERN_CHECKING_CHESSMAN_WITH_PAWN_ATTACK;
      } else {
        ci.forced_capture = pawn_bit;
        if (pawn_col > king_col) {
          ci.origin_options = SOUTH_WEST;
        } else {
          ci.origin_options = SOUTH_EAST
        }
      }
    }

    return ci;
  }

  if (num_checking_chessmen == 1) {
    // TODO: check escape square or discovery
    //
    if (checking_knights) {
      char nb = get_index_of_1st_set_bit(checking_knights);
      if ((get_knight_moves(nb) & (~occupied_squares)) == 0) {
        if (nb < (BOARD_SIDE_LENGTH * (BOARD_SIDE_LENGTH - 1))) {
          ci.code = CHECKED_BY_KNIGHT_NOT_ON_PROMOTION_RANK_NO_PREVIOUS_SQUARE;
          return ci;
        }

        ci.promoted_piece = checking_knights;
        // the king_bit shouldn't matter for this call
        ci.promotion_origin_options = get_previous_pawn_square_possibilities(
            nb, king_bit, occupied_squares);
        if (ci.promotion_origin_options == 0) {
          ci.code = CHECKED_BY_KNIGHT_ON_PROMOTION_RANK_NO_PREVIOUS_SQUARE;
        }
      }

      return ci;
    }

    char aspect = 0;
    char sbit = 0;
    uint64_t ray = 0;
    for (int i = 0; i < NUM_SLIDING_TYPES; i++) {
      aspect += sliding_attacks_by_piece_type[i].aspect;
      sbit += sliding_attacks_by_piece_type[i].bit;
      ray += sliding_attacks_by_piece_type[i].ray;
    }
    aspect -= 2 * MOORE_NEIGHBOURHOOD_SIZE;
    sbit -= 2 * NUM_SQUARES;

    // if a discovery is possible, then it's preferable because it shouldn't
    // require a capture
    uint64_t knights = p.sides[1].pieces[KNIGHT];
    while (knights) {
    }

    // 2. look for an escape square.

    uint64_t other_piece_bb = p.sides[0].pieces[2 * ((aspect + 1) % 2)];
    checking_piece cp =
        get_sliding_attack(p.sides[0].pieces[QUEEN], other_piece_bb,
                           blocking_chessmen - other_piece_bb, sbit, aspect);
    // is there a sliding attack mounted behind 's'?
    if (cp.bit != NUM_SQUARES) {

      bool found_previous_sq = false;
      for (int i = -2; i < 2; i += 2) {
      }

      // TODO: or via pawn discovery

      // TODO: now there has to be an escape square where we weren't checking
      // the king previously, and we have to communicate that this requires a
      // capture. The escape square can potentially be on the same line as
      // checking line, given there's space on the line.
      ci.code = LAYERED_SLIDING_ATTACK_WITH_ONE_CHECKING_PIECE;
    }

    return ci;
  } // else we have 2 checking chessmen

  if (checking_knights) {
    uint64_t ray = sliding_attacks_by_piece_type[BISHOP].ray +
                   sliding_attacks_by_piece_type[SLIDING_QUEEN_INDEX].ray +
                   sliding_attacks_by_piece_type[ROOK].ray;
    if (!(get_knight_moves(get_index_of_1st_set_bit(checking_knights)) & ray)) {
      ci.code = KNIGHT_DOUBLE_CHECK_NO_INTERSECTION;
    }

    return ci;
  }

  // Two rooks are only possible if one rook is on same file as king, king is on
  // 8th rank, and pawn on 7th rank captures and promotes adjacent to king.
  // if (num_sliding_attacks[ROOK] > 1) {
  //   uint64_t rooks = sliding_attacks[ROOK];
  //   char rbit1 = get_index_of_1st_set_bit(rooks);
  //   rooks ^= 1UL << rbit1;
  //   char rbit2 = get_index_of_1st_set_bit(rooks);
  // }
  // (num_sliding_attacks[SLIDING_QUEEN_INDEX] > 1);

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

  checking_piece p0, p1;
  if (num_sliding_attacks[BISHOP] && num_sliding_attacks[ROOK]) {
    typedef char (*offset_func)(char);
    offset_func f = &get_row_num;
    if (get_row_num(sliding_attacks_by_piece_type[ROOK].bit) ==
        get_row_num(king_bit)) {
      f = &get_col_num;
    }
    char ka = f(king_bit);
    char ba = f(sliding_attacks_by_piece_type[BISHOP].bit);
    char ra = f(sliding_attacks_by_piece_type[ROOK].bit);

    if (abs(ra - ka) < abs(ba - ka)) {
      p0 = sliding_attacks_by_piece_type[ROOK];
      p1 = sliding_attacks_by_piece_type[BISHOP];
    } else {
      p0 = sliding_attacks_by_piece_type[BISHOP];
      p1 = sliding_attacks_by_piece_type[ROOK];
    }

  } else {
    p1 = sliding_attacks_by_piece_type[SLIDING_QUEEN_INDEX];
    if (num_sliding_attacks[BISHOP]) {
      p0 = sliding_attacks_by_piece_type[BISHOP];
    } else {
      p0 = sliding_attacks_by_piece_type[ROOK];
    }
  }

  char join_dir;
  if (p0.aspect == WEST && p1.aspect == SOUTH_WEST) {
    join_dir = SOUTH;
  } else if (p0.aspect == SOUTH_WEST && p1.aspect == WEST) {
    join_dir = NORTH_WEST;
  } else {
    if (p0.aspect < p1.aspect) {
      join_dir =
          (p0.aspect + MOORE_NEIGHBOURHOOD_SIZE + 2) % MOORE_NEIGHBOURHOOD_SIZE;
    } else {
      join_dir =
          (p0.aspect + MOORE_NEIGHBOURHOOD_SIZE - 2) % MOORE_NEIGHBOURHOOD_SIZE;
    }
  }

  uint64_t join_to_edge = get_ray_exclusive(p0.bit, join_dir);
  uint64_t intersection = join_to_edge & p1.ray;
  if (intersection) {
    if (!(get_trimmed_ray(join_to_edge, get_index_of_1st_set_bit(intersection),
                          join_dir) &
          blocking_chessmen)) {

      char other = 2 * ((p0.aspect + 1) % 2);
      checking_piece cp = get_sliding_attack(
          p.sides[0].pieces[QUEEN], p.sides[0].pieces[other],
          blocking_chessmen - p.sides[0].pieces[other], p0.bit, p0.aspect);
      if (cp.bit != NUM_SQUARES) {
        ci.code = SLIDING_ATTACK_BEHIND_CHECKING_PIECE_WHICH_JUST_MOVED;
      }

      return ci;
    }

    ci.code = ACUTE_SLIDING_ATTACKS_BLOCKED;
    return ci;
  }

  ci.code = ACUTE_SLIDING_ATTACKS_NO_INTERSECTION;
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
