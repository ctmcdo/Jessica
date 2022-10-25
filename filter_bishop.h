#include "chess_constants.h"
#include "position.h"
#include "prom_slack.h"
#include "util.h"

// NOTE: should add test, sum to -1
#define WHITE_SQ_MASK 0xAA55AA55AA55AA55
#define BLACK_SQ_MASK 0x55AA55AA55AA55AA

typedef struct {
  char slack[NUM_SIDES];
} bishop_slack;

bishop_slack bishop_promotion_slacks(position p) {
  char num_fixed_rooks[NUM_SIDES];
  num_fixed_rooks[0] = _mm_popcnt_u64(p.sides[0].fixed_rooks);
  num_fixed_rooks[1] = _mm_popcnt_u64(p.sides[1].fixed_rooks);

  char num_pawns[NUM_SIDES];
  char promotions[NUM_SIDES] = {0};
  char num_base;
  char total_base_capturable_pieces[NUM_SIDES] = {0};
  for (int i = 0; i < NUM_SIDES; i++) {
    num_pawns[i] = _mm_popcnt_u64(p.sides[i].pawns);
    for (int j = 0; j < NUM_PIECE_TYPES_LESS_KING; j++) {
      num_base = _mm_popcnt_u64(p.sides[i].pieces[j]);
      char d = num_base - BASE_PIECES[0][j];
      if (d > 0) {
        promotions[i] += d;
        num_base = BASE_PIECES[0][j];
      }
      total_base_capturable_pieces[i] += num_base;
    }
    char num_white_sq_bishops =
        _mm_popcnt_u64(p.sides[i].pieces[BISHOP] & BLACK_SQ_MASK);
    char num_black_sq_bishops =
        _mm_popcnt_u64(p.sides[i].pieces[BISHOP] & BLACK_SQ_MASK);
    if (num_black_sq_bishops == 0 || num_white_sq_bishops == 0) {
      promotions[i] += 1;
    }
  }

  slack s =
      promotion_slacks(num_pawns, total_base_capturable_pieces, promotions);
  bishop_slack bs;
  bs.slack[0] = min3(s.pawn_slack[0], s.chessmen_slack[0], s.chessmen_slack[1]);
  bs.slack[1] = min3(s.pawn_slack[1], s.chessmen_slack[0], s.chessmen_slack[1]);

  return bs;
}
