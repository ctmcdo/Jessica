#include "chess_constants.h"
#include "position.h"
#include "prom_slack.h"
#include "util.h"

// NOTE: should add test, sum to -1
#define WHITE_SQ_MASK 0xAA55AA55AA55AA55
#define BLACK_SQ_MASK 0x55AA55AA55AA55AA

typedef struct {
  int slack[NUM_SIDES];
} bishop_slack;

bishop_slack bishop_promotion_slacks(position p) {
  int num_pawns[NUM_SIDES];
  int promotions[NUM_SIDES] = {0};
  int num_base;
  int total_base_capturable_pieces[NUM_SIDES] = {0};
  for (int i = 0; i < NUM_SIDES; i++) {
    num_pawns[i] = _mm_popcnt_u64(p.sides[i].pawns);
    for (int j = 0; j < NUM_PIECE_TYPES_LESS_KING; j++) {
      num_base = _mm_popcnt_u64(p.sides[i].pieces[j]);
      int d = num_base - BASE_PIECES[0][j];
      if (d > 0) {
        promotions[i] += d;
        num_base = BASE_PIECES[0][j];
      }
      total_base_capturable_pieces[i] += num_base;
    }
    int num_white_sq_bishops =
        _mm_popcnt_u64(p.sides[i].pieces[BISHOP] & BLACK_SQ_MASK);
    int num_black_sq_bishops =
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
