#include "chess.h"
#include "prom_slack.h"

#define WHITE_SQ_MASK 0xAA55AA55AA55AA55
#define BLACK_SQ_MASK 0x55AA55AA55AA55AA

slack bishop_affected_promotion_slack(position p) {
  int num_pawns[NUM_SIDES];
  int num_base;
  int promotions[NUM_SIDES] = {0};
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
        _mm_popcnt_u64(p.sides[i].pieces[BISHOP] & WHITE_SQ_MASK);
    int num_black_sq_bishops =
        _mm_popcnt_u64(p.sides[i].pieces[BISHOP] & BLACK_SQ_MASK);
    if (_mm_popcnt_u64(p.sides[i].pieces[BISHOP]) &&
        (num_black_sq_bishops == 0 || num_white_sq_bishops == 0)) {
      promotions[i] += 1;
    }
  }

  return promotion_slack(num_pawns, total_base_capturable_pieces, promotions);
}
