#include "prom_slack.h"
#include "chess.h"

slack promotion_slack(int pawns[NUM_SIDES],
                      int total_base_capturable_pieces[NUM_SIDES],
                      int promotions[NUM_SIDES]) {
  int total_chessmen[NUM_SIDES];
  for (int i = 0; i < NUM_SIDES; i++) {
    total_chessmen[i] =
        pawns[i] + total_base_capturable_pieces[i] + promotions[i];
  }
  slack s;
  for (int i = 0; i < NUM_SIDES; i++) {
    // Naturally, the number of promotions a side has must be less than or
    // equal to the number of pawns that that side is missing
    s.pawn_slack[i] = NUM_PAWNS_PSIDE - pawns[i] - promotions[i];

    // Opposing pawns start on opposite sides of the board s.t. if they were
    // to move forward without capturing they would eventually deadlock.
    // Therefore to promote, one pawn must either (1) be captured - unblocking
    // the file for at least one opposing pawn or (2) make a capture leaving
    // it on its promotion-rank-side of any opposition pawns on the file, or at
    // least freeing the file for the opposition pawn. If a piece is missing we
    // can assume a pawn captured it. If a pawn is missing (and is not
    // considered a promotion) we can assume a pawn captured it and freed
    // another pawn on the capture file
    int opp = (i + 1) % NUM_SIDES;
    s.chessmen_slack[i] =
        2 * (NUM_PAWNS_PSIDE - pawns[opp] - promotions[opp]) +
        (NUM_BASE_CAPTURABLE_PIECES_PSIDE - total_base_capturable_pieces[opp]) +
        (NUM_CHESSMEN_PSIDE_LESS_KING - total_chessmen[i]) - promotions[i];
  }
  return s;
}
