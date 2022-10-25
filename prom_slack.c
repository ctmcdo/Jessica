#include "prom_slack.h"
#include "chess_constants.h"

slack promotion_slacks(char pawns[NUM_SIDES],
                       char total_base_capturable_pieces[NUM_SIDES],
                       char promotions[NUM_SIDES]) {
  slack s;
  for (int i = 0; i < NUM_SIDES; i++) {
    // Naturally, the number of promotions a side has must be less than or equal
    // to the number of pawns that that side is missing
    s.pawn_slack[i] = NUM_PAWNS_PSIDE - pawns[i] - promotions[i];

    // Opposing pawns start on opposite sides of the board s.t. if they were
    // to move forward without capturing they would eventually deadlock.
    // Therefore to promote, one pawn must either (1) be captured - unblocking
    // the file for at least one opposing pawn or (2) make a capture leaving it
    // on its promotion-rank-side of any opposition pawns on the file. If a
    // piece is missing we can assume a pawn captured it. If a pawn is missing
    // (and is not considered a promotion) we can assume a pawn captured it and
    // freed another pawn on the capture file
    int opp = (i + 1) % NUM_SIDES;
    s.chessmen_slack[i] =
        2 * (NUM_PAWNS_PSIDE - pawns[opp] - promotions[opp]) +
        (NUM_BASE_CAPTURABLE_PIECES_PSIDE - total_base_capturable_pieces[opp]) -
        promotions[i];
  }
  return s;
}
