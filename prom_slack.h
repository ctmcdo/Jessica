#pragma once
#include "chess.h"

// The four constraints we adhere to
typedef struct {
  int pawn_slack[NUM_SIDES];
  int chessmen_slack[NUM_SIDES];
} slack;

slack promotion_slack(int pawns[NUM_SIDES],
                      int total_base_capturable_pieces[NUM_SIDES],
                      int promotions[NUM_SIDES]);
