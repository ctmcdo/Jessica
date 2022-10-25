#pragma once

#include "chess_constants.h"

// The four constraints we adhere to
typedef struct slack {
  char pawn_slack[NUM_SIDES];
  char chessmen_slack[NUM_SIDES];
} slack;

slack promotion_slacks(char pawns[NUM_SIDES],
                       char total_base_capturable_pieces[NUM_SIDES],
                       char promotions[NUM_SIDES]);
