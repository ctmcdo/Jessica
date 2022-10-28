#include <assert.h>
#include <gmp.h>
#include <immintrin.h>
#include <inttypes.h>
#include <math.h>
#include <nmmintrin.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "filter_bishop.h"
#include "filter_check.h"
#include "filter_pawn.h"
#include "position.h"
#include "position_sanity.h"
#include "tree_create.h"
#include "tree_search.h"
#include "util.h"

#define NUM_THREADS 8

#define ROOK0_BIT 0
#define ROOK1_BIT (ROOK0_BIT + (BOARD_SIDE_LENGTH - 1))
#define ROOK2_BIT (BOARD_SIDE_LENGTH * (BOARD_SIDE_LENGTH - 1))
#define ROOK3_BIT (ROOK2_BIT + (BOARD_SIDE_LENGTH - 1))
#define a_ASCII_decimal 97

#define TEN_THOUSAND 10000

// TODO: something related to opposite side to move in check? Take a look at
// previous archive

// TODO: use int instead of char in for loops
// and also in functions etc, shouldn't use chars

// TODO: rename sample tree to sample space

// typedef struct {
//  rgns[];
//  sample_size;
//  successes;
//};

// void *filter_samples(void *) {}

long wcount = 0;
long bcount = 0;
long side0toMoveCount = 0;
long side1toMoveCount = 0;
long kcount = 0;
long Kcount = 0;
long qcount = 0;
long Qcount = 0;
long pawnswcount = 0;
long pawnsbcount = 0;

position rotate_position_across_central_rows(position p) {
  position rp;
  rp.side0isBlack = p.side0isBlack;
  rp.enpassant = rotate_bitboard_across_central_rows(p.enpassant);
  for (int i = 0; i < NUM_SIDES; i++) {
    rp.sides[i].fixed_rooks =
        rotate_bitboard_across_central_rows(p.sides[i].fixed_rooks);

    rp.sides[i].pawns = rotate_bitboard_across_central_rows(p.sides[i].pawns);

    for (int j = 0; j < NUM_PIECE_TYPES; j++) {
      rp.sides[i].pieces[j] =
          rotate_bitboard_across_central_rows(p.sides[i].pieces[j]);
    }
  }
  return rp;
}

void fen_helper(uint64_t chessmen, char symbol, int squares[NUM_SQUARES]) {
  int tz = _tzcnt_u64(chessmen);
  while (tz != NUM_SQUARES) {
    squares[tz] = symbol;
    chessmen ^= 1UL << tz;
    tz = _tzcnt_u64(chessmen);
  }
}

void print_fen(position p) {
  int squares[NUM_SQUARES];
  for (int i = 0; i < NUM_SQUARES; i++) {
    squares[i] = '1';
  }

  fen_helper(p.sides[p.side0isBlack].pawns, 'P', squares);
  fen_helper(p.sides[p.side0isBlack].pieces[BISHOP], 'B', squares);
  fen_helper(p.sides[p.side0isBlack].pieces[KNIGHT], 'N', squares);
  fen_helper(p.sides[p.side0isBlack].pieces[ROOK], 'R', squares);
  fen_helper(p.sides[p.side0isBlack].pieces[QUEEN], 'Q', squares);
  fen_helper(p.sides[p.side0isBlack].pieces[KING], 'K', squares);

  fen_helper(p.sides[!p.side0isBlack].pawns, 'p', squares);
  fen_helper(p.sides[!p.side0isBlack].pieces[BISHOP], 'b', squares);
  fen_helper(p.sides[!p.side0isBlack].pieces[KNIGHT], 'n', squares);
  fen_helper(p.sides[!p.side0isBlack].pieces[ROOK], 'r', squares);
  fen_helper(p.sides[!p.side0isBlack].pieces[QUEEN], 'q', squares);
  fen_helper(p.sides[!p.side0isBlack].pieces[KING], 'k', squares);

  for (int i = BOARD_SIDE_LENGTH - 1; i >= 0; i--) {
    for (int j = BOARD_SIDE_LENGTH - 1; j >= 0; j--) {
      printf("%c", squares[i * BOARD_SIDE_LENGTH + j]);
    }
    if (i > 0) {
      printf("/");
    }
  }
  printf(" ");

  char colour = 'w';
  if (p.side0isBlack) {
    colour = 'b';
  }
  printf("%c ", colour);

  bool cr0 = p.sides[0].fixed_rooks & (1UL << ROOK0_BIT);
  bool cr1 = p.sides[0].fixed_rooks & (1UL << ROOK1_BIT);
  bool cr2 = p.sides[1].fixed_rooks & (1UL << ROOK2_BIT);
  bool cr3 = p.sides[1].fixed_rooks & (1UL << ROOK3_BIT);
  if (!cr0 && !cr1 && !cr2 && !cr3) {
    printf("-");

  } else {
    if (cr0) {
      printf("K");
      Kcount += 1;
    }
    if (cr1) {
      printf("Q");
      Qcount += 1;
    }
    if (cr2) {
      printf("k");
      kcount += 1;
    }
    if (cr3) {
      printf("q");
      qcount += 1;
    }
  }
  printf(" ");

  if (p.enpassant) {
    char ep = get_index_of_1st_set_bit(p.enpassant);

    char column_letter =
        a_ASCII_decimal + (BOARD_SIDE_LENGTH - 1) - (ep % BOARD_SIDE_LENGTH);
    assert(a_ASCII_decimal <= column_letter && column_letter <= 104);
    printf("%c", column_letter);

    char row_num = 1 + ((int)(ep / BOARD_SIDE_LENGTH));
    assert(row_num == 4 || row_num == 5);
    printf("%d", row_num);

  } else {
    printf("-");
  }

  printf("\n");
}

int main(int argc, char **argv) {
  compute_binomials();

  printf("Building search space\n");
  position_node *root = malloc(sizeof(position_node));
  build_sample_space(root);
  gmp_printf("Number of positions in sample space: %.10E\n",
             mpz_get_d(root->num_positions));
  // 1.2563633941E+47

  gmp_randstate_t s;
  gmp_randinit_mt(s);
  mpz_t rng;
  mpz_init(rng);

  uint64_t successes = 0;
  long sample_size;
  if (argc > 1) {
    sample_size = strtol(argv[1], NULL, 10);
  } else {
    sample_size = TEN_THOUSAND;
    printf("Using default sample size of ten thousand\n");
  }
  for (int i = 0; i < sample_size; i++) {
    mpz_urandomm(rng, s, root->num_positions);
    position p = retrieve_position(root, rng);

    sanity_check_position(p);

    /*
    checking_info ci = validate_checks(p);
    char check_code = ci.code;
    if (check_code != 0) {
      // printf("Invalid check code %d: ", check_code);

      p.side0isBlack = i % NUM_SIDES == 1;
      if (p.side0isBlack) {
        p = rotate_position_across_central_rows(p);
      }
      // print_fen(p);
      continue;
    }

    bishop_slack bs = bishop_promotion_slacks(p);
    if (bs.slack[0] < 0 || bs.slack[1] < 0) {
      // printf("Invalid promotions code %d: ", promotions_code);

      p.side0isBlack = i % NUM_SIDES == 1;
      if (p.side0isBlack) {
        p = rotate_position_across_central_rows(p);
      }
      // print_fen(p);
      continue;
    }

    // now we assume every piece we take is not a pawn
    int pawn_cost_code = validate_pawn_cost(p, bs.slack[0], bs.slack[1]);
    if (pawn_cost_code != 0) {
      // printf("Invalid pawn cost code %d: ", pawn_cost_code);

      p.side0isBlack = i % NUM_SIDES == 1;
      if (p.side0isBlack) {
        p = rotate_position_across_central_rows(p);
      }
      // print_fen(p);
      continue;
    }
    */

    successes++;

    if (p.side0isBlack) {
      bcount += 1;
    } else {
      wcount += 1;
    }

    if (p.side0toMove) {
      side0toMoveCount += 1;
    } else {
      side1toMoveCount += 1;
    }

    p.side0isBlack = i % NUM_SIDES == 1;
    if (p.side0isBlack) {
      p = rotate_position_across_central_rows(p);
      pawnsbcount += _mm_popcnt_u64(p.sides[0].pawns);
      pawnswcount += _mm_popcnt_u64(p.sides[1].pawns);
    } else {
      pawnswcount += _mm_popcnt_u64(p.sides[0].pawns);
      pawnsbcount += _mm_popcnt_u64(p.sides[1].pawns);
    }
    print_fen(p);
  }

  mpf_t p_hat;
  mpf_init(p_hat);
  mpf_set_ui(p_hat, successes);
  mpf_div_ui(p_hat, p_hat, sample_size);

  mpf_t z995;
  mpf_init(z995);
  mpf_set_d(z995, 2.807);

  mpf_t standard_err;
  mpf_init(standard_err);
  mpf_set_ui(standard_err, 1);
  mpf_sub(standard_err, standard_err, p_hat);
  mpf_mul(standard_err, standard_err, p_hat);
  mpf_div_ui(standard_err, standard_err, sample_size);
  mpf_sqrt(standard_err, standard_err);
  mpf_mul(standard_err, standard_err, z995);
  mpf_t sample_space_size;
  mpf_init(sample_space_size);
  mpf_set_z(sample_space_size, root->num_positions);
  mpf_mul(standard_err, standard_err, sample_space_size);

  mpf_t pbound;
  mpf_init(pbound);
  mpf_set(pbound, sample_space_size);
  mpf_mul(pbound, pbound, p_hat);

  gmp_printf("Probabilistic upperbound (using confidence interval of 99.5%) on "
             "the number of positions in chess is "
             "%.2FE +- %.2FE\n",
             pbound, standard_err);

  printf("bcount: %ld\n", bcount);
  printf("wcount: %ld\n", bcount);
  printf("side0toMoveCount: %ld\n", side0toMoveCount);
  printf("side1toMoveCount: %ld\n", side1toMoveCount);

  printf("kcount: %ld\n", kcount);
  printf("Kcount: %ld\n", Kcount);
  printf("qcount: %ld\n", qcount);
  printf("Qcount: %ld\n", Qcount);

  printf("pawnswcount: %ld\n", pawnswcount);
  printf("pawnsbcount: %ld\n", pawnsbcount);

  return 0;
}
