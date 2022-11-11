#include <assert.h>
#include <immintrin.h>
#include <stdio.h>

#include "fen_print.h"
#include "filter_bishop.h"
#include "filter_check.h"
#include "filter_pawn.h"
#include "position_sanity.h"
#include "tree_create.h"
#include "tree_search.h"
#include "util.h"

#define TEN_THOUSAND 10000

void compute_binomials() {
  for (int n = 0; n <= NUM_SQUARES; n++) {
    binomials[n][0] = 1;
  }

  for (int k = 1; k <= MAX_BISHOPS_PSIDE; k++) {
    for (int n = k; n <= NUM_SQUARES; n++) {
      binomials[n][k] = binomials[n - 1][k - 1] + binomials[n - 1][k];
    }
  }
}

position rotate_position_across_central_rows(position p) {
  position rp;
  rp.side0isBlack = p.side0isBlack;
  rp.side0toMove = p.side0toMove;
  rp.enpassant = rotate_bitboard_across_central_rows(p.enpassant);
  rp.fixed_rooks = rotate_bitboard_across_central_rows(p.fixed_rooks);

  for (int i = 0; i < NUM_SIDES; i++) {
    rp.sides[i].pawns = rotate_bitboard_across_central_rows(p.sides[i].pawns);
    for (int j = 0; j < NUM_PIECE_TYPES; j++) {
      rp.sides[i].pieces[j] =
          rotate_bitboard_across_central_rows(p.sides[i].pieces[j]);
    }
  }
  return rp;
}

int main(int argc, char **argv) {
  compute_binomials();

  printf("Building search space\n");
  position_node *root = malloc(sizeof(position_node));
  build_sample_space(root);
  gmp_printf("Number of positions in sample space: %.10E\n",
             mpz_get_d(root->num_positions));
  // 1.2572648194E+47

  gmp_randstate_t s;
  gmp_randinit_mt(s);
  mpz_t rng;
  mpz_init(rng);

  long successes = 0;
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

#ifndef NDEBUG
    sanity_check_position(p);
#endif

    checking_info ci = validate_checks(p);
    if (ci.code != 0) {
      // if (ci.code == OBTUSE_ASPECTS) {
      //  printf("OBTUSE_ASPECTS follows:\n");
      //  print_fen(p);
      //}
      continue;
    }

    slack bas = bishop_affected_promotion_slack(p);
    if (bas.pawn_slack[0] < 0 || bas.pawn_slack[1] < 0 ||
        bas.chessmen_slack[0] < 0 || bas.chessmen_slack[1] < 0) {
      continue;
    }

    uint64_t pawns[] = {p.sides[0].pawns, p.sides[1].pawns};
    int proms[2] = {0};
    int num_capturable_chessmen[2];
    for (int i = 0; i < NUM_SIDES; i++) {
      num_capturable_chessmen[i] = _mm_popcnt_u64(p.sides[i].pawns);
      for (int j = 0; j < NUM_PIECE_TYPES_LESS_KING; j++) {
        num_capturable_chessmen[i] += _mm_popcnt_u64(p.sides[i].pieces[j]);
      }
    }
    if (!FilterPawn(pawns, p.enpassant, num_capturable_chessmen, proms)) {
      continue;
    }

    // if (p.enpassant && p.side0isBlack) {
    if (p.side0isBlack) {
      p = rotate_position_across_central_rows(p);
    }
    print_fen(p);
    successes++;
  }
  printf("successes: %ld\n", successes);

  // TODO: actually handle this
  // 1b1nn111/Bk11Q1P1/1R111Ppr/qbpr111p/11P111pP/1P11PP11/11N111bR/b11N1bKB ->
  // rook no previous square
  //
  // nbQRB1bq/1111rn11/q11n11P1/QPp11P11/1N1111Pp/1KNb1k1p/11111111/q111Rrr1 w -
  // c6 -> black would have been in check previous move by stationary pawn
  //
  // 1111b11b/111R1KrP/NR11N111/11nPp1pQ/Qq11P111/rRBB1111/1p1pPkbR/1bq11111
  //
  // 1111nK1r/b11B111q/k1pR1Pr1/bR1Q111N/Q11R1111/qR1111Bp/1RNP1bPp/n1111111
  //
  // 111b1b1Q/1Q111111/B1bk1111/111111rK/p1q1ppRr/nn1rB111/rN1RP11Q/N1N1111N
  //
  // nRN1qnk1/111111r1/111n1111/1b11111K/11PQn11r/1qB111p1/n11B11BB/1rbBBNbQ

  // likely unreachable
  // R11111K1/1r111111/1N1QP11q/11qppQ1Q/r11QbqP1/RNbN1111/111nb11R/Q1r1111k b

  // possibly unreachable
  // qB1n111b/1Q11R1Qk/p11N1Pr1/111p111n/11PK1NBN/1R1n1r1q/1BRr111r/1b1b1111

  mpf_t p_hat;
  mpf_init(p_hat);
  mpf_set_ui(p_hat, successes);
  mpf_div_ui(p_hat, p_hat, sample_size);

  mpf_t z95;
  mpf_init(z95);
  mpf_set_d(z95, 1.96);

  mpf_t standard_err;
  mpf_init(standard_err);
  mpf_set_ui(standard_err, 1);
  mpf_sub(standard_err, standard_err, p_hat);
  mpf_mul(standard_err, standard_err, p_hat);
  mpf_div_ui(standard_err, standard_err, sample_size);
  mpf_sqrt(standard_err, standard_err);
  mpf_mul(standard_err, standard_err, z95);
  mpf_t sample_space_size;
  mpf_init(sample_space_size);
  mpf_set_z(sample_space_size, root->num_positions);
  mpf_mul(standard_err, standard_err, sample_space_size);

  mpf_t pbound;
  mpf_init(pbound);
  mpf_set(pbound, sample_space_size);
  mpf_mul(pbound, pbound, p_hat);

  gmp_printf(
      "Probabilistic upperbound using 95%% C.I. on the number of positions "
      "in chess is "
      "%.2FE +- %.2FE\n",
      pbound, standard_err);

  return 0;
}
