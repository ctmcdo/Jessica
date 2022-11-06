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

  gmp_printf("Probabilistic upperbound on the number of positions in chess is "
             "%.2FE +- %.2FE\n",
             pbound, standard_err);

  return 0;
}
