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

#include "fen_print.h"
#include "position_sanity.h"
#include "tree_create.h"
#include "tree_search.h"

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

    if (p.side0isBlack) {
      p = rotate_position_across_central_rows(p);
    }
    print_fen(p);
  }

  return 0;
}
