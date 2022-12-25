#include <assert.h>
#include <immintrin.h>
#include <stdio.h>

#include "fen.h"
#include "filter_bishop.h"
#include "filter_check.h"
#include "filter_pawn.h"
#include "position_sanity.h"
#include "tree_create.h"
#include "tree_search.h"
#include "util.h"

#define NUM_THREADS_TO_SAMPLE_POSITIONS_WITH 8
#define TEN_THOUSAND                         10000

// Use samples from test file instead of generating positions ourselves.
#define USE_TEST_SAMPLES true

long        successes = 0;
FILE*       test_samples_file;
const char* test_samples_filepath
    = "/home/colin/dev/Jessica/samples/processed/sortedRng1mLegalFENsUnique";
int test_samples_file_index = 0;

typedef struct
{
    int              num_positions_to_sample_per_thread;
    position_node*   root;
    gmp_randstate_t* seed;
    pthread_mutex_t* rng_lock;
    pthread_mutex_t* onGen_lock;
    mpz_t            _rng;
} gen_positions_threading_struct;

void compute_binomials()
{
    for (int k = 0; k <= MAX_BISHOPS_PSIDE; k++)
    {
        binomials[0][k] = 0;
    }
    for (int n = 0; n <= NUM_SQUARES; n++)
    {
        binomials[n][0] = 1;
    }

    for (int k = 1; k <= MAX_BISHOPS_PSIDE; k++)
    {
        for (int n = k; n <= NUM_SQUARES; n++)
        {
            binomials[n][k] = binomials[n - 1][k - 1] + binomials[n - 1][k];
        }
    }
}

void process_positions_helper(position* p, pthread_mutex_t* onGen_lock)
{
#ifndef NDEBUG
    sanity_check_position(p);
#endif

    checking_info ci = validate_checks(p);
    if (ci.code != 0)
    {
        return;
    }

    promo_info pi = bishop_affected_promotion_info(p);
    if (pi.slack.pawn_slack[0] < 0 || pi.slack.pawn_slack[1] < 0 || pi.slack.chessmen_slack[0] < 0
        || pi.slack.chessmen_slack[1] < 0)
    {
        return;
    }

    uint64_t pawns[] = { p->sides[0].pawns, p->sides[1].pawns };
    if (!FilterPawn(pawns, p->enpassant, pi.total_captured_base_pieces, pi.promotions) != 0)
    {
        return;
    }
    // successfully generated a potentially reachable position

    if (p->enpassant && p->side0isBlack)
    {
        rotate_position_across_central_rows(p);
    }

    pthread_mutex_lock(onGen_lock);
    successes++;
    print_fen(p);
    pthread_mutex_unlock(onGen_lock);
}

void* process_positions(void* args)
{
    gen_positions_threading_struct gp = *(gen_positions_threading_struct*)args;

#if !USE_TEST_SAMPLES
    for (int i = 0; i < gp.num_positions_to_sample_per_thread; i++)
    {
        pthread_mutex_lock(gp.rng_lock);
        mpz_urandomm(gp._rng, *gp.seed, gp.root->num_positions);
        pthread_mutex_unlock(gp.rng_lock);
        position p = retrieve_position_from_rng(gp.root, gp._rng);
        process_positions_helper(&p, gp.onGen_lock);
    }

#else
    char fen[PROCESSED_FEN_MAX_LENGTH];
    while (fgets(fen, PROCESSED_FEN_MAX_LENGTH, test_samples_file) != NULL)
    {
        position p = retrieve_position_from_fen(fen);
        process_positions_helper(&p, gp.onGen_lock);
    }
#endif
    return NULL;
}

int main(int argc, char** argv)
{

    position_node* root;
#if !USE_TEST_SAMPLES
    compute_binomials();
    printf("Building search space.\n");
    root = malloc(sizeof(position_node));
    build_sample_space(root);
    gmp_printf("Number of positions in sample space: %.10E\n", mpz_get_d(root->num_positions));
    // 1.2572648194E+47
#else
    printf("Sampling from test file.\n");
    test_samples_file = fopen(test_samples_filepath, "r");
    if (test_samples_file == NULL)
    {
        printf("Error opening test sample file: %s\n", test_samples_filepath);
        return 1;
    }
#endif

    long sample_size = 0;
#if !USE_TEST_SAMPLES
    if (argc > 1)
    {
        sample_size = strtol(argv[1], NULL, 10);
    }
    else
    {
        sample_size = TEN_THOUSAND;
        printf("Using default sample size\n");
    }
#endif

    pthread_t thread_ids[NUM_THREADS_TO_SAMPLE_POSITIONS_WITH];
    long num_positions_to_sample_per_thread = sample_size / NUM_THREADS_TO_SAMPLE_POSITIONS_WITH;
    long surplus
        = sample_size - (NUM_THREADS_TO_SAMPLE_POSITIONS_WITH * num_positions_to_sample_per_thread);
    gmp_randstate_t seed;
    gmp_randinit_mt(seed);
    pthread_mutex_t rng_lock;
    pthread_mutex_t onGen_lock;
    if (pthread_mutex_init(&rng_lock, NULL) + pthread_mutex_init(&onGen_lock, NULL) != 0)
    {
        printf("Error initializing mutexes\n");
        return 1;
    }
    printf("Potentially reachable positions:\n");
    for (int i = 0; i < NUM_THREADS_TO_SAMPLE_POSITIONS_WITH; i++)
    {
        gen_positions_threading_struct gp_args = { 0 };
        gp_args.num_positions_to_sample_per_thread
            = num_positions_to_sample_per_thread + (i == 0 ? surplus : 0);

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
        gp_args.root = root;
#pragma GCC diagnostic pop

        gp_args.seed       = &seed;
        gp_args.rng_lock   = &rng_lock;
        gp_args.onGen_lock = &onGen_lock;
        pthread_create(&thread_ids[i], NULL, &process_positions, &gp_args);
    }

    for (int i = 0; i < NUM_THREADS_TO_SAMPLE_POSITIONS_WITH; i++)
    {
        pthread_join(thread_ids[i], NULL);
    }

    pthread_mutex_destroy(&rng_lock);
    pthread_mutex_destroy(&onGen_lock);

#if USE_TEST_SAMPLES
    fclose(test_samples_file);
#endif

    printf("\n");
    printf("Number of potentially reachable positions from sample: %ld\n\n", successes);

#if !USE_TEST_SAMPLES
    mpf_t p_hat;
    mpf_init(p_hat);
    mpf_set_ui(p_hat, successes);
    mpf_div_ui(p_hat, p_hat, sample_size);

    mpf_t zten10;
    mpf_init(zten10);
    mpf_set_d(zten10, 6.81);

    mpf_t standard_err;
    mpf_init(standard_err);
    mpf_set_ui(standard_err, 1);
    mpf_sub(standard_err, standard_err, p_hat);
    mpf_mul(standard_err, standard_err, p_hat);
    mpf_div_ui(standard_err, standard_err, sample_size);
    mpf_sqrt(standard_err, standard_err);
    mpf_mul(standard_err, standard_err, zten10);
    mpf_t sample_space_size;
    mpf_init(sample_space_size);
    mpf_set_z(sample_space_size, root->num_positions);
    mpf_mul(standard_err, standard_err, sample_space_size);

    mpf_t pbound;
    mpf_init(pbound);
    mpf_set(pbound, sample_space_size);
    mpf_mul(pbound, pbound, p_hat);

    gmp_printf("Probabilistic upperbound using a (1-10^10) C.I. on the number of "
               "positions in chess is %.2FE +- %.2FE\n",
               pbound,
               standard_err);
#endif

    return 0;
}
