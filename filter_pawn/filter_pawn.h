#include <inttypes.h>

#define NUM_SIDES 2

#ifdef _cplusplus
extern "C" {
#endif
bool FilterPawn(uint64_t pawns[NUM_SIDES], uint64_t enpassant,
                int num_captuable_chessmen[NUM_SIDES],
                int min_promotions[NUM_SIDES]);
#ifdef _cplusplus
}
#endif
