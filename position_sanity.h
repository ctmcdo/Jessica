#include <assert.h>
#include <nmmintrin.h>
#include <stdlib.h>

#include "chess.h"
#include "util.h"

// This isn't a filter but a sanity check. It should only
// fail if there's a bug in the creation / searching of the
// sample tree.

void sanity_check_position(position p)
{

    // no two chessmen on same square
    uint64_t checkers     = 0;
    uint64_t chessmen_sum = 0;
    for (int i = 0; i < 2; i++)
    {
        checkers += p.sides[i].pawns;
        chessmen_sum += _mm_popcnt_u64(p.sides[i].pawns);
        for (int j = 0; j < 5; j++)
        {
            checkers += p.sides[i].pieces[j];
            chessmen_sum += _mm_popcnt_u64(p.sides[i].pieces[j]);
        }
    }
    assert(chessmen_sum == _mm_popcnt_u64(checkers));

    // only one king per side
    for (int i = 0; i < 2; i++)
    {
        assert(1 == _mm_popcnt_u64(p.sides[i].pieces[4]));
    }

    // pawns aren't on first or last rows
    for (int i = 0; i < 2; i++)
    {
        if (p.sides[i].pawns != 0)
        {
            assert((p.sides[i].pawns >= rcb(1, 0)) && (p.sides[i].pawns < rcb(7, 0)));
        }
    }

    // enpassant
    // Note that we potentially flip the enpassant squares later in main.c
    // so here enpassant square can only be on the 4th rank
    int num_enpassant_pawns = _mm_popcnt_u64(p.enpassant);
    if (num_enpassant_pawns == 1)
    {
        // only on 4th rank
        assert(p.enpassant >= rcb(3, 0) && p.enpassant < rcb(3 + 1, 0));
        // and nothing behind it
        assert(0 == (checkers & (p.enpassant >> 8)));
        assert(0 == (checkers & (p.enpassant >> (2 * 8))));
        // at least one adjacent pawn
#ifndef NDEBUG
        // I use guards here to avoid an unused variable warning. Same below
        uint64_t adjacent_mask;
        if (rcb(3, 0) == p.enpassant)
        {
            adjacent_mask = rcb(3, 1);
        }
        else if (rcb(3, 7) == p.enpassant)
        {
            adjacent_mask = rcb(3, 6);
        }
        else
        {
            adjacent_mask = (p.enpassant << 1) + (p.enpassant >> 1);
        }
        assert(p.sides[1].pawns & adjacent_mask);
        // is accounted for by pawns
        assert(p.enpassant & p.sides[0].pawns);
#endif
    }
    else
    {
        assert(0 == num_enpassant_pawns);
    }

#ifndef NDEBUG
    // rooks with castling rights
    // and kings at home
    if (0 < p.sides[0]._fr)
    {
        uint64_t fixed_rooks_mask = rcb(0, 0) + rcb(0, 7);
        assert(0 == (p.sides[0]._fr & (~fixed_rooks_mask)));
        assert(0 != (p.sides[0].pieces[4] & rcb(0, 3)));
        assert(0 != (p.sides[0].pieces[2] & fixed_rooks_mask));
        // is accounted for by rooks
        assert(p.sides[0]._fr & p.sides[0].pieces[2]);
    }
    if (0 < p.sides[1]._fr)
    {
        uint64_t fixed_rooks_mask = rcb(7, 0) + rcb(7, 7);
        assert(0 == (p.sides[1]._fr & (~fixed_rooks_mask)));
        assert(0 != (p.sides[1].pieces[4] & rcb(7, 3)));
        assert(0 != (p.sides[1].pieces[2] & fixed_rooks_mask));
        assert(p.sides[1]._fr & p.sides[1].pieces[2]);
    }
#endif
}
