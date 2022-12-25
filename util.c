#include "chess.h"

// (row, column) to bitboard
uint64_t rcb(int row, int col)
{
    return 1UL << ((row * BOARD_SIDE_LENGTH) + col);
}

// taken from chessprogrammingwiki->flipping,
// mirroring and rotating.
uint64_t rotate_bitboard_across_central_rows(uint64_t x)
{
    return ((x << 56)) | ((x << 40) & (uint64_t)(0x00ff000000000000))
           | ((x << 24) & (uint64_t)(0x0000ff0000000000))
           | ((x << 8) & (uint64_t)(0x000000ff00000000))
           | ((x >> 8) & (uint64_t)(0x00000000ff000000))
           | ((x >> 24) & (uint64_t)(0x0000000000ff0000))
           | ((x >> 40) & (uint64_t)(0x000000000000ff00)) | ((x >> 56));
}

void rotate_position_across_central_rows(position* p)
{
    p->enpassant   = rotate_bitboard_across_central_rows(p->enpassant);
    p->fixed_rooks = rotate_bitboard_across_central_rows(p->fixed_rooks);

    for (int i = 0; i < NUM_SIDES; i++)
    {
        p->sides[i].pawns = rotate_bitboard_across_central_rows(p->sides[i].pawns);
        for (int j = 0; j < NUM_PIECE_TYPES; j++)
        {
            p->sides[i].pieces[j] = rotate_bitboard_across_central_rows(p->sides[i].pieces[j]);
        }
        p->sides[i]._fr = rotate_bitboard_across_central_rows(p->sides[i]._fr);
    }
}
