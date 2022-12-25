#include <ctype.h>

#include "chess.h"
#include "util.h"

#define a_ASCII_decimal       97
#define NUM_CORNER_SQUARES    4
#define PRINT_FORWARD_SLASHES true

void fen_helper(uint64_t chessmen, char symbol, int squares[NUM_SQUARES])
{
    int tz = _tzcnt_u64(chessmen);
    while (tz != NUM_SQUARES)
    {
        squares[tz] = symbol;
        chessmen ^= 1UL << tz;
        tz = _tzcnt_u64(chessmen);
    }
}

char get_index_of_1st_set_bit(uint64_t n)
{
    for (int i = 0; i < NUM_SQUARES; i++)
    {
        if (n & 1)
            return i;
        n = n >> 1;
    }
    return -1;
}

void print_fen(position* p)
{
    int squares[NUM_SQUARES];
    for (int i = 0; i < NUM_SQUARES; i++)
    {
        squares[i] = '1';
    }

    fen_helper(p->sides[p->side0isBlack].pawns, 'P', squares);
    fen_helper(p->sides[p->side0isBlack].pieces[BISHOP], 'B', squares);
    fen_helper(p->sides[p->side0isBlack].pieces[KNIGHT], 'N', squares);
    fen_helper(p->sides[p->side0isBlack].pieces[ROOK], 'R', squares);
    fen_helper(p->sides[p->side0isBlack].pieces[QUEEN], 'Q', squares);
    fen_helper(p->sides[p->side0isBlack].pieces[KING], 'K', squares);

    fen_helper(p->sides[!p->side0isBlack].pawns, 'p', squares);
    fen_helper(p->sides[!p->side0isBlack].pieces[BISHOP], 'b', squares);
    fen_helper(p->sides[!p->side0isBlack].pieces[KNIGHT], 'n', squares);
    fen_helper(p->sides[!p->side0isBlack].pieces[ROOK], 'r', squares);
    fen_helper(p->sides[!p->side0isBlack].pieces[QUEEN], 'q', squares);
    fen_helper(p->sides[!p->side0isBlack].pieces[KING], 'k', squares);

    for (int i = BOARD_SIDE_LENGTH - 1; i >= 0; i--)
    {
        for (int j = BOARD_SIDE_LENGTH - 1; j >= 0; j--)
        {
            printf("%c", squares[i * BOARD_SIDE_LENGTH + j]);
        }
#if PRINT_FORWARD_SLASHES
        if (i > 0)
        {
            printf("/");
        }
#endif
    }
    printf(" ");

    char colour = 'w';
    if (p->side0isBlack == p->side0toMove)
    {
        colour = 'b';
    }
    printf("%c ", colour);

    bool K = p->fixed_rooks & (1UL << WHITE_KINGSIDE_ROOK_BIT);
    bool Q = p->fixed_rooks & (1UL << WHITE_QUEENSIDE_ROOK_BIT);
    bool k = p->fixed_rooks & (1UL << BLACK_KINGSIDE_ROOK_BIT);
    bool q = p->fixed_rooks & (1UL << BLACK_QUEENSIDE_ROOK_BIT);
    if (!K && !Q && !k && !q)
    {
        printf("-");
    }
    else
    {
        if (K)
        {
            printf("K");
        }
        if (Q)
        {
            printf("Q");
        }
        if (k)
        {
            printf("k");
        }
        if (q)
        {
            printf("q");
        }
    }
    printf(" ");

    if (p->enpassant)
    {
        char ep = get_index_of_1st_set_bit(p->enpassant);

        char column_letter = a_ASCII_decimal + (BOARD_SIDE_LENGTH - 1) - (ep % BOARD_SIDE_LENGTH);
        assert(a_ASCII_decimal <= column_letter && column_letter <= 104);
        printf("%c", column_letter);

        char row_num = 1 + ((int)(ep / BOARD_SIDE_LENGTH));
        assert(row_num == 4 || row_num == 5);
        if (row_num == 4)
        {
            printf("3");
        }
        else
        {
            printf("6");
        }
    }
    else
    {
        printf("-");
    }

    printf("\n");
}

position retrieve_position_from_fen(char fen[PROCESSED_FEN_MAX_LENGTH])
{
    int  curr        = (NUM_SQUARES - 1) + 2;
    bool blackToMove = fen[curr] == 'b';
    assert(blackToMove || fen[curr] == 'w');
    position p    = { 0 };
    p.side0toMove = !blackToMove;

    curr += 2;
    if (fen[curr] == '-')
    {
        curr += 2;
    }
    else
    {
        for (int i = 0; i < NUM_CORNER_SQUARES; i++)
        {
            switch (fen[curr])
            {
            case 'K': p.sides[0]._fr += 1UL << (WHITE_KINGSIDE_ROOK_BIT); break;
            case 'Q': p.sides[0]._fr += 1UL << (WHITE_QUEENSIDE_ROOK_BIT); break;
            case 'k': p.sides[1]._fr += 1UL << (BLACK_KINGSIDE_ROOK_BIT); break;
            case 'q': p.sides[1]._fr += 1UL << (BLACK_QUEENSIDE_ROOK_BIT); break;
            default: printf("Unexpected castling option\n"); exit(1);
            }
            curr++;
            if (fen[curr] == ' ')
            {
                curr++;
                break;
            }
        }
        p.fixed_rooks = p.sides[0]._fr + p.sides[1]._fr;
    }
    if (fen[curr] != '-')
    {
        int enpassant_col = (BOARD_SIDE_LENGTH - 1) - (fen[curr] - a_ASCII_decimal);
        curr++;
        int enpassant_row;
        switch (fen[curr])
        {
        case '3': enpassant_row = 4; break;
        case '6':
            enpassant_row  = 5;
            p.side0isBlack = true;
            p.side0toMove  = blackToMove;
            uint64_t swp   = p.sides[0]._fr;
            p.sides[0]._fr = p.sides[1]._fr;
            p.sides[1]._fr = swp;
            break;
        default: printf("Unexpected enpassant row\n"); exit(1);
        }
        p.enpassant = enpassant_row * (BOARD_SIDE_LENGTH) + enpassant_col;
    }

    for (int i = 0; i < NUM_SQUARES; i++)
    {
        if (fen[i] == '1')
        {
            continue;
        }

        bool     white  = isupper(fen[i]);
        bool     side   = white == p.side0isBlack;
        uint64_t square = 1UL << ((NUM_SQUARES - 1) - i);
        char     lower  = tolower(fen[i]);
        switch (lower)
        {
        case 'b': p.sides[side].pieces[BISHOP] += square; break;
        case 'k': p.sides[side].pieces[KING] = square; break;
        case 'n': p.sides[side].pieces[KNIGHT] += square; break;
        case 'p': p.sides[side].pawns += square; break;
        case 'q': p.sides[side].pieces[QUEEN] += square; break;
        case 'r': p.sides[side].pieces[ROOK] += square; break;
        default: printf("Unexpected piece: %c\n", fen[i]); exit(1);
        }
    }

    if (p.enpassant && p.side0isBlack) // i.e. enpassant_row == 5
    {
        rotate_position_across_central_rows(&p);
    }
    return p;
}
