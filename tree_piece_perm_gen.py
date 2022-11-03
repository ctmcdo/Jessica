import itertools as it
import math
import numpy as np

# Can have 0, 1, or 2 rooks with castling rights
NUM_FIXED_ROOK_SCENARIOS = 3

# Dimension 0 represents the number of rooks with castling rights.
# Dimension 1 is the piece type, which varies for 2 fixed rooks because
# we purposedly keep the array values decreasing
NUM_NON_FIXED_CAPTURABLE_BASE_PIECES = [[2, 2, 2, 1], [2, 2, 1, 1], [2, 2, 1, 0]]
############################# b, n, r, q    b, n, r, q    b, n, q, r
############################################################### <-->
########################### Notice how we switch queens and rooks ^

# We define a covered set (for a fixed rook scenario) to be of the form [i, j, k, l],
# 2 >= i >= j >= k >= l >= 0, l < 2 and (a) k < 2 if we're considering a scenario
# with at least 1 fixed rook, (b) l == 0 if 2 fixed rooks, which corresponds to the
# base piece scenarios above. Every side has an associated covered set which represents
# its base pieces. We consider permutations of covered sets ignoring the constraint
# i >= j >= k >= l. Each permutation has an associated increase in minimum promotions.
# There are far less covered sets than sides and hence we drastically reduce the
# computation behind generating promotion-feasible permutations.
# The number of covered sets is sum_{i=1}^{4}(i + 1) = 14
NUM_COVERED_SETS = 14

# We can only start with at most 2 of any piece type (of the same colour), such as bishops
MAX_OF_ANY_ONE_BASE_PIECE = 2

# |{b, n, r, q}| == 4
NUM_PIECE_TYPES_LESS_KING = 4

# (2, 2, 1, 0) gives rise to permutation with greatest cost: (1, 0, 2, 2)
# which is base cost + 3. There are also base + 1, base + 2, ... so 4 in total
NUM_UNIQUE_PERM_COSTS = 4

CS_INDEX_SHAPE = (
    NUM_FIXED_ROOK_SCENARIOS,
    MAX_OF_ANY_ONE_BASE_PIECE + 1,
    MAX_OF_ANY_ONE_BASE_PIECE + 1,
    MAX_OF_ANY_ONE_BASE_PIECE + 1,
    MAX_OF_ANY_ONE_BASE_PIECE + 1,
)
PERMS_SHAPE = (
    NUM_FIXED_ROOK_SCENARIOS,
    NUM_COVERED_SETS,
    math.factorial(NUM_PIECE_TYPES_LESS_KING),
    NUM_PIECE_TYPES_LESS_KING,
)
ADDN_COST_TO_NPERMS_SHAPE = (
    NUM_FIXED_ROOK_SCENARIOS,
    NUM_COVERED_SETS,
    NUM_UNIQUE_PERM_COSTS,
)

# covered set indices -> for the fixed rook scenario and [i, j, k, l], get an
# associated index
fr_coveredSet_index = np.negative(np.ones(CS_INDEX_SHAPE, np.int32))

# for a fixed rook scenario, a covered set index, and a permutation number (which
# will be in [0, 24), return the permutation. The permutations are partially ordered
# by the minimum promotion cost
fr_coveredSetIndex_permIndex_perm = np.negative(np.ones(PERMS_SHAPE, np.int32))

# for a fixed rook scenario and a covered set index, get the number of permutations
# with minimum additional promotion cost of 0, 1, 2, 3. It's cumulative
fr_coveredSetIndex_permAddnCost_numPerms = np.negative(
    np.ones(ADDN_COST_TO_NPERMS_SHAPE, np.int32)
)


def covered_sets(fr_case_num):
    b = NUM_NON_FIXED_CAPTURABLE_BASE_PIECES[fr_case_num]
    s = b.copy()

    cset_num = 0
    fr_coveredSet_index[fr_case_num][s[0]][s[1]][s[2]][s[3]] = cset_num
    yield s, 0

    cset_num += 1
    while s[0] != 0:
        # generate next covered set
        for i in range(NUM_PIECE_TYPES_LESS_KING - 1, -1, -1):
            if s[i] != 0:
                s[i] -= 1
                for j in range(i + 1, NUM_PIECE_TYPES_LESS_KING):
                    s[j] = min(s[i], b[j])

                # and assign it an index
                fr_coveredSet_index[fr_case_num][s[0]][s[1]][s[2]][s[3]] = cset_num
                yield s, cset_num

                cset_num += 1
                break


def cost(s, b):
    c = 0
    for i in range(NUM_PIECE_TYPES_LESS_KING):
        diff = s[i] - b[i]
        if diff > 0:
            c += diff
    return c


def c_arr_literal_str_helper(a, depth):
    s = "{"
    if len(a.shape) == 1:
        for i in range(a.shape[0] - 1):
            s += str(a[i]) + ", "
        s += str(a[a.shape[0] - 1])
    else:
        for i in range(a.shape[0]):
            s += c_arr_literal_str_helper(a[i], depth + 1)
    s += "}"
    if depth != 0:
        s += ", "
    return s


def c_arr_literal_str(a):
    return c_arr_literal_str_helper(a, 0)


def c_dimensions_str(shape):
    s = ""
    for dim in shape:
        s += "[" + str(dim) + "]"
    return s


def format_c_arr_str(a, npname):
    s = "int "
    s += npname
    s += c_dimensions_str(np.shape(a))
    s += " = "
    s += c_arr_literal_str(a)
    s += ";"
    return s


for i in range(NUM_FIXED_ROOK_SCENARIOS):
    for (s, cset_num) in covered_sets(i):
        cost_to_perms = {}
        for p in it.permutations((0, 1, 2, 3)):
            permuted_set = [s[pe] for pe in p]
            c = cost(permuted_set, NUM_NON_FIXED_CAPTURABLE_BASE_PIECES[i])
            if c not in cost_to_perms:
                cost_to_perms[c] = []
            else:
                cost_to_perms[c].append(tuple(p))
        c = 0
        pi = 0
        for c, perms in sorted(cost_to_perms.items()):
            for e in perms:
                for j in range(NUM_PIECE_TYPES_LESS_KING):
                    fr_coveredSetIndex_permIndex_perm[i][cset_num][pi][j] = e[j]
                pi += 1

            fr_coveredSetIndex_permAddnCost_numPerms[i][cset_num][c] = pi

        to_be_copied = fr_coveredSetIndex_permAddnCost_numPerms[i][cset_num][c]
        for j in range(c + 1, NUM_UNIQUE_PERM_COSTS):
            fr_coveredSetIndex_permAddnCost_numPerms[i][cset_num][j] = to_be_copied


f = open("tree_piece_perm.c", "w")
f.write("// Generated by gen_prom_perms.py.\n")
f.write("// Formatted by clang thereafter\n\n")
f.write('#include "tree_common.h"\n\n')
f.write(format_c_arr_str(fr_coveredSet_index, "fr_coveredSet_index"))
f.write("\n")
f.write(
    format_c_arr_str(
        fr_coveredSetIndex_permIndex_perm, "fr_coveredSetIndex_permIndex_perm"
    )
)
f.write("\n")
f.write(
    format_c_arr_str(
        fr_coveredSetIndex_permAddnCost_numPerms,
        "fr_coveredSetIndex_permAddnCost_numPerms",
    )
)
