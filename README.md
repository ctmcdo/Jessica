# Jessica
Generate random chess positions quickly.

We generate positions from a sample space approximately 10^47 large. The criterion for whether a board is included in our space is defined in prom_slack.c.

We use threading to build the tree, but not for sampling. To generate a large number of samples it would be useful to spawn multiple threads in main and sample within each thread. On a 4-core laptop with hyper-threading, I can build and sample 10,000 positions in .4 seconds. The reason for being so quick is documented in tree_create.c. It's mostly due to:
1. Not including duplicate chess sides where the only difference is piece type. For instance, two sides which are the same except for one side's number of bishops is the other side's number of knights and vice versa shouldn't be included twice in the tree. Instead we generalise the piece type, and later assign the generic pieces actual piece types, i.e, rook, bishop, knight, queen. The caveat is rooks with castling rights complicate the matter somewhat.
2. Noticing that there are only a limited number of permutations of piece types for a generic side. Hence we don't need to constantly compute the number of permutations. To understand this better please take a look at tree_piece_perm_gen.py.
3. Writing the program in c.

# Dependencies
An x86 processor which supports -mabm -mbmi -mbmi2
GMP
