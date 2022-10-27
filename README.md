# Jessica
Computing an approximation of the number of chess positions.

I'm leaving this project unfinished. The filters are probably buggy and fen printing may need attention.
However, I think it's impressive that we build the search tree so quickly (in 1/4 of a second on my 4-core laptop), that's where the majority of the effort has gone.
The sample space is approximately 10^47 large. The computation is so quick for a number of reasons. Most of the reasons should be clear from comments in tree_create.h but the permutations logic is mostly undocumented. You'll have to figure it out for yourself by reading the python 😊

# Dependencies
GMP
