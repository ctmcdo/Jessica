EXECUTABLE_NAME = Jessica

CC = gcc
LIBS = -lgmp -lm -L filter_pawn/build/lib -l:libfilter_pawn.so -Wl,-rpath filter_pawn/build/lib
CFLAGS = -pthread -mabm -mbmi -mbmi2 -Wall
OPTIMISATION = -O3

INCLUDE = -I filter_pawn
SRC = $(filter-out $(wildcard *test.c unused_check_filter.c), $(wildcard *.c))

all:
	$(CC) ${INCLUDE} ${SRC} -o ${EXECUTABLE_NAME} $^ ${LIBS} ${CFLAGS} ${OPTIMISATION} -DNDEBUG

debug:
	$(CC) ${INCLUDE} ${SRC} -o ${EXECUTABLE_NAME} $^ ${LIBS} ${CFLAGS} -g

format:
	find . -name '*.h' -o -name '*.c' -o -name '.*.cc' | xargs clang-format -i && black tree_piece_perm_gen.py
