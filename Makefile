EXECUTABLE_NAME = Jess

CC = gcc
LIBS = -lgmp -lm -L filter_pawn/build/lib -l:libfilter_pawn -Wl,-rpath filter_pawn/build/lib
CFLAGS = -pthread -mabm -mbmi -mbmi2 -Wall
OPTIMISATION = -O3

SRC = $(filter-out $(wildcard *test.c unused_check_filter.c), $(wildcard *.c))

all:
	$(CC) ${SRC} -o ${EXECUTABLE_NAME} $^ ${LIBS} ${CFLAGS} ${OPTIMISATION} -DNDEBUG

debug:
	$(CC) ${SRC} -o ${EXECUTABLE_NAME} $^ ${LIBS} ${CFLAGS} -g

format:
	find . -name '*.h' ! -path './dependencies/*' \
		-o -name '*.c' ! -path './dependencies/*' | xargs clang-format -i \
		&& black tree_piece_perm_gen.py
