EXECUTABLE_NAME = Jess

CC = gcc
LIBS = -lgmp -lm 
CFLAGS = -fPIC -pthread -mabm -mbmi -mbmi2
OPTIMISATION = -O3

#SRC = $(filter-out $(wildcard *test.c), $(wildcard *.c))
#SRC += -I dependencies/hungarian dependencies/hungarian/hungarian.c

SRC = main.c prom_slack.c tree_piece_perm.c util.c -I dependencies/hungarian dependencies/hungarian/hungarian.c

all:
	$(CC) ${SRC} -o ${EXECUTABLE_NAME} $^ ${LIBS} ${CFLAGS} ${OPTIMISATION}

debug:
	$(CC) ${SRC} -o ${EXECUTABLE_NAME} $^ ${LIBS} ${CFLAGS} -g

format:
	find . -name '*.h' ! -path './dependencies/*' \
		-o -name '*.c' ! -path './dependencies/*' | xargs clang-format -i \
		&& black tree_piece_perm_gen.py
