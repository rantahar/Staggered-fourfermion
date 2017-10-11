CC=gcc 
CFLAGS=-march=native -Wall -Wextra -std=c99 -O3 -llapack -lm #-DDEBUG

DEPS=Makefile Staggered.h mersenne.h

default: Staggered exact_test

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

Staggered: Staggered.o mersenne_inline.o fermion_matrix.o worm_update.o local_update.o invert_cg.o $(DEPS)
	$(CC) -o Staggered Staggered.o mersenne_inline.o fermion_matrix.o worm_update.o local_update.o invert_cg.o  $(CFLAGS)

exact_test: exact_test.o mersenne_inline.o fermion_matrix.o worm_update.o local_update.o $(DEPS)
	$(CC) -o exact_test exact_test.o mersenne_inline.o fermion_matrix.o worm_update.o local_update.o  $(CFLAGS)


clean:
	rm -f *.o  

