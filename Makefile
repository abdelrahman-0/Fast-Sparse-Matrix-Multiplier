# Optimize, turn on additional warnings
CFLAGS=-O3 -g -Wall -Wextra -no-pie

.PHONY: all
all: main
main: main.c coos.c jds.c utility.c matr_mult_coos.S matr_mult_jds.S store_coos.S eliminate_zeros.S coos_to_jds.S sort_triples_col.S sort.S
	$(CC) $(CFLAGS) -o $@ $^
.PHONY: clean
clean:
	rm -f main
