#ifndef COOS_H_INCLUDED
#define COOS_H_INCLUDED

#include "utility.h"

extern void matr_mult_coos(const void* a, const void* b, void* result);

void matr_mult_coos_c(const coos_matrix* a, const coos_matrix* b, coos_matrix* res);
void matr_mult_coos_map_c(const map* a, const map* b, coos_matrix* result);

void benchmark_coos(coos_matrix* c1, coos_matrix* c2, uint64_t iterations);


#endif
