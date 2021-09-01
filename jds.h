#ifndef JDS_H_INCLUDED
#define JDS_H_INCLUDED

#include "utility.h"

extern void matr_mult_jds(const void* a, const void* b, void* result);
void matr_mult_jds_c(jds_matrix* a, jds_matrix* b, jds_matrix* result);

void benchmark_jds(jds_matrix* j1, jds_matrix* j2, uint64_t iterations);
void coos_to_jds_c(coos_matrix* coos, uint64_t* nnz_per_row, jds_matrix* jds);

#endif
