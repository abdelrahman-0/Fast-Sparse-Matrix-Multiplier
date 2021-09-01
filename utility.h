#ifndef UTILITY_H_INCLUDED
#define UTILITY_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>

typedef struct coos_triple
{
    uint64_t row_idx;
    uint64_t col_idx;
    float val;
} coos_triple;

typedef struct coos_matrix
{
    uint64_t num_rows;
    uint64_t num_cols;
    uint64_t num_entries;
    coos_triple* entries;
} coos_matrix;

typedef struct jds_matrix{
    uint64_t num_rows;
    uint64_t num_cols;
    uint64_t length_values;
    float* values;
    uint64_t* col_idx;
    uint64_t length_col_ptr;
    uint64_t* col_ptr;
    uint64_t* permutation;
} jds_matrix;

typedef struct pair
{
    uint64_t first;
    float second;
}pair;
typedef struct map_entry
{
    uint64_t idx;
    uint64_t num_elements;
    pair* pairs;
}map_entry;
typedef struct map
{
    uint64_t num_elements;
    map_entry* array;
}map;

extern void store_coos(uint64_t row, uint64_t col, float product, coos_matrix* c);
extern void eliminate_zeros(coos_matrix* coos);

void store_coos_c(uint64_t row, uint64_t col, float value, coos_matrix* coos);
void eliminate_zeros_c(coos_matrix* coos);

void sort_triples_col_c(coos_matrix* coos);
void sort_triples_row_c(coos_matrix* coos);
void swap_triples(coos_triple* array, uint64_t i, uint64_t j);

void sort_c(uint64_t* array_1, uint64_t size, uint64_t* array_2);
void swap(uint64_t* array, uint64_t i, uint64_t j);

void print_help();
int check_malloc(void* ptr);

float** coos_to_2d(coos_matrix* coos);
float** matr_mul_2d(coos_matrix* a, coos_matrix* b, float** f1, float** f2);

coos_matrix* read_input_coos(char* file_name);
jds_matrix* read_input_jds(char* file_name);

void coos_output(coos_matrix* matr);
void jds_output(jds_matrix* matr);

coos_matrix* initialize_coos_result(coos_matrix* c1, coos_matrix* c2);
jds_matrix* initialize_jds_result(jds_matrix* j1,jds_matrix* j2);

map* to_map_col(coos_matrix* in);
map* to_map_row(coos_matrix* in);

void free_coos(coos_matrix* c1);
void free_map(map* map1);
void free_jds(jds_matrix* j1);
void free_2d_matr(float** f1, uint64_t num_rows);

#endif
