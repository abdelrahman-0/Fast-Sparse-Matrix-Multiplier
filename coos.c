#include "coos.h"

/*
 * Multiplies two COOS matrices and stores the result in res
 */
void matr_mult_coos_c(const coos_matrix* a, const coos_matrix* b, coos_matrix* res)
{
    /*
     * loop through non-zero elements of matrix A
     */
    for(uint64_t i = 0; i < a->num_entries; ++i)
    {
        uint64_t col_idx = a->entries[i].col_idx;
        /*
         * loop through non-zero elements of matrix B
         */
        for(uint64_t j = 0; j < b->num_entries; ++j)
        {
            /*
             * multiply the two non-zero elements ony if the column index of the first element
             * is equal to the column index of the second element
             */
            if(col_idx == b->entries[j].row_idx)
            {
                /*
                 * if the condition is passed,
                 * the elements are multiplied and the result is added to the correct position in the result matrix
                 */
                float product = a->entries[i].val * b->entries[j].val;
                /*
                 * function to store
                 */
                store_coos_c(a->entries[i].row_idx, b->entries[j].col_idx, product, res);
            }
        }
    }
    /*
     * Remove all entires with the value 0
     */
    eliminate_zeros_c(res);
}

/*
 * Multiplies two COOS matrices stored in a Map
 */
void matr_mult_coos_map_c(const map* a, const map* b, coos_matrix* result)
{
    int found = 0;
    uint64_t start_here = 0;
    /*
     * loop through cols of A
     */
    for(uint64_t j=0; j < a->num_elements; ++j)
    {
        uint64_t current_col = a->array[j].idx;
        /*
         * check if container B has row equal to current column of A (if that's not the case, we can simply skip
         * over to the next column in A, keeping track of where we stopped in the rows of B (since the containers
         * of the cols of A and the rows of B are sorted))
         */
        for(uint64_t s=start_here; s < b->num_elements; ++s)
        {
            /*
             * no need to continue looking for row in B if row > current_col since rows are sorted in B
             */
            if(current_col < b->array[s].idx)
            {
                break;
            }
            /*
             * check if B contains row equal to current column of A
             */
            if(current_col == b->array[s].idx)
            {
                found = 1;
            }
            /*
             * update where to start looking in the next iteration (since both containers of A and B are sorted)
             */
            start_here++;
        }
        if(found)
        {
            /*
             * loop through rows of A for a fixed column (current_col)
             */
            for(uint64_t i=0; i < a->array[j].num_elements; ++i)
            {
                uint64_t row_A = a->array[j].pairs[i].first;
                float value_A = a->array[j].pairs[i].second;
                /*
                 * loop through the columns of B for a fixed row
                 */
                for(uint64_t k=0; k < b->array[start_here-1].num_elements; ++k)
                {
                    /*
                     * multiply the 2 floats and add them to results
                     */
                    uint64_t col_B = b->array[start_here-1].pairs[k].first;
                    float value_B = b->array[start_here-1].pairs[k].second;
                    float product = value_A * value_B;
                    store_coos_c(row_A,col_B,product,result);
                }
            }
        }
        found = 0;
    }
    /*
     * Remove all entires with the value 0
     */
    eliminate_zeros_c(result);
}

/*
 * multiplies the COOS matrices c1 and c2, n number of times and measures the time it takes.
 * This is done for all implementations of COOS, as well as multiplying them as 2D float-arrays.
 * The resulting time and average time are printed to the console for each implementation
 */
void benchmark_coos(coos_matrix* c1, coos_matrix* c2, uint64_t iterations)
{
    if(c1->num_rows != c2->num_cols)
    {
        perror("No benchmarking done. Invalid multiplication");
        return;
    }
    if(c1->num_entries == 0 || c2->num_entries == 0)
    {
        perror("No benchmarking done. Input matrix is empty");
        return;
    }

    struct timespec start;
    struct timespec end;
    double time = 0.0;
    double avg_time = 0.0;

    printf("\n-Benchmark COOS matrix multiplication(for %ld Iterations):\n", iterations);

    for (uint64_t i = 0; i < iterations; ++i)
    {
        coos_matrix *res = initialize_coos_result(c1,c2);
        if(res == NULL) return;
        clock_gettime(CLOCK_MONOTONIC, &start);
        matr_mult_coos_c(c1, c2, res);
        clock_gettime(CLOCK_MONOTONIC, &end);
        time += end.tv_sec - start.tv_sec + 1e-9 * (end.tv_nsec - start.tv_nsec);
        free_coos(res);
    }

    avg_time = (time / iterations);
    printf("\tC Implementation:\n\t\ttotal: \t\t%.8f seconds\n\t\taverage: \t%.8f seconds\n", time, avg_time);

    map *map1 = to_map_col(c1);
    map *map2 = to_map_row(c2);
    if(map1 == NULL || map2 == NULL) return;

    time = 0.0;
    for (uint64_t i = 0; i < iterations; ++i)
    {
        coos_matrix *res = initialize_coos_result(c1,c2);
        if(res == NULL) return;
        clock_gettime(CLOCK_MONOTONIC, &start);
        matr_mult_coos_map_c(map1, map1, res);
        clock_gettime(CLOCK_MONOTONIC, &end);
        time += end.tv_sec - start.tv_sec + 1e-9 * (end.tv_nsec - start.tv_nsec);
        free_coos(res);
    }
    avg_time = (time / iterations);
    printf("\tMap C Implementation:\n\t\ttotal: \t\t%.8f seconds\n\t\taverage: \t%.8f seconds\n", time, avg_time);

    free_map(map1);
    free_map(map2);

    time = 0.0;
    for (uint64_t i = 0; i < iterations; ++i)
    {
        coos_matrix *res = initialize_coos_result(c1,c2);
        if(res == NULL) return;
        clock_gettime(CLOCK_MONOTONIC, &start);
        matr_mult_coos(c1, c2, res);
        clock_gettime(CLOCK_MONOTONIC, &end);
        time += end.tv_sec - start.tv_sec + 1e-9 * (end.tv_nsec - start.tv_nsec);
        free_coos(res);
    }
    avg_time = (time / iterations);
    printf("\tAssembly Implementation:\n\t\ttotal: \t\t%.8f seconds\n\t\taverage: \t%.8f seconds\n", time, avg_time);

    float** f1 = coos_to_2d(c1);
    float** f2 = coos_to_2d(c2);
    if(f1 == NULL || f2 == NULL) return;
    float** matr_2d = NULL;

    time = 0.0;
    for (uint64_t i = 0; i < iterations; ++i)
    {
        clock_gettime(CLOCK_MONOTONIC, &start);
        matr_2d = matr_mul_2d(c1,c2,f1,f2);
        clock_gettime(CLOCK_MONOTONIC, &end);
        time += end.tv_sec - start.tv_sec + 1e-9 * (end.tv_nsec - start.tv_nsec);
        if(matr_2d == NULL) return;
        free_2d_matr(matr_2d, c1->num_rows);
    }

    avg_time = (time / iterations);
    printf("\t2D Matrix Multiplication  Implementation:\n\t\ttotal: \t\t%.8f seconds\n\t\taverage: \t%.8f seconds\n", time, avg_time);

    free_2d_matr(f1, c1->num_rows);
    free_2d_matr(f2, c2->num_rows);
}
