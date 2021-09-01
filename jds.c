#include "jds.h"

/*
 * Multiplies two JDS matrices and stores the result in result
 */
void matr_mult_jds_c(jds_matrix* a, jds_matrix* b, jds_matrix* result)
{
    /*
     * create COOS intermediate result
     */
    coos_matrix* coos = malloc(sizeof(coos_matrix));
    coos->num_entries = 0;
    coos->num_rows= a->num_rows;
    coos->num_cols = b->num_cols;
    if(a->length_values * b->length_values > a->num_rows * b->num_cols)
    {
        coos->entries = calloc((a->num_rows * b->num_cols), sizeof(coos_triple));
    }
    else
    {
        coos->entries = calloc((a->length_values * b->length_values), sizeof(coos_triple));
    }
    /*
     * Stores the number of non-zero elements stored at each row index
     */
    uint64_t* nnz_per_row = calloc(a->num_rows, sizeof(uint64_t));

    /*
     * loop through each jagged diagonal of the matrix b
     */
    for(uint64_t jd = 0; jd < b->length_col_ptr-1; ++jd)
    {
        uint64_t perm_iterator_b = 0;
        /*
         * for each element of the jagged diagonal iterate through the first matrix
         */
        for(uint64_t element=b->col_ptr[jd]; element < b->col_ptr[jd+1]; ++element)
        {
            //current element only affects its column in the result
            uint64_t col = b->col_idx[element];
            //loop through jagged diagonals of first matrix
            for(uint64_t j=0; j < a->length_col_ptr-1; ++j)
            {
                //loop through non-zero elements of each jagged diagonal of first matrix
                uint64_t permu_iterator_a = 0;
                for(uint64_t k=a->col_ptr[j]; k < a->col_ptr[j+1]; ++k){
                    //To multiply 2 elements p of matrix A with q of matrix B, the original row index
                    //of q must be equal to the original column index of p
                    if(a->col_idx[k] == b->permutation[perm_iterator_b])
                    {
                        float product = a->values[k] * b->values[element];
                        uint64_t row = a->permutation[permu_iterator_a];
                        uint64_t nnz_before = coos->num_entries;
                        store_coos_c(row, col, product, coos);

                        if(coos->num_entries != nnz_before)
                        {
                            nnz_per_row[row] += 1;
                        }
                    }
                    permu_iterator_a += 1;
                }
            }
            perm_iterator_b += 1;
        }
    }
    /*
     * This function eliminates all the COOS tripels that contain a value of 0. This could happen when
     * equally valued positive and negative values are added to the same row and column. If this function
     * is unseccussful, it sets the number of rows of the COOS paramter matrix to 0.
     */
    eliminate_zeros_c(coos);
    /*
     * If the above function is successfully executed, we can proceed to transform the COOS intermediate
     * matrix into the JDS format.
     */
    if(coos->num_rows != 0)
    {
        coos_to_jds_c(coos, nnz_per_row, result);
    }
    free_coos(coos);
    free(nnz_per_row);
}

/*
 * multiplies the JDS matrices j1 and jc2, n number of times and measures the time it takes.
 * This is done for both implementations of JDS.
 * The resulting time and average time are printed to the console for each implementation
 */
void benchmark_jds(jds_matrix* j1, jds_matrix* j2, uint64_t iterations)
{
    if(j1->num_rows != j2->num_cols)
    {
        perror("No benchmarking done. Invalid multiplication");
        return;
    }
    if(j1->length_values == 0 || j2->length_values == 0)
    {
        perror("No benchmarking done. Input matrix is empty");
        return;
    }



    struct timespec start;
    struct timespec end;
    double time = 0.0;
    double avg_time = 0.0;

    printf("\n-Benchmark JDS matrix multiplication(for %ld Iterations):\n", iterations);

    for (uint64_t i = 0; i < iterations; ++i) {
        jds_matrix* result = initialize_jds_result(j1,j2);
        if(result == NULL) return;
        clock_gettime(CLOCK_MONOTONIC, &start);
        matr_mult_jds_c(j1, j2, result);
        clock_gettime(CLOCK_MONOTONIC, &end);
        time += end.tv_sec - start.tv_sec + 1e-9 * (end.tv_nsec - start.tv_nsec);
        free_jds(result);
    }

    avg_time = (time / iterations);
    printf("\tC Implementation:\n\t\ttotal: \t\t%.8f seconds\n\t\taverage: \t%.8f seconds\n", time, avg_time);

    time = 0.0;
    for (uint64_t i = 0; i < iterations; ++i) {
        jds_matrix* result = initialize_jds_result(j1,j2);
        if(result == NULL) return;
        clock_gettime(CLOCK_MONOTONIC, &start);
        matr_mult_jds(j1, j2, result);
        clock_gettime(CLOCK_MONOTONIC, &end);
        time += end.tv_sec - start.tv_sec + 1e-9 * (end.tv_nsec - start.tv_nsec);
        free_jds(result);
    }
    avg_time = (time / iterations);
    printf("\tAssembly Implementation:\n\t\ttotal: \t\t%.8f seconds\n\t\taverage: \t%.8f seconds\n\n", time, avg_time);
}

/*
 * Takes a COOS matrix and the number of non-zero elements at each row index
 * and then outputs the same matrix as a JDS matrix
 */
void coos_to_jds_c(coos_matrix* coos, uint64_t* nnz_per_row, jds_matrix* jds)
{
    /*
     * have an array containing 0, 1, ..., num_rows-1
     * sort this array with respect to nnz_per_row
     */
    uint64_t num_rows = coos->num_rows;

    sort_c(nnz_per_row, num_rows, jds->permutation);
    sort_triples_col_c(coos);

    uint64_t num_non_zeros = coos->num_entries;
    jds->length_values = num_non_zeros;
    uint64_t num_of_JD = nnz_per_row[0];

    jds->length_col_ptr = num_of_JD + 1;
    jds->col_ptr[0] = 0;

    uint64_t values_iterator = 0;
    /*
     * Iterate through each Jagged diagonal to assign the values in it
     */
    for(uint64_t i = 0; i < num_of_JD; ++i)
    {
        /*
         * Iterate through each row of the matrix that has non-zero elements in it
         */
        for(uint64_t j = 0; j < num_rows && nnz_per_row[j] != 0; ++j)
        {
            uint64_t skip = 0;
            /*
             * Iterate through each entry in the COOS matrix.
             */
            for(uint64_t k = 0; k < num_non_zeros; ++k)
            {
                if(coos->entries[k].row_idx == jds->permutation[j])
                {
                    if(skip == i)
                    {
                        jds->values[values_iterator] = coos->entries[k].val;
                        jds->col_idx[values_iterator] = coos->entries[k].col_idx;
                        values_iterator++;
                        break;
                    }
                    skip += 1;
                }
            }

            nnz_per_row[j] -= 1;
        }
        jds->col_ptr[i+1] = values_iterator;
    }
}

