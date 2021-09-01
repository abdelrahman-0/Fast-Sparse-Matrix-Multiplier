#include "utility.h"

/*
 * Stores a new coos_triple with the given parameters.
 * If a triple already exists with the given coordinates,
 * the given value is added to the one in the matrix
 */
void store_coos_c(uint64_t row, uint64_t col, float value, coos_matrix* coos)
{
    uint64_t i = 0;
    for(; i < coos->num_entries; ++i)
    {
        if(coos->entries[i].row_idx == row && coos->entries[i].col_idx == col)
        {
            coos->entries[i].val += value;
            return;
        }
    }
    coos_triple triple = {.row_idx = row, .col_idx = col, .val = value};
    coos->entries[i] = triple;
    coos->num_entries++;
}

/*
 * Goes through the given matrix and removes all elements with the value of 0
 */
void eliminate_zeros_c(coos_matrix* coos)
{
    coos_triple* new_entries = malloc(coos->num_entries * sizeof(coos_triple));
    if(check_malloc(new_entries))
    {
        coos->num_rows = 0;
        return;
    }
    uint64_t new_entries_iterator = 0;
    for(uint64_t i = 0; i < coos->num_entries; ++i)
    {
        if(coos->entries[i].val != 0.0)
        {
            new_entries[new_entries_iterator] = coos->entries[i];
            coos->num_entries++;
            new_entries_iterator++;
        }
        coos->num_entries--;
    }
    free(coos->entries);
    coos->entries =  new_entries;
}

/*
 * Sorts the COOS matrix triples by column
 */
void sort_triples_col_c(coos_matrix* coos)
{
    coos_triple* entries = coos->entries;
    for(uint64_t i = 1; i < coos->num_entries; ++i)
    {
        for(uint64_t j = 0; j < i; ++j)
        {
            if(entries[i].col_idx < entries[j].col_idx)
            {
                swap_triples(entries, i, j);
            }
        }
    }
}
/*
 * Sorts the COOS matrix triples by row
 */
void sort_triples_row_c(coos_matrix* coos)
{
    coos_triple* entries = coos->entries;
    for(uint64_t i = 1; i < coos->num_entries; ++i)
    {
        for(uint64_t j = 0; j < i; ++j)
        {
            if(entries[i].row_idx < entries[j].row_idx)
            {
                swap_triples(entries, i, j);
            }
        }
    }

}
/*
 * Swaps two COOS triples
 */
void swap_triples(coos_triple* array, uint64_t i, uint64_t j)
{
    coos_triple temp = array[i];
    array[i] = array[j];
    array[j] = temp;
}

/*
 * Takes two arrays and sorts both in respect to the values of the first one.
 * This is done in descending order.
 */
void sort_c(uint64_t* array_1, uint64_t size, uint64_t* array_2)
{
    for(uint64_t i = 1; i < size; ++i)
    {
        for(uint64_t j = 0; j < i; ++j)
        {
            if(array_1[i] > array_1[j])
            {
                swap(array_1, i, j);
                swap(array_2, i, j);
            }
        }
    }
}
void swap(uint64_t* array, uint64_t i, uint64_t j)
{
    uint64_t temp = array[i];
    array[i] = array[j];
    array[j] = temp;
}

/*
 * Prints the guide for the program arguments (-h/--help).
 */
void print_help()
{
    printf("Please use any combination of the follwing flags:\n");

    printf("-h (--help)\n");
    printf("\tShows available program arguments\n\n");

    printf("-c file1.txt file2.txt\n");
    printf("\tThis is the COOS flag. It takes two input .txt files, each containing the input matrices in the COOS format. An example of a 5x5 identity matrix is found in coos_example_input.txt.\n");
    printf("\tThe files should have the following structure:\n");
    printf("\t\t(num_rows) (num_cols)\n");
    printf("\t\t(num_of_nonzero_elements)\n");
    printf("\t\t[list of space-separated float values]\n");
    printf("\t\t[list of space-separated column indexes]\n");
    printf("\t\t[list of space-separated row indexes]\n");
    printf("\n\tThe output is found in coos_output.txt\n\n");

    printf("-j file1.txt file2.txt\n");
    printf("\tThis is the JDS flag. It takes two input .txt files, each containing the input matrices in the JDS format. An example of a 5x5 identity matrix is found in jds_example_input.txt\n");
    printf("\tThe files should have the following structure:\n");
    printf("\t\t(num_rows) (num_cols)\n");
    printf("\t\t(num_of_nonzero_elements)\n");
    printf("\t\t[list of space-separated float values]\n");
    printf("\t\t[list of space-separated column indexes]\n");
    printf("\t\t(length_col_ptr)\n");
    printf("\t\t[list of space-separated column pointer indexes]\n");
    printf("\t\t[list of space-separated row permutations]\n");

    printf("\n\tThe output is found in jds_output.txt\n\n");

    printf("-b NUM\n");
    printf("\tThis is the Benchmarking flag. The multiplication of the two matrices is done for NUM many iterations. The total and average running times of the iterations are then output to the console\n\n");
}
/*
 *  Checks if pointer is NULL.
 *  If it is, it outputs an error.
 */
int check_malloc(void* ptr)
{
    if(ptr == NULL)
    {
        perror("Unable to allocate enough memory using malloc\n");
        return 1;
    }
    return 0;
}

/*
 * Converts a COOS matrix to a matrix stored in a 2D float-array.
 */
float** coos_to_2d(coos_matrix* coos)
{
    float** result = malloc(coos->num_rows * 8);
    if(check_malloc(result))return NULL;
    for (uint64_t j = 0; j < coos->num_rows; ++j) {
        result[j] = calloc(coos->num_cols, sizeof(float));
        if(check_malloc(result[j]))return NULL;
    }

    for (uint64_t i = 0; i < coos->num_entries; ++i) {
        result[coos->entries[i].row_idx][coos->entries[i].col_idx] = coos->entries[i].val;
    }
    return result;
}
/*
 * Multiplies two 2D float-array matrices and returns the value as a 2D float-array.
 */
float** matr_mul_2d(coos_matrix* a, coos_matrix* b, float** f1, float** f2)
{
    float** result = malloc(a->num_rows  * 8);
    if(check_malloc(result))return NULL;
    for (uint64_t j = 0; j < a->num_rows; ++j) {
        result[j] = calloc(b->num_cols, sizeof(float));
        if(check_malloc(result[j]))return NULL;
    }

    for(uint64_t i=0;i<a->num_rows;i++)
    {
        for(uint64_t j=0;j<b->num_cols;j++)
        {
            result[i][j]=0;
            for(uint64_t k=0;k<a->num_cols;k++)
            {
                result[i][j]+=f1[i][k]*f2[k][j];
            }
        }
    }
    return result;
}

/*
 * Reads a file in the COOS input format and returns it as a coos_matrix.
 */
coos_matrix* read_input_coos(char* file_name)
{
    FILE* file = fopen(file_name, "r");
    //exit program if file cannot be opened
    if (file == NULL)
    {
        perror("Cannot open file");
        return NULL;
    }

    uint64_t nnz, row, col;
    char c;
    if(fscanf(file, "%ld", &row )!= 1) {perror("Invalid input file format"); return NULL;}
    if(fscanf(file, "%ld", &col )!= 1) {perror("Invalid input file format"); return NULL;}
    if(fscanf(file, "%ld", &nnz )!= 1) {perror("Invalid input file format"); return NULL;}
    coos_matrix* coos = malloc(sizeof(coos_matrix));
    if(check_malloc(coos))return NULL;

    coos->num_entries = nnz;
    coos->num_cols = col;
    coos->num_rows = row;
    if(nnz == 0)return coos;

    float* values = calloc(nnz, sizeof(float));
    if(check_malloc(values))return NULL;

    uint64_t* rows = calloc(nnz, sizeof(uint64_t));
    if(check_malloc(rows))return NULL;

    uint64_t* cols = calloc(nnz, sizeof(uint64_t));
    if(check_malloc(cols))return NULL;

    //move to the next line
    while ((c = fgetc(file)) != EOF && (c!= '\n')){}
    for (uint64_t i = 0; i < nnz; i++)
    {
        if(fscanf(file, "%f", &(values[i]))!= 1) {perror("Invalid input file format"); return NULL;}
    }
    //move to the next line
    while ((c = fgetc(file)) != EOF && (c!= '\n')){}
    for (uint64_t i = 0; i < nnz; i++)
    {
        if(fscanf(file, "%ld", &(cols[i]))!= 1) {perror("Invalid input file format"); return NULL;}
        if(cols[i] >= coos->num_cols) {perror("Error in COOS Input: Row index Out of Bounds"); return NULL;}
    }
    //move to the next line
    while ((c = fgetc(file)) != EOF && (c!= '\n')){}
    for (uint64_t i = 0; i < nnz; i++)
    {
        if(fscanf(file, "%ld", &(rows[i]))!= 1) {perror("Invalid input file format"); return NULL;}
        if(rows[i] >= coos->num_rows) {perror("Error in COOS Input: Column index Out of Bounds"); return NULL;}
    }

    coos->entries = calloc(nnz, sizeof(coos_triple));
    if(check_malloc(coos->entries))return NULL;

    for (uint64_t i = 0; i < nnz; ++i)
    {
        coos_triple entry;
        entry.val = values[i];
        entry.row_idx = rows[i];
        entry.col_idx = cols[i];
        coos->entries[i] = entry;
    }
    free(values);
    free(rows);
    free(cols);
    fclose(file);
    return coos;
}
/*
 * Reads a file in the JDS input format and returns it as a jds_matrix.
 */
jds_matrix* read_input_jds(char* file_name)
{
    FILE* file = fopen(file_name, "r");
    //exit program if file cannot be opened
    if (file == NULL)
    {
        perror("Cannot open file");
        return NULL;
    }

    uint64_t NNZ, row, col, length_col_ptr;
    char c;
    if(fscanf(file, "%ld", &row )!= 1) {perror("Invalid input file format"); return NULL;}
    if(fscanf(file, "%ld", &col )!= 1) {perror("Invalid input file format"); return NULL;}
    if(fscanf(file, "%ld", &NNZ )!= 1) {perror("Invalid input file format"); return NULL;}

    if(row == 0 || col == 0){perror("Number of rows and columns in input matrix needs to be >0"); return NULL;}

    //initialize jds
    jds_matrix *jds = malloc(sizeof(jds_matrix));
    if(check_malloc(jds))return NULL;

    jds->num_cols = col;
    jds->num_rows = row;
    jds->length_values = NNZ;
    if(NNZ == 0)return jds;

    jds->values = malloc(NNZ * sizeof(float));
    if(check_malloc(jds->values))return NULL;

    //move to the next line
    while ((c = fgetc(file)) != EOF && (c!= '\n')){}
    for (uint64_t i = 0; i < NNZ; i++)
    {
        if(fscanf(file, "%f", &(jds->values[i]))!= 1) {perror("Invalid input file format"); return NULL;}
    }

    jds->col_idx = malloc(NNZ * sizeof(uint64_t));
    if(check_malloc(jds->col_idx))return NULL;

    //move to the next line
    while ((c = fgetc(file)) != EOF && (c!= '\n')){}
    for (uint64_t i = 0; i < NNZ; i++)
    {
        if(fscanf(file, "%ld", &(jds->col_idx[i]))!= 1) {perror("Invalid input file format"); return NULL;}
        if(jds->col_idx[i] >= jds->num_cols) {perror("Error in JDS Input: Column index Out of Bounds"); return NULL;}
    }

    while ((c = fgetc(file)) != EOF && (c!= '\n')){}
    if(fscanf(file, "%ld", &length_col_ptr )!= 1) {perror("Invalid input file format"); return NULL;}

    jds->length_col_ptr = length_col_ptr;
    jds->col_ptr = malloc(length_col_ptr * sizeof(uint64_t));
    if(check_malloc(jds->col_ptr))return NULL;

    //move to the next line
    while ((c = fgetc(file)) != EOF && (c!= '\n')){}
    for (uint64_t i = 0; i < length_col_ptr; i++)
    {
        if(fscanf(file, "%ld", &(jds->col_ptr[i]))!= 1) {perror("Invalid input file format"); return NULL;}
        if(jds->col_ptr[i] > jds->length_values) {perror("Error in JDS Input: Column pointer Out of Bounds"); return NULL;}
    }

    jds->permutation  = malloc(row * sizeof(uint64_t));
    if(check_malloc(jds->permutation))return NULL;

    //move to the next line
    while ((c = fgetc(file)) != EOF && (c!= '\n')){}
    for (uint64_t i = 0; i < row; i++)
    {
        if(fscanf(file, "%ld ", &(jds->permutation[i]))!= 1) {perror("Invalid input file format"); return NULL;}
        if(jds->permutation[i] >= jds->num_rows) {perror("Error in JDS Input: Permutation Out of Bounds"); return NULL;}
    }
    while ((c = fgetc(file)) != EOF){}

    fclose(file);
    return jds;
}

/*
 * Takes a coos_matrix and stores in coos_output.txt
 * The file is written in the same format as the COOS input
 */
void coos_output(coos_matrix* matr)
{
    FILE* file = fopen("coos_output.txt", "w");
    if (file == NULL)
    {
        perror("Unable to create output file");
        exit(EXIT_FAILURE);
    }

    fprintf(file, "%ld ", matr->num_rows);
    fprintf(file, "%ld \n", matr->num_cols);
    fprintf(file, "%ld\n", matr->num_entries);
    if(matr->num_entries == 0)
    {
        fclose(file);
        return;
    }
    for (uint64_t i = 0; i < matr->num_entries; ++i)
    {
        fprintf(file, "%f ", matr->entries[i].val);
    }
    fprintf(file, "\n");
    for (uint64_t i = 0; i < matr->num_entries; ++i)
    {
        fprintf(file, "%ld ", matr->entries[i].row_idx);
    }
    fprintf(file, "\n");
    for (uint64_t i = 0; i < matr->num_entries; ++i)
    {
        fprintf(file, "%ld ", matr->entries[i].col_idx);
    }
    fprintf(file, "\n");
    fclose(file);
}
/*
 * Takes a jds_matrix and stores in jds_output.txt
 * The file is written in the same format as the JDS input
 */
void jds_output(jds_matrix* matr)
{
    FILE* file = fopen("jds_output.txt", "w");

    if (file == NULL)
    {
        perror("Unable to create output file");
        exit(EXIT_FAILURE);
    }

    fprintf(file, "%ld ", matr->num_rows);
    fprintf(file, "%ld\n", matr->num_cols);
    fprintf(file, "%ld\n", matr->length_values);

    for (uint64_t i = 0; i < matr->length_values; ++i)
    {
        fprintf(file, "%f ", matr->values[i]);
    }
    fprintf(file, "\n");
    for (uint64_t i = 0; i < matr->length_values; ++i)
    {
        fprintf(file, "%ld ", matr->col_idx[i]);
    }
    fprintf(file, "\n");

    fprintf(file, "%ld\n", matr->length_col_ptr);
    for (uint64_t i = 0; i < matr->length_col_ptr; ++i)
    {
        fprintf(file, "%ld ", matr->col_ptr[i]);
    }
    fprintf(file, "\n");
    for (uint64_t i = 0; i < matr->num_rows; ++i)
    {
        fprintf(file, "%ld ", matr->permutation[i]);
    }
    fclose(file);
}

/*
 * Initializes the necessary memory and values for a COOS matrix
 * to store the result of multiplying c1 and c2.
 */
coos_matrix* initialize_coos_result(coos_matrix* c1, coos_matrix* c2)
{
    uint64_t nnz = 0;
    coos_matrix* res =  malloc(sizeof(coos_matrix));
    if(check_malloc(res))return NULL;
    if(c1->num_entries * c2->num_entries > c1->num_rows * c2->num_cols)
    {
        nnz = c1->num_rows * c2->num_cols;
    }
    else
    {
        nnz = c1->num_entries * c2->num_entries;
    }
    res->entries = malloc(nnz * sizeof(coos_triple));
    if(check_malloc(res->entries))return NULL;
    res->num_rows = c1->num_rows;
    res->num_cols = c2->num_cols;
    res->num_entries = 0;

    return res;
}
/*
 * Initializes the necessary memory and values for a JDS matrix
 * to store the result of multiplying j1 and j2
 */
jds_matrix* initialize_jds_result(jds_matrix* j1,jds_matrix* j2)
{
    uint64_t nnz = 0;
    jds_matrix* res = malloc(sizeof(jds_matrix));
    if(check_malloc(res))return NULL;
    res->num_rows = j1->num_rows;
    res->num_cols = j2->num_cols;
    res->length_values = 0;
    res->length_col_ptr = 0;
    if(j1->length_values * j2->length_values > j1->num_rows * j2->num_cols)
    {
        nnz = j1->num_rows * j2->num_cols;
    }
    else
    {
        nnz = j1->length_values * j2->length_values;
    }
    res->values = malloc(nnz * sizeof(float));
    if(check_malloc(res->values))return NULL;
    res->col_idx = malloc(nnz * sizeof(uint64_t));
    if(check_malloc(res->col_idx))return NULL;
    res->col_ptr = malloc((res->num_cols + 1)  * sizeof(uint64_t));
    if(check_malloc(res->col_ptr))return NULL;
    res->permutation  = malloc(j1->num_rows * sizeof(uint64_t));
    if(check_malloc(res->permutation))return NULL;

    for(uint64_t i = 0; i < res->num_rows; ++i)
    {
        res->permutation[i] = i;
    }

    return res;
}

/*
 * Creates a COOS map using the given COOS matrix.
 * It uses the row index of each element as the key to store the pair with the column index and value
 */
map* to_map_col(coos_matrix* coos)
{
    sort_triples_col_c(coos);

    uint64_t current_col = 0;
    uint64_t columns_in_use = 0;

    if(coos->num_entries != 0)
    {
        current_col = coos->entries[0].col_idx;
        columns_in_use = 1;
    }
    for (uint64_t i = 0; i < coos->num_entries; ++i)
    {
        if(coos->entries[i].col_idx != current_col)
        {
            current_col = coos->entries[i].col_idx;
            columns_in_use += 1;
        }
    }

    map* res = malloc(sizeof(map));
    if(check_malloc(res))return NULL;

    res->num_elements = columns_in_use;
    res->array = malloc(sizeof(map_entry)*columns_in_use);
    if(check_malloc(res->array))return NULL;
    /*
     * Stores the number of elements in the COOS matrix at each column index
     */
    uint64_t* elements_in_column = calloc(columns_in_use, sizeof(uint64_t));
    if(check_malloc(elements_in_column))return NULL;
    current_col = 0;
    uint64_t array_iterator = 0;
    if(coos->num_entries != 0)
    {
        current_col = coos->entries[0].col_idx;
        res->array[array_iterator].idx = current_col;
    }

    for (uint64_t i = 0; i < coos->num_entries; ++i)
    {
        if(coos->entries[i].col_idx != current_col)
        {
            array_iterator += 1;
            current_col = coos->entries[i].col_idx;
            res->array[array_iterator].idx = current_col;
            elements_in_column[array_iterator] += 1;
        }
        else
        {
            elements_in_column[array_iterator] += 1;
        }
    }
    /*
     * This is used to separate the coos->entries array into columns
     * so that they can be iterated through in the inner for-Loop.
     * The value is equal to the index in coos->entries where a new column starts.
     */
    uint64_t entry_offset = 0;
    for (uint64_t i = 0; i < columns_in_use; ++i)
    {
        res->array[i].num_elements = elements_in_column[i];
        res->array[i].pairs = malloc(elements_in_column[i]* sizeof(pair));
        if(check_malloc(res->array[i].pairs))return NULL;
        for (uint64_t j = 0; j < elements_in_column[i]; ++j)
        {
            res->array[i].pairs[j].first = coos->entries[entry_offset+j].row_idx;
            res->array[i].pairs[j].second = coos->entries[entry_offset+j].val;
        }
        entry_offset += elements_in_column[i];
    }

    return res;
}
/*
 * Creates a COOS map using the given COOS matrix.
 * It uses the column index of each element as the key to store the pair with the row index and value
 */
map* to_map_row(coos_matrix* coos)
{
    sort_triples_row_c(coos);

    uint64_t current_row = 0;
    uint64_t rows_in_use = 0;

    if(coos->num_entries != 0)
    {
        current_row = coos->entries[0].row_idx;
        rows_in_use = 1;
    }
    for (uint64_t i = 0; i < coos->num_entries; ++i)
    {
        if(coos->entries[i].row_idx != current_row)
        {
            current_row = coos->entries[i].row_idx;
            rows_in_use += 1;
        }
    }

    map* res = malloc(sizeof(map));
    if(check_malloc(res))return NULL;

    res->num_elements = rows_in_use;
    res->array = malloc(sizeof(map_entry)*rows_in_use);
    if(check_malloc(res->array))return NULL;
    /*
     * Stores the number of elements in the COOS matrix at each row index
     */
    uint64_t* elements_in_row = calloc(rows_in_use, sizeof(uint64_t));
    if(check_malloc(elements_in_row))return NULL;
    current_row = 0;
    uint64_t array_iterator = 0;
    if(coos->num_entries != 0)
    {
        current_row = coos->entries[0].row_idx;
        res->array[array_iterator].idx = current_row;
    }

    for (uint64_t i = 0; i < coos->num_entries; ++i)
    {
        if(coos->entries[i].row_idx != current_row)
        {
            array_iterator += 1;
            current_row = coos->entries[i].row_idx;
            res->array[array_iterator].idx = current_row;
            elements_in_row[array_iterator] += 1;
        }
        else
        {
            elements_in_row[array_iterator] += 1;
        }
    }
    /*
     * This is used to separate the coos->entries array into rows
     * so that they can be iterated through in the inner for-Loop.
     * The value is equal to the index in coos->entries where a new row starts.
     */
    uint64_t entry_offset = 0;
    for (uint64_t i = 0; i < rows_in_use; ++i)
    {
        res->array[i].num_elements = elements_in_row[i];
        res->array[i].pairs = malloc(elements_in_row[i]* sizeof(pair));
        if(check_malloc(res->array[i].pairs))return NULL;
        for (uint64_t j = 0; j < elements_in_row[i]; ++j)
        {
            res->array[i].pairs[j].first = coos->entries[entry_offset+j].col_idx;
            res->array[i].pairs[j].second = coos->entries[entry_offset+j].val;
        }
        entry_offset += elements_in_row[i];
    }

    return res;
}

/*
 * Frees the memory of a coos_matrix
 */
void free_coos(coos_matrix* c1)
{
    if(c1 != NULL)
    {
        free(c1->entries);
        free(c1);
    }
}
/*
 * Frees the memory of a COOS map
 */
void free_map(map* map1)
{
    for (uint64_t i = 0; i < map1->num_elements; ++i)
    {
        free(map1->array[i].pairs);
    }
    free(map1->array);
    free(map1);
}
/*
 * Frees the memory of a jds_matrix
 */
void free_jds(jds_matrix* j1)
{
    if(j1 != NULL)
    {
        free(j1->values);
        free(j1->col_idx);
        free(j1->col_ptr);
        free(j1->permutation);
        free(j1);
    }
}
/*
 * Frees the memory of a 2D float-array
 */
void free_2d_matr(float** f1, uint64_t num_rows)
{
    for (uint64_t i = 0; i < num_rows; ++i) {
        free(f1[i]);
    }
    free(f1);
}