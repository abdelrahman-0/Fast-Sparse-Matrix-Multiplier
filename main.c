#include "coos.h"
#include "jds.h"

int main(int argc, char** argv)
{
    /*
     * Check if the given argument is (-h/--help) or if there are no arguments at all
     */
    if(argc == 1 || (argc == 2 && (strcmp(argv[1],"-h") == 0 || (strcmp(argv[1],"--help") == 0))))
    {
        print_help();
        return 0;
    }

    char c;

    int c_flag = 0;
    int j_flag = 0;
    int b_flag = 0;

    uint64_t iterations = 1;

    coos_matrix* c1 = NULL,*c2 = NULL;
    jds_matrix* j1 = NULL,*j2 = NULL;
    /*
     * Parse the given arguments
     *  -c To multiply two COOS matrices
     *  -j To multiply two JDS matrices
     *  -b Tto output a benchmark if its also given a multiplication argument
     *     (-c, -j, or both)
     */
    while ((c = getopt(argc, argv, "cjb")) != -1)
    {
        switch (c)
        {
            case 'c':
                if(optind + 1 < argc)
                {
                    c_flag = 1;
                    /*
                     * COOS matrix multiplication
                     */
                    if((c1 = read_input_coos(argv[optind])) == NULL) break;
                    if((c2 = read_input_coos(argv[optind+1])) == NULL) break;

                    if(c1->num_rows != c2->num_cols)
                    {
                        perror("Invalid multiplication");
                        break;
                    }

                    coos_matrix *res = initialize_coos_result(c1,c2);

                    if(c1->num_entries == 0 || c2->num_entries == 0)
                    {
                        coos_output(res);
                        break;
                    }

                    if(res != NULL)
                    {
                        matr_mult_coos(c1, c2, res);
                        if(res->num_rows != 0)
                        {
                            coos_output(res);
                        }
                    }

                    free_coos(res);
                }
                else
                {
                    perror("not enough Arguments for -c please use -h/--help for more information");
                }
                break;
            case 'j':
                if(optind + 1 < argc)
                {
                    j_flag = 1;
                    /*
                     * JDS matrix multiplication
                     */
                    if((j1 = read_input_jds(argv[optind])) == NULL) break;
                    if((j2 = read_input_jds(argv[optind+1])) == NULL) break;

                    if(j1->num_rows != j2->num_cols)
                    {
                        perror("Invalid multiplication");
                        break;
                    }
                    jds_matrix* res = initialize_jds_result(j1,j2);
                    if(j1->length_values == 0 || j2->length_values == 0)
                    {
                        jds_output(res);
                        break;
                    }
                    if(res != NULL)
                    {
                        matr_mult_jds(j1,j2,res);
                        if(res->num_rows != 0)
                        {
                            jds_output(res);
                        }
                    }

                    free_jds(res);
                }
                else
                {
                    perror("not enough Arguments for -j please use -h/--help for more information");
                }
                break;
            case 'b':
                b_flag = 1;
                char *c;
                if(optind < argc && strtol(argv[optind], &c, 10) > 0)
                {
                    iterations = strtol(argv[optind], &c, 10);
                }
                break;
            case '?':
                printf("unknown option: %c\n", optopt);
                break;
        }
    }
    /*
     * If the b_flag is set Benchmark the multiplication methods that are also set
     */
    if(b_flag == 1)
    {
        if(c_flag == 1)
        {
            if(c1 != NULL && c2 != NULL)
            {
                benchmark_coos(c1,c2,iterations);
            }
        }
        if(j_flag == 1)
        {
            if(j1 != NULL && j2 != NULL)
            {
                benchmark_jds(j1,j2,iterations);
            }
        }
    }

    free_coos(c1);
    free_coos(c2);
    free_jds(j1);
    free_jds(j2);

    return 0;
}


