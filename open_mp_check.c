#define _GNU_SOURCE

#include <stdio.h>
#include <omp.h>
#include <sched.h>


int main(int argc, char *argv[])
{
#pragma omp parallel
    {
        int number_of_cores, thread_id, number_of_threads, max_threads, core_id;
        number_of_cores = omp_get_num_procs();
        thread_id = omp_get_thread_num();
        number_of_threads = omp_get_num_threads();
        max_threads = omp_get_max_threads();
        core_id = sched_getcpu();

        if (thread_id == 0) {
            printf("%i : number_of_cores = %i\n", thread_id, number_of_cores);
            printf("%i : max_threads = %i\n", thread_id, max_threads);
            printf("%i : thread_id = %i\n", thread_id, number_of_threads);
        }
        printf("%i : Thread %i out of %i, %i\n", thread_id, thread_id, number_of_threads, core_id);
    }
    return(0);
}