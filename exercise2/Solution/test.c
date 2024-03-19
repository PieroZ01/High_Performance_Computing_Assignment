// Test to see if openMP is working and the program scales with the number of threads
// We simply sum the numbers from 1 to 1000000000

#include <stdio.h>
#include <omp.h>

int main() {
    int i;
    double sum = 0;
    double start, end;

    start = omp_get_wtime();
    #pragma omp parallel for reduction(+:sum)
    for (i = 1; i <= 1000000000; i++) {
        sum += i;
    }
    end = omp_get_wtime();
    
    printf("Time: %f\n", end - start);
    return 0;
}


