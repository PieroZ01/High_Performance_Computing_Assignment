// Libraries
#include <mpi.h>

// Header
#include "mandelbrot.h"

// Main function
int main(int argc, char *argv[])
{
  // Initialize MPI
  int mpi_provided_thread_level;
  MPI_Init_threads(&argc, &argv, MPI_THREAD_FUNNELED, &mpi_provided_thread_level);
  if (mpi_provided_thread_level < MPI_THREAD_FUNNELED)
  {
    printf("A problem arise when asking for MPI_THREAD_FUNNELED level\n");
    MPI_Finalize();
    exit(1);
  }

  // Get the number of processes and the rank of the process
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Variables to be read from command line arguments, with default values
  int n_x = argc > 1 ? atoi(argv[1]) : XWIDTH;
  int n_y = argc > 2 ? atoi(argv[2]) : YWIDTH;
  double x_L = argc > 3 ? atof(argv[3]) : -2.0;
  double y_L = argc > 4 ? atof(argv[4]) : -2.0;
  double x_R = argc > 5 ? atof(argv[5]) : 2.0;
  double y_R = argc > 6 ? atof(argv[6]) : 2.0;
  int I_max = argc > 7 ? atoi(argv[7]) : MAXVAL;

  // Delta x and y
  double dx = (x_R - x_L) / n_x;
  double dy = (y_R - y_L) / n_y;

  // Get the local amount of rows to be computed by each process
  int rows_per_process = n_y / size;
  int start_row = rank * rows_per_process;
  int end_row = (rank == size - 1) ? n_y : start_row + rows_per_process;
  int local_rows = end_row - start_row;

  // Define the 2D matrix M of integers (short int) whose entries [j][i] are the image's pixels
  // (Allocate only the memory for the local part of the matrix M on each process)
  short int *local_M = (short int *)malloc(n_x * local_rows * sizeof(short int));

  // Compute the mandelbrot set
  #pragma omp parallel for schedule(dynamic)
    for (int j = 0; j < local_rows; ++j)
    {
        double y = y_L + (start_row + j) * dy;
        for (int i = 0; i < n_x; ++i)
        {
          double complex c = x_L + i * dx + y * I;
          local_M[j * n_x + i] = mandelbrot(c, I_max);
        }
    }

  // Define the global matrix M to gather the results from all the processes
  short int *global_M = NULL;

  // Gather the results to the master process
  if (rank == 0)
  {
    global_M = (short int *)malloc(n_x * n_y * sizeof(short int));
  }
  MPI_Gather(local_M, n_x * local_rows, MPI_SHORT, global_M, n_x * local_rows, MPI_SHORT, 0, MPI_COMM_WORLD);

  // The master process writes the image to a pgm file
  if (rank == 0)
  {
    write_pgm_image(global_M, I_max, n_x, n_y, "mandelbrot.pgm");
    free(global_M);
  }

  // Free the memory
  free(local_M);

  // Finalize MPI
  MPI_Finalize();

  return 0;
}