// Libraries
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

// Header
#include "mandelbrot.h"

// Main function
int main(int argc, char *argv[])
{
  // Initialize MPI
  int mpi_provided_thread_level;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mpi_provided_thread_level);
  if (mpi_provided_thread_level < MPI_THREAD_FUNNELED)
  {
    printf("A problem arised when asking for MPI_THREAD_FUNNELED level\n");
    MPI_Finalize();
    exit(1);
  }

  // Get the number of processes and the rank of the process
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Variables to be read from command line arguments, with default values
  const int n_x = argc > 1 ? atoi(argv[1]) : XWIDTH;
  const int n_y = argc > 2 ? atoi(argv[2]) : YWIDTH;
  const double x_L = argc > 3 ? atof(argv[3]) : -2.0;
  const double y_L = argc > 4 ? atof(argv[4]) : -2.0;
  const double x_R = argc > 5 ? atof(argv[5]) : 2.0;
  const double y_R = argc > 6 ? atof(argv[6]) : 2.0;
  const int I_max = argc > 7 ? atoi(argv[7]) : MAXVAL;

  // Number of iterations and timer variable
  const int n_iter = 6;
  double timer = 0.0;

  // Define the array containing the time measurements
  double *time_measures = (double *)malloc(n_iter * sizeof(double));

  // Open a csv file to write the time measurements (NOTE: change the file name accordingly)
  FILE *file = fopen("../../Results/omp_scaling_mandelbrot.csv", "a+");
  if (file == NULL)
  {
    printf("Error opening the file\n");
    exit(1);
  }
  if (rank == 0)
  {
    fseek(file, 0, SEEK_END);
    if (ftell(file) == 0)
    {
      fprintf(file, "\"n_processes\", \"n_threads\", \"n_x\", \"n_y\", \"I_max\", \"Average time (s)\", \"Std. deviation (s)\"\n");
      fflush(file);
    }
  }

  // Delta x and y
  const double dx = (x_R - x_L) / n_x;
  const double dy = (y_R - y_L) / n_y;

  // Get the local amount of rows to be computed by each process
  const int rows_per_process = n_y / size;
  const int start_row = rank * rows_per_process;
  const int end_row = (rank == size - 1) ? n_y : start_row + rows_per_process;
  const int local_rows = end_row - start_row;

  // Define the 2D matrix M of integers (short int) whose entries [j][i] are the image's pixels
  // (Allocate only the memory for the local part of the matrix M on each process)
  short int *local_M = (short int *)malloc(n_x * local_rows * sizeof(short int));

  // Loop over the number of iterations
  for (int p = 0; p < n_iter; ++p)
  {
    // Sincrhonize all the processes before starting the computation
    MPI_Barrier(MPI_COMM_WORLD);

    // The master process measures the time (start the timer)
    if (rank == 0)
    {
      timer = MPI_Wtime();
    }

    // Compute the mandelbrot set
    #pragma omp parallel for schedule(dynamic)
      for (int j = 0; j < local_rows; ++j)
      {
        double y = (y_L + (start_row + j) * dy) * I;
        for (int i = 0; i < n_x; ++i)
        {
          double complex c = x_L + i * dx + y;
          local_M[j * n_x + i] = mandelbrot(c, I_max);
        }
     }
  
    // Wait for all the processes to finish the computation of the mandelbrot set
    MPI_Barrier(MPI_COMM_WORLD);

    // The master process measures the time (stop the timer)
    if (rank == 0)
    {
      time_measures[p] = MPI_Wtime() - timer;
    }
  }

  // The master process writes the results to the csv file
  if (rank == 0)
  {
    #pragma omp parallel
    {
      #pragma omp master
      {
      int n_threads = omp_get_num_threads(); // Get the number of threads
      fprintf(file, "%d, %d, %d, %d, %d, %f, %f\n", size, n_threads, n_x, n_y, I_max, mean(time_measures, n_iter), std_dev(time_measures, n_iter));
      fflush(file);
      }
    }
  }

  // Close the csv file
  fclose(file);

  // Free the memory for the time measurements
  free(time_measures);

  // Define the global matrix M to gather the results from all the processes
  short int *global_M = NULL;

  // Gather the results to the master process
  if (rank == 0)
  {
    global_M = (short int *)malloc(n_x * n_y * sizeof(short int));
  }

  MPI_Gather(local_M, n_x * local_rows, MPI_SHORT, global_M, n_x * local_rows, MPI_SHORT, 0, MPI_COMM_WORLD);

  // Free the memory for the local part of the matrix M on each process
  free(local_M);

  // The master process writes the image to a pgm file and frees the memory
  if (rank == 0)
  {
    write_pgm_image(global_M, I_max, n_x, n_y, "mandelbrot.pgm");
  }

  // Free the memory for the global matrix M
  free(global_M);

  // Finalize MPI
  MPI_Finalize();

  return 0;
}
