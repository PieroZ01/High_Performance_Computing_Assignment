// Libraries
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

// Header
#include "mandelbrot.h"

// Define the function to compute the mandelbrot set
static inline int mandelbrot(const double complex c, const int max_iter)
{
  double complex z = 0.0;
  int k = 0;
  while (creal(z)*creal(z) + cimag(z)*cimag(z) < 4.0 && k < max_iter)
  {
    z = z*z + c;
    k++;
  }
  return k;
}

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

  // Open a csv file to write the time measurements (NOTE: change the file's name accordingly)
  FILE *file = fopen("../../Results/mpi_scaling_mandelbrot.csv", "a+");
  if (file == NULL)
  {
    printf("Error opening the file\n");
    MPI_Finalize();
    exit(1);
  }
  if (rank == 0)
  {
    // Move the file pointer to the end of the file
    fseek(file, 0, SEEK_END);
    // If the file is empty, write the header
    if (ftell(file) == 0)
    {
      fprintf(file, "\"n_processes\",\"n_threads\",\"n_x\",\"n_y\",\"I_max\",\"Average time\",\"Standard deviation\",\"Communication time\"\n");
      fflush(file);
    }
  }

  // Variables to be read from command line arguments, with default values
  const int n_x = argc > 1 ? atoi(argv[1]) : XWIDTH;
  const int n_y = argc > 2 ? atoi(argv[2]) : YWIDTH;
  const double x_L = argc > 3 ? atof(argv[3]) : -2.75;
  const double y_L = argc > 4 ? atof(argv[4]) : -2.0;
  const double x_R = argc > 5 ? atof(argv[5]) : 1.25;
  const double y_R = argc > 6 ? atof(argv[6]) : 2.0;
  const int I_max = argc > 7 ? atoi(argv[7]) : MAXVAL;

  // Number of iterations to measure the average time and the standard deviation
  const int n_iterations = 6;

  // Time measurements' variables
  double timer = 0.0;
  double *time_taken = (double *)malloc(n_iterations * sizeof(double));
  double communication_time = 0.0;

  // Delta x and y
  const double dx = (x_R - x_L) / n_x;
  const double dy = (y_R - y_L) / n_y;

  // Get the local amount of columns to be computed by each process
  const int columns_per_process = n_x / size;
  const int start_col = rank * columns_per_process;
  const int end_col = (rank == size - 1) ? n_x : start_col + columns_per_process;
  const int local_columns = end_col - start_col;

  // Define the 2D matrix M of integers (short int) whose entries [j][i] are the image's pixels
  // (Allocate only the memory for the local part of the matrix M on each process)
  short int *local_M = (short int *)malloc(n_y * local_columns * sizeof(short int));
  // (Each thread will compute a part of the local part of the matrix M)
  // (The number of threads is defined by the environment variable OMP_NUM_THREADS)

  // Loop over the number of iterations to measure the average time and the standard deviation
  for (int t = 0; t < n_iterations; ++t)
  {
    // Sinchronize all the processes before starting the computation
    if (size > 1)
    {
      MPI_Barrier(MPI_COMM_WORLD);
    }

    // Measure the time (start the timer)
    if (rank == 0)
    {
      timer = MPI_Wtime();
    }

    #pragma omp parallel for schedule(dynamic)
      for (int i = 0; i < local_columns; ++i)
      {
        const double x = x_L + (start_col + i) * dx;
        const int index = i * n_y;
        for (int j = 0; j < n_y; ++j)
        {
          double complex c = x + (y_L + j * dy) * I;
          local_M[index + j] = mandelbrot(c, I_max);
        }
      }

    // Sinchronize all the processes after the computation
    if (size > 1)
    {
      MPI_Barrier(MPI_COMM_WORLD);
    }

    // Measure the time (stop the timer)
    if (rank == 0)
    {
      time_taken[t] = MPI_Wtime() - timer;
    }
  }

  // Define the global matrix M to gather the results from all the processes
  short int *global_M = NULL;

  // Gather the results to the master process and start the timer to measure the communication time
  if (rank == 0)
  {
    global_M = (short int *)malloc(n_x * n_y * sizeof(short int));
    timer = MPI_Wtime();
  }

  MPI_Gather(local_M, n_y * local_columns, MPI_SHORT, global_M, n_y * local_columns, MPI_SHORT, 0, MPI_COMM_WORLD);

  // Stop the timer to measure the communication time
  if (rank == 0)
  {
    communication_time = MPI_Wtime() - timer;
  }

  // Free the memory for the local part of the matrix M on each process
  free(local_M);

  // The master process writes the results to the csv file
  if (rank == 0)
  {
    #pragma omp parallel
    {
      #pragma omp master
      {
      int n_threads = omp_get_num_threads(); // Get the number of threads
      fprintf(file, "%d, %d, %d, %d, %d, %f, %f, %f\n", size, n_threads, n_x, n_y, I_max, mean(time_taken, n_iterations),\
      std_dev(time_taken, n_iterations), communication_time);
      fflush(file);
      }
    }
  }

  // Free the memory for the time measurements
  free(time_taken);

  // Close the csv file
  fclose(file);

  // The master process writes the image to a pgm file in the build/bin directory and frees the memory
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
