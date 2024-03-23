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

  // Get the local amount of rows to be computed by each process
  // (The rows are distributed among the processes in a round robin fashion with the first
  // processes getting one more row if the number of rows is not divisible by the number of processes)
  const int rows_per_process = n_y / size;
  const int start_row = rank;
  const int remaining_rows = n_y % size;
  const int local_rows = (rank < remaining_rows) ? rows_per_process + 1 : rows_per_process;

  // Define the 2D matrix M of integers (short int) whose entries [j][i] are the image's pixels
  // (Allocate only the memory for the local part of the matrix M on each process)
  short int *local_M = (short int *)malloc(n_x * local_rows * sizeof(short int));
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

    // Compute the local part of the matrix M on each process
    // (Each thread will compute a part of the local part of the matrix M;
    // each process is assigned its rows to be computed in a round robin fashion)
    #pragma omp parallel for schedule(dynamic)
      for (int j = 0; j < local_rows; ++j)
      {
        const double y = y_L + (start_row + j * size) * dy;
        const int index = j * n_x;
        for (int i = 0; i < n_x; ++i)
        {
          double complex c = x_L + i * dx + y * I;
          local_M[index + i] = mandelbrot(c, I_max);
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
  
  // Define the variables receivedcounts and displs to be used in the MPI_Gatherv function
  int *receivedcounts = (int *)malloc(size * sizeof(int));
  int *displs = (int *)malloc(size * sizeof(int));
  for (int r = 0; r < size; ++r)
  {
    receivedcounts[r] = (r < remaining_rows) ? (rows_per_process + 1) * n_x : rows_per_process * n_x;
    displs[r] = r * rows_per_process * n_x + (r < remaining_rows ? r : remaining_rows) * n_x;
  }

  // Gather the results to the master process and start the timer to measure the communication time
  if (rank == 0)
  {
    global_M = (short int *)malloc(n_x * n_y * sizeof(short int));
    timer = MPI_Wtime();
  }

  // Gather the local parts of the matrix M from all the processes to the master process
  // (We need to use MPI_Gatherv because the number of rows computed by each process could
  // be different and the master process needs to know the number of rows computed by each process:
  // some process could have computed one more row than the others)

  MPI_Gatherv(local_M, n_x * local_rows, MPI_SHORT, global_M, receivedcounts, displs, MPI_SHORT, 0, MPI_COMM_WORLD);

  // Stop the timer to measure the communication time
  if (rank == 0)
  {
    communication_time = MPI_Wtime() - timer;
  }

  // Free the memory for the local part of the matrix M on each process
  free(local_M);

  // Free the memory for the receivedcounts and displs arrays
  free(receivedcounts);
  free(displs);

  // The master process reorders the matrix M to write the image to a pgm file: the global matrix M is
  // now composed of blocks of rows computed by each process in a round robin fashion; we need to reorder
  // the rows to write the image correctly: each row has to be placed in the correct position in the image, 
  // so the first row of process 0 has to be placed in the first row of the image, the first row of process 1
  // has to be placed in the second row of the image, and so on
  // (This is to be done only if the number of processes is greater than 1, otherwise the matrix M is already ordered)
  if (size>1)
  {
    if (rank == 0)
    {
      // Reorder the matrix M
      short int *reordered_M = (short int *)malloc(n_x * n_y * sizeof(short int));
      for (int r = 0; r < size; ++r)
      {
        const int start_row = r * rows_per_process + min(r, remaining_rows);
        const int local_rows = (r < remaining_rows) ? rows_per_process + 1 : rows_per_process;
       for (int j = 0; j < local_rows; ++j)
       {
        for (int i = 0; i < n_x; ++i)
        {
          reordered_M[(start_row + j) * n_x + i] = global_M[(r * rows_per_process + j) * n_x + i];
        }
       }
      }

      // Free the memory for the global matrix M
      free(global_M);

      global_M = reordered_M;

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
    // Write the image to a pgm file
    write_pgm_image(global_M, I_max, n_x, n_y, "mandelbrot.pgm");

    // Free the memory for the reordered matrix M
    free(global_M);
  }

  // Finalize MPI
  MPI_Finalize();

  return 0;
}
