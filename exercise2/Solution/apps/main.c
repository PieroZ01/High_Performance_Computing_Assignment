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
  short int **global_M = NULL;
  // The master process allocates the memory for the global matrix M
  if (rank == 0)
  {
    // Alllocate a multidimensional array of integers (short int)
    global_M = (short int **)malloc(n_y * sizeof(short int *));
    for (int p = 0; p < n_y; ++p)
    {
      global_M[p] = (short int *)malloc(n_x * sizeof(short int));
    }
  }

  // Define the local part of the matrix M on each process
  // (Allocate only the memory for the local part of the matrix M on each process)
  short int **local_M = (short int **)malloc(local_rows * sizeof(short int *));
  for (int q = 0; q < local_rows; ++q)
  {
    local_M[q] = (short int *)malloc(n_x * sizeof(short int));
  }
  // (Each thread will compute a part of the local part of the global matrix M)
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
        for (int i = 0; i < n_x; ++i)
        {
          double complex c = x_L + i * dx + y * I;
          local_M[j][i] = mandelbrot(c, I_max);
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

  // Gather the results to the master process and start the timer to measure the communication time
  if (rank == 0)
  {
    timer = MPI_Wtime();
  }

  // Each process sends its local part of the matrix M to the master process
  
  // Loop over the number of round
  for (int r = 0; r < rows_per_process; ++r)
  {
    // Position variable
    int position = 0;

    // Each process sends to the master process it's r-th row of the local part of the matrix M
    MPI_Send(local_M[r], n_x, MPI_SHORT, 0, r, MPI_COMM_WORLD);

    // The master process receives the r-th row of the local part of the matrix M from each process
    // and stores it in the correct position in the global matrix M
    if (rank == 0)
    {
      for (int p = 0; p < size; ++p)
      {
        MPI_Recv(global_M[position + p], n_x, MPI_SHORT, p, r, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }

    // Increment the position
    position += size;
      
  }

  // Manage case where the number of rows is not divisible by the number of processes
  if (remaining_rows > 0)
  {
    int position = n_y - remaining_rows;

    if (rank < remaining_rows)
    {
      // The process sends to the master process it's last row of the local part of the matrix M
      MPI_Send(local_M[local_rows - 1], n_x, MPI_SHORT, 0, rows_per_process, MPI_COMM_WORLD);
    }

    // The master process receives the last row of the local part of the matrix M from each process
    // and stores it in the correct position in the global matrix M
    if (rank == 0)
    {
      for (int p = 0; p < remaining_rows; ++p)
      {
        MPI_Recv(global_M[position + p], n_x, MPI_SHORT, p, rows_per_process, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }
  }

  // Stop the timer to measure the communication time
  if (rank == 0)
  {
    communication_time = MPI_Wtime() - timer;
  }

  // Free the memory for the local part of the matrix M on each process
  for (int q = 0; q < local_rows; ++q)
  {
    free(local_M[q]);
  }
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
    // Flatten the global matrix M into a 1D array to write the image to a pgm file
    short int *global_M_flat = (short int *)malloc(n_x * n_y * sizeof(short int));
    for (int j = 0; j < n_y; ++j)
    {
      for (int i = 0; i < n_x; ++i)
      {
        global_M_flat[j * n_x + i] = global_M[j][i];
      }
    }

    write_pgm_image(global_M_flat, I_max, n_x, n_y, "mandelbrot.pgm");

    // Free the memory for the global matrix M
    for (int j = 0; j < n_y; ++j)
    {
      free(global_M[j]);
    }
    free(global_M);
  }

  // Finalize MPI
  MPI_Finalize();

  return 0;
}
