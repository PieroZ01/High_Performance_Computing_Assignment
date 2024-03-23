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

  // Gather the results to the master process and start the timer to measure the communication time
  if (rank == 0)
  {
    global_M = (short int *)malloc(n_x * n_y * sizeof(short int));
    timer = MPI_Wtime();
  }

  // Each process sends its local part of the matrix M to the master process
  // row by row, in order to make sure that the data is received in the correct order by the master process
  // and that the global matrix M is correctly reconstructed (in fact, each process' rows are not contiguous
  // in the global matrix M. They are distributed in a round robin fashion among the processes, so the master
  // process needs to receive the rows in the correct order to reconstruct the global matrix M correctly, row by row,
  // following the precise order used to distribute the rows among the processes)
  // (We use MPI_Ssend to make sure that the data is sent in the correct order)
  
  // Loop over the number of rounds - 1
  for (int r = 0; r < rows_per_process; ++r)
  {
    // Loop over the number of processes (each process, identified by its rank, sends its row to the master process)
    for (int p = 0; p < size; ++p)
    {
      // If the process' rank is equal to the current p, it sends its r row to the master process
      // (The r row of the matrix local_M is selected by the formula: r * n_x)
      if (rank == p)
      {
        // Print the rank of the process and the row to be sent
        printf("Rank: %d, Row: %d\n", rank, r);

        MPI_Ssend(local_M + r * n_x, n_x, MPI_SHORT, 0, r * p, MPI_COMM_WORLD);
      }

      // If the process' rank is equal to 0, it receives the r row from the process p and stores it in the global matrix M,
      // in the correct position, which is the r + p row of the matrix global_M
      if (rank == 0)
      {
        MPI_Recv(global_M + (r + p) * n_x, n_x, MPI_SHORT, p, r * p, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }
  }

  // Loop over the number of remaining rows: the first remaining_rows processes send their (rows_per_process + 1) row
  // to the master process
  if (remaining_rows != 0)
  {
    for (int q = 0; q < remaining_rows; ++q)
    {
      if (rank == q)
      {
        MPI_Ssend(local_M + rows_per_process * n_x, n_x, MPI_SHORT, 0, q, MPI_COMM_WORLD);
      }

      if (rank == 0)
      {
        MPI_Recv(global_M + (rows_per_process + q) * n_x, n_x, MPI_SHORT, q, q, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }
  }
  
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

    // Free the memory for the global matrix M
    free(global_M);
  }

  // Finalize MPI
  MPI_Finalize();

  return 0;
}
