#!/bin/bash

#SBATCH --job-name=omp_scaling_mandelbrot
#SBATCH --nodes=1
#SBATCH --partition=EPYC
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64

# Load the required module
module load openMPI/4.1.5/gnu/12.2.1

# Navigate to the /build/bin directory
cd ./build/bin

# Define variables
max_threads=128
step=2

# Run the scaling test with threads from 2 to 128 with a step of 2:
# (Run the program with 1 MPI task and 2 to 128 OMP threads)
for ((threads_n=2; threads_n<=$max_threads; threads_n+=$step))
do
    echo "----------------------------------------------------------------------------------------------------------------------------------"
    echo "Running mandelbrot with $threads_n threads"
    # Export OpenMP environment variables
    export OMP_NUM_THREADS=$threads_n
    # Choose the places and the binding policy (threads, cores, sockets - close, spread)
    export OMP_PLACES=cores
    export OMP_PROC_BIND=spread
    # Run the program (pass the desired arguments to the program) with a single MPI task
    mpirun -np 1 --map-by socket --bind-to socket ./main 1000 1000 -2.0 -2.0 2.0 2.0 1000
    echo "----------------------------------------------------------------------------------------------------------------------------------"
done

# Print the completion message
echo "----------------------------------------------------------------------------------------------------------------------------------"
echo "Job completed!"
echo "----------------------------------------------------------------------------------------------------------------------------------"
echo "Run on partition: $SLURM_JOB_PARTITION"
echo "Host: $SLURM_JOB_NODELIST"
echo "----------------------------------------------------------------------------------------------------------------------------------"
