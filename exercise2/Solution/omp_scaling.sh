#!/bin/bash

#SBATCH --job-name=omp_scaling_mandelbrot
#SBATCH --nodes=1
#SBATCH --partition=EPYC
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --ntasks=1

# Load the required OpenMP module
module load openMPI/4.1.5/gnu/12.2.1

# Navigate to the /build/bin directory
cd ./build/bin

# Define variables
max_threads=128
step=2

# Run the scaling test with threads from 2 to 128 with a step of 2:
# (Run the program with 1 MPI task and 2 to 128 OMP threads)
for ((threads=2; threads<=$max_threads; threads+=$step))
do
    echo "----------------------------------------------------------------------------------------------------------------------------------"
    echo "Running mandelbrot with $threads threads"
    # Export OpenMP environment variables
    export OMP_NUM_THREADS=$threads
    #export OMP_PROC_BIND=spread
    #export OMP_PLACES=threads
    # Run the program (pass the desired arguments to the program) with a single MPI task
    mpirun ./main 1000 1000 -2.0 -2.0 2.0 2.0 1000
    echo "----------------------------------------------------------------------------------------------------------------------------------"
done

# Print the completion message
echo "----------------------------------------------------------------------------------------------------------------------------------"
echo "Job completed!"
echo "----------------------------------------------------------------------------------------------------------------------------------"
echo "Run on partition: $SLURM_JOB_PARTITION"
echo "Host: $SLURM_JOB_NODELIST"
echo "----------------------------------------------------------------------------------------------------------------------------------"
