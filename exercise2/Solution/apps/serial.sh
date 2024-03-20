#!/bin/bash

#SBATCH --job-name=serial_mandelbrot
#SBATCH --nodes=1
#SBATCH --partition=EPYC
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --ntasks=1

# Load the required module
module load openMPI/4.1.5/gnu/12.2.1

# Navigate to the /build/bin directory
cd ./build/bin

# Run the code with 1 MPI task and 1 OMP thread (serial execution)
echo "----------------------------------------------------------------------------------------------------------------------------------"
echo "Running mandelbrot with 1 MPI task and 1 OMP thread"
export OMP_NUM_THREADS=1
mpirun -np 1 --map-by socket --bind-to socket ./main 1000 1000 -2.0 -2.0 2.0 2.0 1000
echo "----------------------------------------------------------------------------------------------------------------------------------"

# Print the completion message
echo "----------------------------------------------------------------------------------------------------------------------------------"
echo "Job completed!"
echo "----------------------------------------------------------------------------------------------------------------------------------"
echo "Run on partition: $SLURM_JOB_PARTITION"
echo "Host: $SLURM_JOB_NODELIST"
echo "----------------------------------------------------------------------------------------------------------------------------------"
