#!/bin/bash

#SBATCH --job-name=mpi_scaling_mandelbrot
#SBATCH --nodes=2
#SBATCH --partition=EPYC
#SBATCH --exclusive
#SBATCH --time=02:00:00

# Load the required module
module load openMPI/4.1.5/gnu/12.2.1

# Navigate to the /build/bin directory
cd ./build/bin

# Define variables
max_tasks=8
step=4

# Run the scaling test with tasks from 4 to 254 with a step of 4:
# (Run the program with 4 to 254 MPI tasks and 1 OMP thread per task)
for ((tasks=4; tasks<=$max_tasks; tasks+=$step))
do
    echo "----------------------------------------------------------------------------------------------------------------------------------"
    echo "Running mandelbrot with $tasks tasks"
    # Set the number of OMP threads per MPI task to 1
    export OMP_NUM_THREADS=1
    # Run the program (pass the desired arguments to the program)
    # (Choose the mapping policy (core, socket, node))
    mpirun -np $tasks --map-by core ./main 1000 1000 -2.75 -2.0 1.25 2.0 65535
    echo "----------------------------------------------------------------------------------------------------------------------------------"
done

# Print the completion message
echo "----------------------------------------------------------------------------------------------------------------------------------"
echo "Job completed!"
echo "----------------------------------------------------------------------------------------------------------------------------------"
echo "Run on partition: $SLURM_JOB_PARTITION"
echo "Host: $SLURM_JOB_NODELIST"
echo "----------------------------------------------------------------------------------------------------------------------------------"
