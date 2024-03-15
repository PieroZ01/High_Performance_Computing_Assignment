#!/bin/bash

#SBATCH --job-name=mpi_scaling_mandelbrot
#SBATCH --nodes=2
#SBATCH --partition=EPYC
#SBATCH --exclusive
#SBATCH --time=02:00:00

# Load the required OpenMP module
module load openMPI/4.1.5/gnu/12.2.1

# Navigate to the /build/bin directory
cd ./build/bin

# Define the output file path
output_csv="../../Results/mpi_scaling_mandelbrot.csv"

# Go to the csv output file and write the header
echo "Tasks,Time" > $output_csv

# Define variables
max_tasks=256
step=4

# Run the scaling test with tasks from 4 to 256 with a step of 4:
# write on the csv file the number of tasks and the time taken to execute the program
for ((tasks=4; tasks<=$max_tasks; tasks+=$step))
do
    echo "----------------------------------------------------------------------------------------------------------------------------------"
    echo "Running mandelbrot with $tasks tasks"
    export OMP_NUM_THREADS=1
    start_time=$(date +%s.%N)
    mpirun -np $tasks ./main
    end_time=$(date +%s.%N)
    elapsed_time=$(echo "$end_time - $start_time" | bc)
    # Write the results to the csv file
    echo "$tasks,$elapsed_time" >> $output_csv
done

# Print the completion message
echo "----------------------------------------------------------------------------------------------------------------------------------"
echo "Job completed!"
echo "----------------------------------------------------------------------------------------------------------------------------------"
echo "Run on partition: $SLURM_JOB_PARTITION"
echo "Host: $SLURM_JOB_NODELIST"
echo "----------------------------------------------------------------------------------------------------------------------------------"
