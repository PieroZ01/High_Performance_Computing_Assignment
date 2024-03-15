#!/bin/bash

#SBATCH --job-name=omp_scaling_mandelbrot
#SBATCH --nodes=1
#SBATCH --partition=EPYC
#SBATCH --exclusive
#SBATCH --time=02:00:00

# Load the required OpenMP module
module load openMPI/4.1.5/gnu/12.2.1

# Navigate to the /build/bin directory
cd ./build/bin

# Define the output file path
output_csv="../../Results/omp_scaling_mandelbrot.csv"

# Go to the csv output file and write the header
echo "Threads,Time" > $output_csv

# Define variables
max_threads=128
step=2

# Run the scaling test with threads from 2 to 128 with a step of 2:
# write on the csv file the number of threads and the time taken to execute the program
for ((threads=2; threads<=$max_threads; threads+=$step))
do
    echo "----------------------------------------------------------------------------------------------------------------------------------"
    echo "Running mandelbrot with $threads threads"
    export OMP_NUM_THREADS=$threads
    start_time=$(date +%s.%N)
    mpirun -np 1 ./main
    end_time=$(date +%s.%N)
    elapsed_time=$(echo "$end_time - $start_time" | bc)
    # Write the results to the csv file
    echo "$threads,$elapsed_time" >> $output_csv
done

# Print the completion message
echo "----------------------------------------------------------------------------------------------------------------------------------"
echo "Job completed!"
echo "----------------------------------------------------------------------------------------------------------------------------------"
echo "Run on partition: $SLURM_JOB_PARTITION"
echo "Host: $SLURM_JOB_NODELIST"
echo "----------------------------------------------------------------------------------------------------------------------------------"
