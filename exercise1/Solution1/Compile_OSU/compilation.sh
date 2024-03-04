#!/bin/bash

#SBATCH --job-name=compilation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH -A dssc
#SBATCH -p EPYC
#SBATCH --output=compilation.out
#SBATCH --error=compilation.err
#SBATCH --mem-per-cpu=1500MB

# Load the required modules
module load openMPI/4.1.5/gnu/12.2.1

# Enter the OSU benchmark directory
cd ./osu-micro-benchmarks-7.3
make

# Compile the OSU benchmark
./configure CC=$(which mpicc) CXX=$(which mpicxx)
make

# Write on the output file for which partition the job was submitted
echo "-------------------------------------------------"
echo "Compilation of the OSU benchmark completed"
echo "-------------------------------------------------"
echo "Run on partition: $SLURM_JOB_PARTITION"
echo "Host: $SLURM_JOB_NODELIST"
