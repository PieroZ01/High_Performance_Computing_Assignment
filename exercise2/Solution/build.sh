#!/bin/bash

#SBATCH --job-name=compilation_ex2C
#SBATCH --nodes=1
#SBATCH --time=00:10:00
#SBATCH --partition=EPYC

# set -x

# Load the required modules
module load architecture/AMD
module load openMPI/4.1.5/gnu/12.2.1
module load cmake/3.28.1

# Build the project
cmake -S . -B build/
make -C build/

# set +x

# Check if the compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful :)"
else
    echo "Compilation failed :("
    exit 1
fi
