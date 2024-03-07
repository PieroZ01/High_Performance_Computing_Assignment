#!/bin/bash

#SBATCH --job-name=bcast_default
#SBATCH --nodes=2
#SBATCH --partition=EPYC
#SBATCH --exclusive
#SBATCH --time=02:00:00

# Load the required MPI module
module load openMPI/4.1.5/gnu/12.2.1

# Navigate to the OSU Micro-Benchmarks directory
cd ./Compile_OSU/osu-micro-benchmarks-7.3/c/mpi/collective/blocking/

# Define the output file path
output_csv="../../bcast/Results/bcast_default.csv"

# Define the range of cores values
n_cores=$(seq 2 4 256)
# Define the number of iterations
iter=1000
# Define the map types
map="core socket node"

# Go to the csv output file and write the header
echo "Algorithm,Mapping,Processes,MessageSize,Latency" > $output_csv

# Run the benchmark test selecting the default algorithm 0
for map in $map
do
    for cores in $n_cores
    do
        echo "----------------------------------------------------------------------------------------------------------------------------------"
        echo "Benchmarking Bcast with $cores processes and $map mapping"
        mpirun -np $cores --map-by $map \
        --mca coll_tuned_use_dynamic_rules true \
        --mca coll_tuned_bcast_algorithm 0 \
        osu_bcast -i $iter -f -z \
        | tail -n 21 | awk -v np="$np" -v map="$map" '{printf "Default,%s,%s,%s,%s\n",map,np,$1, $2}' \
        | sed 's/,$//' >> $output_csv
    done
done

# Print the completion message
echo "----------------------------------------------------------------------------------------------------------------------------------"
echo "Job completed!"
echo "----------------------------------------------------------------------------------------------------------------------------------"
echo "Run on partition: $SLURM_JOB_PARTITION"
echo "Host: $SLURM_JOB_NODELIST"
echo "----------------------------------------------------------------------------------------------------------------------------------"
