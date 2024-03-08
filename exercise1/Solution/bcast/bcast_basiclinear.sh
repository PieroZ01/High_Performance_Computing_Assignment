#!/bin/bash

#SBATCH --job-name=bcast_basiclinear
#SBATCH --nodes=2
#SBATCH --partition=EPYC
#SBATCH --exclusive
#SBATCH --time=02:00:00

# Load the required MPI module
module load openMPI/4.1.5/gnu/12.2.1

# Navigate to the OSU Micro-Benchmarks directory
cd ./Compile_OSU/osu-micro-benchmarks-7.3/c/mpi/collective/blocking/

# Define the output file path
output_csv="../../../../../../bcast/Results/bcast_basiclinear.csv"

# Define the range of cores values
step=4
min_cores=2
max_cores=258
# Define the number of iterations
iter=1000
# Define the map types
maps="core socket node"

# Go to the csv output file and write the header
echo "Algorithm,Mapping,Processes,MessageSize,Latency" > $output_csv

# Run the benchmark test selecting the basic linear algorithm 1
for map in $maps
do
    for cores in $(seq $min_cores $step $max_cores)
    do
        echo "----------------------------------------------------------------------------------------------------------------------------------"
        echo "Benchmarking Bcast with $cores processes and $map mapping"
        mpirun -np $cores --map-by $map \
        --mca coll_tuned_use_dynamic_rules true \
        --mca coll_tuned_bcast_algorithm 1 \
        osu_bcast -i $iter -f -z \
        | tail -n 21 | awk -v cores="$cores" -v map="$map" '{printf "Default,%s,%s,%s,%s\n",map,cores,$1,$2}' \
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
