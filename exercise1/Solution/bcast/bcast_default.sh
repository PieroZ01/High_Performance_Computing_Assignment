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
output_csv="../../../../../../bcast/Results/bcast_default.csv"

# Define the range of cores values
step=4
min_cores=2
max_cores=256
# Define the number of iterations
iter=3000
# Define the map types
maps="core socket"

# Go to the csv output file and write the header
echo "Algorithm,Mapping,Processes,MessageSize,Latency" > $output_csv

# Run the benchmark test selecting the default algorithm 0
for map in $maps
do
    for cores in $(seq $min_cores $step $max_cores)
    do
        echo "----------------------------------------------------------------------------------------------------------------------------------"
        echo "Benchmarking Bcast with $cores processes and $map mapping"
        mpirun -np $cores --map-by $map \
        --mca coll_tuned_use_dynamic_rules true \
        --mca coll_tuned_bcast_algorithm 0 \
        osu_bcast -x 1000 -i $iter -f -z \
        | tail -n 21 | awk -v cores="$cores" -v map="$map" '{printf "Default,%s,%s,%s,%s\n",map,cores,$1,$2}' \
        | sed 's/,$//' >> $output_csv
    done
    # Run the test with 256 cores as well
    cores_final=256
    echo "----------------------------------------------------------------------------------------------------------------------------------"
    echo "Benchmarking Bcast with $cores_final processes and $map mapping"
    mpirun -np $cores_final --map-by $map \
    --mca coll_tuned_use_dynamic_rules true \
    --mca coll_tuned_bcast_algorithm 0 \
    osu_bcast -x 1000 -i $iter -f -z \
    | tail -n 21 | awk -v cores="$cores_final" -v map="$map" '{printf "Default,%s,%s,%s,%s\n",map,cores,$1,$2}' \
    | sed 's/,$//' >> $output_csv
done

# Print the completion message
echo "----------------------------------------------------------------------------------------------------------------------------------"
echo "Job completed!"
echo "----------------------------------------------------------------------------------------------------------------------------------"
echo "Run on partition: $SLURM_JOB_PARTITION"
echo "Host: $SLURM_JOB_NODELIST"
echo "----------------------------------------------------------------------------------------------------------------------------------"
