#!/bin/bash

#SBATCH --job-name=barrier_linear
#SBATCH --nodes=2
#SBATCH --partition=EPYC
#SBATCH --exclusive
#SBATCH --time=02:00:00

# Load the required MPI module
module load openMPI/4.1.5/gnu/12.2.1

# Navigate to the OSU Micro-Benchmarks directory
cd ./Compile_OSU/osu-micro-benchmarks-7.3/c/mpi/collective/blocking/

# Define the output file path
output_csv="../../../../../../barrier/Results/barrier_linear.csv"

# Define the range of cores values
step=4
min_cores=2
max_cores=256
# Define the number of iterations
iter=5000
# Define the map types
maps="core socket"

# Go to the csv output file and write the header
echo "Algorithm,Mapping,Processes,Latency" > $output_csv

# Run the benchmark test selecting the linear algorithm 1
for map in $maps
do
    for cores in $(seq $min_cores $step $max_cores)
    do
        echo "----------------------------------------------------------------------------------------------------------------------------------"
        echo "Benchmarking barrier with $cores processes and $map mapping"
        mpirun -np $cores --map-by $map \
        --mca coll_tuned_use_dynamic_rules true \
        --mca coll_tuned_barrier_algorithm 1 \
        osu_barrier -x 1000 -i $iter -f -z \
        | tail -n 1 | awk -v cores="$cores" -v map="$map" '{printf "Linear,%s,%s,%s\n",map,cores,$1}' >> $output_csv
    done
    # Run the test with 256 cores as well
    cores_final=256
    echo "----------------------------------------------------------------------------------------------------------------------------------"
    echo "Benchmarking barrier with $cores_final processes and $map mapping"
    mpirun -np $cores_final --map-by $map \
    --mca coll_tuned_use_dynamic_rules true \
    --mca coll_tuned_barrier_algorithm 1 \
    osu_barrier -x 1000 -i $iter -f -z \
    | tail -n 1 | awk -v cores="$cores_final" -v map="$map" '{printf "Linear,%s,%s,%s\n",map,cores,$1}' >> $output_csv
done

# Print the completion message
echo "----------------------------------------------------------------------------------------------------------------------------------"
echo "Job completed!"
echo "----------------------------------------------------------------------------------------------------------------------------------"
echo "Run on partition: $SLURM_JOB_PARTITION"
echo "Host: $SLURM_JOB_NODELIST"
echo "----------------------------------------------------------------------------------------------------------------------------------"
