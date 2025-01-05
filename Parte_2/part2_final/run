#!/bin/sh
#
#SBATCH --exclusive     # exclusive node for the job
#SBATCH --time=05:00    # allocation for 5 minutes

for threads in {8..33}
do
    export OMP_NUM_THREADS=$threads
    echo "Running with OMP_NUM_THREADS=$threads"
    perf stat -r 5 -e instructions,cycles ./fluid_sim 
done