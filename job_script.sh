#!/bin/bash

#SBATCH --job-name="binara job"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --time=0-04:00:00

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "Job shell script started."
srun cmake-build-release-zaratan1/binara_exe 110602878 34 1 1 0 0