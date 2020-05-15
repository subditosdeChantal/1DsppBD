#!/bin/bash

#SBATCH --job-name=L2000
#SBATCH --output=jobs/L2000_%A_%a.out
#SBATCH --error=jobs/L2000_%A_%a.err
#SBATCH --array=0-3
#SBATCH --partition=normal
#SBATCH --ntasks=1

# Print the task id.
b=$(seq -f '%02g' $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_ID)
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "My PADDED_SLURM_ARRAY_TASK_ID: " $b

# Add lines here to run your computations.
jobs/BDp_02_L2000 < "jobs/L2000_"$b".input"
