#!/bin/bash

#SBATCH --job-name=downsampling_array
#SBATCH --output=/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/out/downsampling_array_%A_%a.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/out/downsampling_array_%A_%a.err
#SBATCH --cpus-per-task=30
#SBATCH --time=05:00:00
#SBATCH --array=1-9



module load R

file=/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Downsampling/tissues.txt
export tissue=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $file | xargs)


echo $tissue

Rscript Scripts/14.9.Model_downsampling.R -t $tissue




