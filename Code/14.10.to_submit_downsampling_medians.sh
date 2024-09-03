#!/bin/bash

#SBATCH --job-name=downsampling_medians
#SBATCH --output=/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/out/downsampling_medians_%A.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/out/downsampling_medians_%A.err
#SBATCH --cpus-per-task=48
#SBATCH --time=01:00:00
#SBATCH --qos=debug


module load R

Rscript Scripts/14.10.Downsampling_medians.R




