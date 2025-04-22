#!/bin/bash
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Script to run on the cluster the modeling of cis-driven effects

#SBATCH --job-name=cis-driven-DMP
#SBATCH --output=/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/out/CisDriven_%A.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/err/CisDriven_%A.err
#SBATCH --cpus-per-task=48
#SBATCH --ntasks=2
#SBATCH --qos=gp_bscls
#SBATCH --account=bsc83
#SBATCH --constraint=highmem
#SBATCH --time=48:00:00

module load icu R

echo $1
Rscript Scripts/08.Model.cis_driven.R -t $1
