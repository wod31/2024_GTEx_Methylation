#!/bin/bash

#SBATCH --job-name=plink_filter
#SBATCH --output=/gpfs/scratch/bsc83/bsc83535/GTEx/v8/genotype_data/plink_filter2.out
#SBATCH --error=/gpfs/scratch/bsc83/bsc83535/GTEx/v8/genotype_data/plink_filter2.err
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=2
#SBATCH --time=05:00:00

module load plink


#plink --bfile GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased --indep-pairwise 50 10 0.1

#plink --bfile GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased --extract plink.prune.in --make-bed --out GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phasedLDPrune

#plink --bfile GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased --chr 21 --make-bed --out GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phasedLDPrunechr21