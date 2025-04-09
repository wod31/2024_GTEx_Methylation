#!/bin/bash
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Get donor genotypes of methylation QTLs (mQTLs) to run in a cluster

#SBATCH --job-name=imQTL_Donor
#SBATCH --output=/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/out/imQTL_Donor_p_t.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/err/imQTL_Donor_p_t.err
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH --time=48:00:00

#/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Scripts/08.imQTL.Donor_Genotypes.sh

vcf_file="GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001.vcf.gz"
pathIn=/gpfs/projects/bsc83/Projects/ribosomal_proteins/Raquel/00_data/genotyped_variants/

# Variables ####
pathOut=/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/Fst/imQTL_GT/
mkdir -p ${pathOut}

module load htslib #/1.10.2 
module load vcftools

#tissue=Prostate
tissue_id=`echo $line | awk '{print $3}'`
vcftools --gzvcf ${pathIn}/${vcf_file} --snps <(zcat /gpfs/scratch/bsc83/bsc83535/GTEx/v9/mQTLs/${tissue}.mQTLs.conditional.txt.gz  | cut -f2 | awk '{if(NR>1)print}')  --recode --stdout | gzip -c > ${pathOut}/${tissue}.imQTLs.Donor_genotypes.vcf.gz
vcf-query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%SAMPLE=%GTR]\n' ${pathOut}/${tissue}.imQTLs.Donor_genotypes.vcf.gz | gzip -c > ${pathOut}/${tissue}.imQTLs.Donor_genotypes.txt.gz 

