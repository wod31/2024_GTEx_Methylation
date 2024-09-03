#!/bin/bash

vcf_file="GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001.vcf.gz"
pathIn=/gpfs/projects/bsc83/Projects/ribosomal_proteins/Raquel/00_data/genotyped_variants/

# Variables ####
pathOut=/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/Fst/imQTL_GT/
mkdir -p ${pathOut}

module load htslib #/1.10.2 
module load vcftools
while read line
do
	tissue=`echo $line | awk '{print $1}'`
	tissue_id=`echo $line | awk '{print $3}'`
	#vcftools --gzvcf ${pathIn}/${vcf_file} --snps <(zcat /gpfs/scratch/bsc83/bsc83535/GTEx/v9/mQTLs/${tissue}.mQTLs.conditional.txt.gz  | cut -f2 | awk '{if(NR>1)print}')  --recode --stdout | gzip -c > ${pathOut}/${tissue}.imQTLs.Donor_genotypes.vcf.gz
	vcf-query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%SAMPLE=%GTR]\n' ${pathOut}/${tissue}.imQTLs.Donor_genotypes.vcf.gz | gzip -c > ${pathOut}/${tissue}.imQTLs.Donor_genotypes.txt.gz 
done < Scripts/Tissue_names.tab
