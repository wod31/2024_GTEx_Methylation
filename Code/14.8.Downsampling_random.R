#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to downsample
# @software version: R=4.2.2


#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

first_dir <- "~/Documents/mn4/"
# first_dir <- "/gpfs/"
admixture_ancestry <- read.table(paste0(first_dir, '/scratch/bsc83/bsc83535/GTEx/v8/genotype_data/admixture_inferred_ancestry.txt'))
colnames(admixture_ancestry) <- c('SUBJID','AFRv1','EURv1','inferred_ancestry','AFRv2','EURv2')


tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", 
             "KidneyCortex", "Testis", "WholeBlood", "MuscleSkeletal")

set.seed(123)
# for(number in c(30, 35, 40, 50, 100)){
for(number in c(35)){
    # print(number)
  for(tissue in tissues){
    dir.create(paste0("Downsampling/", number, "/", tissue, "/"), showWarnings = F, recursive = T)
    metadata <- readRDS(paste0("Tissues/", tissue, "/metadata.rds"))
    metadata <- merge(metadata, admixture_ancestry[,c("SUBJID","EURv1")], by='SUBJID')
    if(nrow(metadata)<number){next}
    print(tissue)
    for(i in 1:100){
      metadata_subset <- metadata[sample(1:nrow(metadata), number, replace=F),]
      saveRDS(metadata_subset, paste0("Downsampling/", number, "/", tissue, "/metadata_downsampling_", i, ".rds"))
    }
  }
}

