#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Process TFBS data
# @software version: R=4.2.2

suppressMessages(library(dplyr))
library(tidyr)
suppressMessages(library(data.table))
library(valr)

#Set path 
setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")

print("Reading CpG annotation") #3 minutes
Sys.time()
first_dir <- "/gpfs/"
# first_dir <- "~/Documents/mn4/"
data_path <- paste0(first_dir, "scratch/bsc83/bsc83535/GTEx/v9/Oliva/")
annotation <- read.csv(paste0(data_path, "GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv"))
Sys.time()

print("Preprocessing CpG annotation")

ann_bed <- annotation[!is.na(annotation$MAPINFO) & !is.na(annotation$CHR),] %>%
  dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>% distinct()
ann_bed$chrom <- paste0('chr',ann_bed$chrom)
ann_bed$start <- ann_bed$start-1


#for(tissue in c("Lng", "Myo", "Kid", "Utr", "Bld", "Brs", "Prs", "Dig", "ALL")){
#for(tissue in c("Lng", "Myo", "Kid", "Utr", "Bld", "Brs", "Prs", "Dig")){
#for(tissue in c("Myo", "Kid", "Utr", "Bld", "Brs", "Prs", "Dig", "ALL")){
for(tissue in c("ALL")){
    print(tissue) #Myo, Kid and UTR are each 5 minutes. Bld is more than 30 min, from Brs to Dig is fine, but ALL needs more memory. So, in 2 hours I could do all but ALL
  
  print("Reading chipatlas annotation")
  Sys.time()
  chip <- fread(paste0(data_path, paste0("ChipAtlas_data/Oth.", tissue, ".05.AllAg.AllCell.bed")), header = F)
  Sys.time()
  
  dim(chip)
  
  #preprocess data
  colnames(chip) <- c('chrom','start', 'end', 'to_split', 'score', 'strand', 'start_summit', 'end_summit', 'coord')
  chip$TF <- gsub('.*=','', gsub('%.*', '', chip$to_split))
  chip <- chip[,c(1:3, 10)]
  
  Sys.time()
  print("About to merge TFBS and CpGs tested")
  chip_cpgs <- bed_intersect(ann_bed, chip, suffix = c("_ann", "_tfbs"))
  chip_cpgs <- chip_cpgs[chip_cpgs$.overlap ==1,]
  
  dim(chip_cpgs)
  Sys.time()
  
  print("Removing duplicates")
  chip_cpgs <- chip_cpgs[,c(1:4,7)]
  chip_cpgs <- unique(chip_cpgs)
  dim(chip_cpgs)
  Sys.time()

  print(tissue)
  write.csv(chip_cpgs, paste0(first_dir, "projects/bsc83/Projects/GTEx_v8/Methylation/TFBS/chip_seq_", tissue, ".csv"))
  print("done")
}
