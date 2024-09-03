###################################################
#### get residuals methylation for each trait #####
###################################################

#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez and Winona Oliveros
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to obtain residuals for methylation
# @software version: R=4.2.2
rm(list=ls())

# Load libraries ####
suppressMessages(library(edgeR))
library(ggplot2)
library(ggpubr)

#Set path 
# setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")
first_dir <- "/gpfs/"
setwd(paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Methylation/"))

# -------------- #
print(Sys.time())
#-------------- #

tissues <- list.dirs("Tissues/", full.names = F)[-1]
# tissues <- tissues[tissues!="KidneyCortex"]
# tissues <- tissues[33:length(tissues)]
tissues <- 'Ovary'

sexual_tissues <- c("Prostate", "Testis", "Ovary")
#Do a for loop 
for(i in 1:length(tissues)){
  tissue <- tissues[i]
  print(paste0("Computing residuals for ", tissue))
  metadata <- readRDS(paste0("Tissues/", tissue, "/metadata.rds"))
  
  # if(is.null(exprs_genes)){
  #   exprs_genes <- rownames(dea_res$BMI)
  # }
  # Counts
  beta <- readRDS(paste0("Tissues/", tissue, "/data.rds")) #Beta value of CpGs assayed
  
  print("Reading Admixture results")
  admixture_ancestry <- read.table('/gpfs/scratch/bsc83/bsc83535/GTEx/v8/genotype_data/admixture_inferred_ancestry.txt')
  colnames(admixture_ancestry) <- c('SUBJID','AFRv1','EURv1','inferred_ancestry','AFRv2','EURv2')
  metadata <- merge(metadata, admixture_ancestry[,c("SUBJID","EURv1")], by='SUBJID')
  
  ### Filter metadata and betas
  # metadata$Ancestry <- as.factor(metadata$Ancestry)
  # if(sum(metadata$Ancestry=="AMR")<5){
  #   print("Not enough AMRs")
  #   metadata <- metadata[metadata$Ancestry!="AMR",]
  #   beta <- beta[,metadata$SUBJID]
  # }
  # 
  # print("Filtering ASN")
  # metadata <- metadata[metadata$Ancestry!="ASN",]
  # metadata$Ancestry <- droplevels(metadata$Ancestry)
  # beta <- beta[,colnames(beta) %in% metadata$SUBJID]
  # 
  metadata$SEX <- as.factor(metadata$SEX)
  if(tissue %in% sexual_tissues){
    print("Sexual tissue")
    metadata <- metadata[,-which(names(metadata) == "SEX")]
    individual_variables <- c("EURv1", "AGE", "BMI", "TRISCHD", "DTHHRDY")
    traits <- c("EURv1", "AGE", "BMI")
  } else{
    individual_variables <- c("EURv1", "SEX", "AGE", "BMI", "TRISCHD", "DTHHRDY")
    traits <- c("EURv1", "SEX", "AGE", "BMI")
  }
  metadata$DTHHRDY <- as.factor(metadata$DTHHRDY)
  rownames(metadata) <- metadata$SUBJID
  metadata$SUBJID <- NULL
  # if(sum(metadata$Ancestry=="AFR")<5){
  #   print("Not enough AFRs")
  #   metadata <- metadata[, !colnames(metadata) %in% c("Ancestry")]
  #   individual_variables <- individual_variables[!individual_variables %in% "Ancestry"]
  # }
  
  ### mke sure order is the same
  beta <- beta[,rownames(metadata)]
  probes <- rownames(beta)
  metadata_2 <- metadata[,c("PEER1", "PEER2", "PEER3", "PEER4", "PEER5", individual_variables)]
  metadata_2 <- metadata[,c(individual_variables)]
  
  # #Ading Chip, Array and Plate 
  # technical <- read.csv(paste0(first_dir, "/scratch/bsc83/bsc83535/GTEx/v9/Oliva/MethylationEPIC_GTEX_Batch01_11_SampleSheet.15-12-2021.1000_samples.csv"))
  # technical <- technical[,c("Sample.ID.for.data.sharing.and.public.release", "Sample_Plate", "Sentrix_ID", "Sentrix_Position")]
  # colnames(technical) <- c("Sample_ID", "Plate", "Chip", "Array")
  # technical$Plate <- as.factor(make.names(technical$Plate))
  # technical$Chip <- as.factor(make.names(technical$Chip))
  # technical$Array <- as.factor(technical$Array)
  # #We want the samples from "technical" that belong to the tissue of interest
  # tissues <- list.dirs("Tissues/", full.names = F)[-1]
  # tissue_info <- cbind(tissues, c("Breast - Mammary Tissue", "Colon - Transverse", "Kidney - Cortex",        
  #                                 "Lung", "Muscle - Skeletal", "Ovary", "Prostate", "Testis", "Whole Blood"))
  # sample_tissue <- read.delim(paste0(first_dir, "/scratch/bsc83/bsc83535/GTEx/v9/Oliva/eGTExDNA_Pierce_Jan18.09-11-2021.tsv"))
  # ids <- sample_tissue$Sample.ID.for.data.sharing.and.public.release[sample_tissue$Tissue.Site.Detail==tissue_info[tissue_info[,1]==tissue, 2]]
  # technical_subset <- technical[technical$Sample_ID %in% ids,]
  # technical_subset$Sample_ID <- as.character(technical_subset$Sample_ID)
  # technical_subset$SUBJID <- sapply(technical_subset$Sample_ID, function(i) paste(unlist(strsplit(i,split = "-"))[1:2],collapse="-" ) )
  # technical_subset$Sample_ID <- NULL
  # 
  # metadata_2 <- merge(metadata_2, technical_subset, by.x="row.names", by.y="SUBJID")
  # rownames(metadata_2) <- metadata_2$Row.names
  # metadata_2$Row.names <- NULL
  
  #Dropping levels that our subset do not contain:
  #metadata_2$Plate <- droplevels(metadata_2$Plate)
  metadata_2$DTHHRDY <- droplevels(metadata_2$DTHHRDY)
  # metadata_2$Chip <- droplevels(metadata_2$Chip)
  # metadata_2$Array <- droplevels(metadata_2$Array)
  
  print("metadata is prepared")
  
  # Create M values
  beta <- sapply(beta, as.numeric)
  rownames(beta) <- probes
  M <- log2(beta/(1-beta)) #-> most papers say M is better for differential although beta should be plotted for interpretation
  
  #mod_2 <- model.matrix( as.formula(paste0("~", paste0(colnames(metadata_2), collapse="+"))), data =  metadata_2)
  
  for (trait in traits) {
    covariates <- names(metadata_2)[!names(metadata_2) %in% c("Donor", "Sample", trait,"Chip")]
    
    # Expression ~ covariates 
    fml_args_mod <- paste(c(covariates), collapse = " + ")
    mod <- model.matrix( as.formula(paste(" ~  ", paste(fml_args_mod,collapse = " "))), data =  metadata_2)
    
    # Limma fit
    fit <- lmFit(M, mod)
    
    # Calculate expression residuals
    exprs_residuals <- resid(fit, M)
    
    saveRDS(exprs_residuals, paste0("Tissues/", tissue, "/",trait,"_methylation_residuals.continous.noPEERs.rds"))
  }
}
  

