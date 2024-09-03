###################################################
#### get residuals expressiom for each trait #####
###################################################

#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez and Raquel Garcia-Perez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to run residuals on gene expression
# @software version: R=4.2.2
rm(list=ls())

# Load libraries ####
suppressMessages(library(edgeR))
library(ggplot2)
library(ggpubr)

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
  metadata <- readRDS(paste0("Tissues/", tissue, "/metadata_expression.rds"))

  print("Reading Admixture results")
  admixture_ancestry <- read.table('/gpfs/scratch/bsc83/bsc83535/GTEx/v8/genotype_data/admixture_inferred_ancestry.txt')
  colnames(admixture_ancestry) <- c('Donor','AFRv1','EURv1','inferred_ancestry','AFRv2','EURv2')
  metadata <- merge(metadata, admixture_ancestry[,c("Donor","EURv1")], by='Donor')
  
  metadata$Sex <- as.factor(metadata$Sex)
  if(tissue %in% sexual_tissues){
    print("Sexual tissue")
    metadata <- metadata[,-which(names(metadata) == "Sex")]
    individual_variables <- c("EURv1", "Age", "BMI", "IschemicTime", "HardyScale")
    traits <- c("EURv1", "Age", "BMI")
  } else{
    individual_variables <- c("EURv1", "Sex", "Age", "BMI", "IschemicTime", "HardyScale")
    traits <- c("EURv1", "Sex", "Age", "BMI")
  }
  metadata$HardyScale <- as.factor(metadata$HardyScale)
  rownames(metadata) <- metadata$Donor
  metadata$Donor <- NULL
  
  # Counts
  counts <- readRDS(paste0("Tissues/", tissue, "/counts.rds")) #Counts of expressed genes in samples of interest for the given tissue
  
  # Create DGEList object
  dge <- DGEList(counts)
  
  # Calculate normalization factors (does not do the normalization, only computes the factors)
  dge <- calcNormFactors(dge)
  
  # Voom
  v <- voom(dge, design = NULL, normalize="quantile", save.plot=F, plot = F) # samples are treated as replicates
  
  rownames(metadata) <- metadata$Sample
  metadata_2 <- metadata[,c("PEER1", "PEER2","RIN","ExonicRate", individual_variables)]
  metadata_2 <- metadata[,c("RIN","ExonicRate", individual_variables)]
  metadata_2$HardyScale <- droplevels(metadata_2$HardyScale)
  metadata_2$IschemicTime <- as.numeric(metadata_2$IschemicTime)

  for (trait in traits) {
    covariates <- names(metadata_2)[!names(metadata_2) %in% c("Donor", "Sample", trait)]
    
    # Expression ~ covariates 
    fml_args_mod <- paste(c(covariates), collapse = " + ")
    mod <- model.matrix( as.formula(paste(" ~  ", paste(fml_args_mod,collapse = " "))), data =  metadata_2)
    
    # Limma fit
    fit <- lmFit(v, mod)
    
    # Calculate expression residuals
    exprs_residuals <- resid(fit, v)
    
    saveRDS(exprs_residuals, paste0("Tissues/", tissue, "/",trait,"_expression_residuals.continous.noPEERs.rds"))
  }
  
  covariates <- names(metadata_2)[!names(metadata_2) %in% c("Donor", "Sample")]

  Expression ~ covariates
  fml_args_mod <- paste(c(covariates), collapse = " + ")
  mod <- model.matrix( as.formula(paste(" ~  ", paste(fml_args_mod,collapse = " "))), data =  metadata_2)

  # Limma fit
  fit <- lmFit(v, mod)
  #library(limma)

  limma_function <- function(fit, x){
    covariate <<- x #makeContrast does not read the function's environment, so I add covariate to the general environment in my session
    contrast.matrix <- suppressWarnings(makeContrasts(covariate, levels=fit$design)) #Warnings due to change of name from (Intercept) to Intercept
    fitConstrasts <- suppressWarnings(contrasts.fit(fit, contrast.matrix)) #Warning due to Intercept name
    eb = eBayes(fitConstrasts)
    tt.smart.sv <- topTable(eb,adjust.method = "BH",number=Inf)
    return(tt.smart.sv)
  }

  to_run <- 'Age'
  print(to_run)
  res_M <- lapply(to_run, function(x) limma_function(fit, x) )
  names(res_M) <- to_run

  saveRDS(res_M, paste0("Tissues/", tissue, "/",to_run,"_differential_expression.continous.noPEERs.rds"))

}
