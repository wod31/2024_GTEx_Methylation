#!/usr/bin/env Rscript

# Parsing
library(optparse)
parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-t", "--tissue"), type="character",
                     dest="tissue",
                     help="Tissue")
options=parse_args(parser)
tissue=options$tissue
# tissue <- "Lung"

print(tissue)
#first_dir <- "~/marenostrum/"
first_dir <- "/gpfs/"
project_path <- paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Methylation/")
#project_path <- paste0(first_dir, "/Projects/GTEx_v8/Methylation/")


print("Reading data")
Sys.time()
data <- readRDS(paste0(project_path,'Tissues/', tissue, "/data.rds")) #From whole compressed data in 5.6G to compressed 1.4G/1.1Gb only in Lung (the highest number of samples)
Sys.time() #12 minutes to load 15 Gb
beta <- data

print("Reading metadata")
metadata <- readRDS(paste0(project_path, 'Tissues/', tissue, "/metadata.rds"))

print("Reading Admixture results")
admixture_ancestry <- read.table('/gpfs/scratch/bsc83/MN4/bsc83/bsc83535/GTEx/v8/genotype_data/admixture_inferred_ancestry.txt')
colnames(admixture_ancestry) <- c('SUBJID','AFRv1','EURv1','inferred_ancestry','AFRv2','EURv2')
metadata <- merge(metadata, admixture_ancestry[,c("SUBJID","EURv1")], by='SUBJID')

print("Filtering Sex")

metadata_female <- metadata[metadata$SEX == 2,-which(names(metadata) == "SEX")]
metadata_male <- metadata[metadata$SEX == 1,-which(names(metadata) == "SEX")]
individual_variables <- c("EURv1", "AGE", "BMI", "TRISCHD", "DTHHRDY")
beta_male <- beta[,colnames(beta) %in% metadata_male$SUBJID]
beta_female <- beta[,colnames(beta) %in% metadata_female$SUBJID]

metadata_female$DTHHRDY <- as.factor(metadata_female$DTHHRDY)
rownames(metadata_female) <- metadata_female$SUBJID
metadata_female$SUBJID <- NULL

metadata_male$DTHHRDY <- as.factor(metadata_male$DTHHRDY)
rownames(metadata_male) <- metadata_male$SUBJID
metadata_male$SUBJID <- NULL

### make sure order is the same
model_function <- function(mod){
  print("Modelling")
  Sys.time()
  fit_M <- lmFit(M, mod)

  to_run <- c("AGE")

  print(to_run)
  res_M <- lapply(to_run, function(x) limma_function(fit_M, x) )
  names(res_M) <- to_run
  
  to_return <- res_M
  Sys.time()
  return(to_return)
}

library(limma)
limma_function <- function(fit, x){
  covariate <<- x #makeContrast does not read the function's environment, so I add covariate to the general environment in my session
  contrast.matrix <- suppressWarnings(makeContrasts(covariate, levels=fit$design)) #Warnings due to change of name from (Intercept) to Intercept
  fitConstrasts <- suppressWarnings(contrasts.fit(fit, contrast.matrix)) #Warning due to Intercept name
  eb = eBayes(fitConstrasts)
  tt.smart.sv <- topTable(eb,adjust.method = "BH",number=Inf)
  return(tt.smart.sv)
}


beta_female <- beta_female[,rownames(metadata_female)]
beta_male <- beta_male[,rownames(metadata_male)]

probes_male <- rownames(beta_male)
probes_female <- rownames(beta_female)

metadata_female <- metadata_female[,c("PEER1", "PEER2", "PEER3", "PEER4", "PEER5", individual_variables)]
metadata_male <- metadata_male[,c("PEER1", "PEER2", "PEER3", "PEER4", "PEER5", individual_variables)]

print("metadata is prepared")

#male
beta_male <- sapply(beta_male, as.numeric)
rownames(beta_male) <- probes_male
M <- log2(beta_male/(1-beta_male)) #-> most papers say M is better for differential although beta should be plotted for interpretation

mod_2 <- model.matrix( as.formula(paste0("~", paste0(colnames(metadata_male), collapse="+"))), data =  metadata_male)

res_2 <- model_function(mod_2) #I would use 5 PEERs

saveRDS(res_2, paste0(project_path, "/Tissues/", tissue, "/DML_results_5_PEERs_continous.peer.AgeMale.rds")) 
print("Using 5 PEERs:")
print(paste0("  Age: ", sum(res_2$AGE$adj.P.Val<0.05)))

#female
beta_female <- sapply(beta_female, as.numeric)
rownames(beta_female) <- probes_female
M <- log2(beta_female/(1-beta_female)) #-> most papers say M is better for differential although beta should be plotted for interpretation

mod_2 <- model.matrix( as.formula(paste0("~", paste0(colnames(metadata_female), collapse="+"))), data =  metadata_female)

res_2 <- model_function(mod_2) #I would use 5 PEERs

saveRDS(res_2, paste0(project_path, "/Tissues/", tissue, "/DML_results_5_PEERs_continous.peer.AgeFemale.rds")) 
print("Using 5 PEERs:")
print(paste0("  Sex: ", sum(res_2$AGE$adj.P.Val<0.05)))


