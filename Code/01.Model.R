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
#### only run with categorical ancestry
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

metadata$SEX <- as.factor(metadata$SEX)
if(length(levels(metadata$SEX))==1){
  print("Sexual tissue")
  metadata <- metadata[,-which(names(metadata) == "SEX")]
  individual_variables <- c("EURv1", "AGE", "BMI", "TRISCHD", "DTHHRDY")
  # individual_variables <- c("EURv1", "AGE", "TRISCHD", "DTHHRDY")
} else{
  individual_variables <- c("EURv1", "SEX", "AGE", "BMI", "TRISCHD", "DTHHRDY")
  # individual_variables <- c("EURv1", "SEX", "AGE", "TRISCHD", "DTHHRDY")
}
metadata$DTHHRDY <- as.factor(metadata$DTHHRDY)
rownames(metadata) <- metadata$SUBJID
### filter ovary samples ####
#metadata <- metadata[!metadata$SUBJID %in% c('GTEX-S4UY','GTEX-15FZZ','GTEX-14LZ3','GTEX-1JMOU','GTEX-R55G','GTEX-13QJC'),]
metadata$SUBJID <- NULL
#### only run with categorical ancestry
# if(sum(metadata$Ancestry=="AFR")<5){
#   print("Not enough AFRs")
#   metadata <- metadata[, !colnames(metadata) %in% c("Ancestry")]
#   individual_variables <- individual_variables[!individual_variables %in% "Ancestry"]
# }

### make sure order is the same
beta <- beta[,rownames(metadata)]

probes <- rownames(beta)

metadata_2 <- metadata[,c("PEER1", "PEER2", "PEER3", "PEER4", "PEER5","PEER18","PEER9","PEER8", individual_variables)]

#Ading Chip, Array and Plate? No, PEER values are very correlated to them
# technical <- read.csv(paste0(first_dir, "/scratch/bsc83/bsc83535/GTEx/v9/Oliva/MethylationEPIC_GTEX_Batch01_11_SampleSheet.15-12-2021.1000_samples.csv"))
# technical <- technical[,c("Sample.ID.for.data.sharing.and.public.release", "Sample_Plate", "Sentrix_ID", "Sentrix_Position")]
# colnames(technical) <- c("Sample_ID", "Plate", "Chip", "Array")
# technical$Plate <- as.factor(make.names(technical$Plate))
# technical$Chip <- as.factor(make.names(technical$Chip))
# technical$Array <- as.factor(technical$Array)
# #We want the samples from "technical" that belong to the tissue of interest
# tissues <- list.dirs(paste0(project_path, "Tissues"), full.names = F)[-1]
# tissue_info <- cbind(tissues, c("Breast - Mammary Tissue", "Colon - Transverse", "Kidney - Cortex",        
#                                 "Lung", "Muscle - Skeletal", "Ovary", "Prostate", "Testis", "Whole Blood"))
# sample_tissue <- read.delim(paste0(first_dir, "/scratch/bsc83/bsc83535/GTEx/v9/Oliva/eGTExDNA_Pierce_Jan18.09-11-2021.tsv"))
# ids <- sample_tissue$Sample.ID.for.data.sharing.and.public.release[sample_tissue$Tissue.Site.Detail==tissue_info[tissue_info[,1]==tissue, 2]]
# technical_subset <- technical[technical$Sample_ID %in% ids,]
# technical_subset$Sample_ID <- as.character(technical_subset$Sample_ID)
# technical_subset$SUBJID <- sapply(technical_subset$Sample_ID, function(i) paste(unlist(strsplit(i,split = "-"))[1:2],collapse="-" ) )
# technical_subset$Sample_ID <- NULL

# metadata_2 <- merge(metadata_2, technical_subset, by.x="row.names", by.y="SUBJID")
# rownames(metadata_2) <- metadata_2$Row.names
# metadata_2$Row.names <- NULL

#Dropping levels that our subset do not contain:
# metadata_2$Plate <- droplevels(metadata_2$Plate)
# metadata_2$DTHHRDY <- droplevels(metadata_2$DTHHRDY)
# metadata_2$Chip <- droplevels(metadata_2$Chip)

print("metadata is prepared")

library(limma)

limma_function <- function(fit, x){
  covariate <<- x #makeContrast does not read the function's environment, so I add covariate to the general environment in my session
  contrast.matrix <- suppressWarnings(makeContrasts(covariate, levels=fit$design)) #Warnings due to change of name from (Intercept) to Intercept
  fitConstrasts <- suppressWarnings(contrasts.fit(fit, contrast.matrix)) #Warning due to Intercept name
  eb = eBayes(fitConstrasts)
  tt.smart.sv <- topTable(eb,adjust.method = "BH",number=Inf)
  return(tt.smart.sv)
}

beta <- sapply(beta, as.numeric)
rownames(beta) <- probes
M <- log2(beta/(1-beta)) #-> most papers say M is better for differential although beta should be plotted for interpretation

mod_2 <- model.matrix( as.formula(paste0("~", paste0(colnames(metadata_2), collapse="+"))), data =  metadata_2)

model_function <- function(mod){
  print("Modelling")
  Sys.time()
  fit_M <- lmFit(M, mod)
  # summary(decideTests(fit_M))
  # to_run <- c("EURv1", "SEX2", "AGE")
  to_run <- c("EURv1", "SEX2", "AGE", "BMI")
  ### run with categorical Ancestry 
  # if(length(levels(metadata$Ancestry))>2){ #Do we have Ancestry AMR?
  #   to_run <- c("AncestryAMR", "AncestryEUR", "SEX2", "AGE", "BMI", "AncestryAMR-AncestryEUR")
  # } else if(is.null(metadata$Ancestry)){
  #   to_run <- c("SEX2", "AGE", "BMI")
  # } else{
  #   to_run <- c("AncestryEUR", "SEX2", "AGE", "BMI")
  # }
  if(!"SEX" %in% colnames(metadata)){
    to_run <- to_run[!to_run=="SEX2"]
  }
  print(to_run)
  res_M <- lapply(to_run, function(x) limma_function(fit_M, x) )
  names(res_M) <- to_run

  to_return <- res_M
  Sys.time()
  return(to_return)
}

res_2 <- model_function(mod_2) #I would use 5 PEERs
# saveRDS(res_2, paste0(project_path, "/Tissues/", tissue, "/DML_results.rds")) #Jose used the wrong name at the beginning, DMR instead of DML
# saveRDS(res_2, paste0(project_path, "/Tissues/", tissue, "/DML_results_batch_5_PEERs.rds"))
saveRDS(res_2, paste0(project_path, "/Tissues/", tissue, "/DML_results_5_PEERs_continous.peer.rds")) #This is the final model we are using
# saveRDS(res_2, paste0(project_path, "/Tissues/", tissue, "/DML_results_no_BMI.rds")) 
print("Using 5 PEERs:")
print(paste0("  EURv1: ", sum(res_2$EURv1$adj.P.Val<0.05)))
print(paste0("  Age: ", sum(res_2$AGE$adj.P.Val<0.05)))
# print(paste0("  BMI: ", sum(res_2$BMI$adj.P.Val<0.05)))

