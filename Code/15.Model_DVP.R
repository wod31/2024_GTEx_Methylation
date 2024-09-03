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
# first_dir <- "~/marenostrum/"
# first_dir <- "~/Documents/mn4/"
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
admixture_ancestry <- read.table(paste0(first_dir, '/scratch/bsc83/MN4/bsc83/bsc83535/GTEx/v8/genotype_data/admixture_inferred_ancestry.txt'))
colnames(admixture_ancestry) <- c('SUBJID','AFRv1','EURv1','inferred_ancestry','AFRv2','EURv2')
metadata <- merge(metadata, admixture_ancestry[,c("SUBJID","EURv1")], by='SUBJID')

metadata$SEX <- as.factor(metadata$SEX)
if(length(levels(metadata$SEX))==1){
  sexual <- T
  print("Sexual tissue")
  metadata <- metadata[,-which(names(metadata) == "SEX")]
  individual_variables <- c("EURv1", "AGE", "BMI", "TRISCHD", "DTHHRDY")
} else{
  sexual <- F
  individual_variables <- c("EURv1", "SEX", "AGE", "BMI", "TRISCHD", "DTHHRDY")
}
metadata$DTHHRDY <- as.factor(metadata$DTHHRDY)
rownames(metadata) <- metadata$SUBJID
metadata$SUBJID <- NULL

### make sure order is the same
beta <- beta[,rownames(metadata)]

probes <- rownames(beta)

metadata_2 <- metadata[,c("PEER1", "PEER2", "PEER3", "PEER4", "PEER5", individual_variables)]

print("metadata is prepared")

beta <- sapply(beta, as.numeric)
rownames(beta) <- probes
M <- log2(beta/(1-beta)) 

mod_2 <- model.matrix( as.formula(paste0("~", paste0(colnames(metadata_2), collapse="+"))), data =  metadata_2)


# #DVPs:
library(missMethyl)
if(sexual){
  fitvar <- varFit(M, design = mod_2, coef = c(1, 7, 8, 9))
} else{
  fitvar <- varFit(M, design = mod_2, coef = c(1, 7, 8, 9, 10))
}
library(limma)
summary(decideTests(fitvar))

topDV <- topVar(fitvar, coef=7, number=150000)
#topDV <- topDV[topDV$Adj.P.Value<0.05,]
saveRDS(topDV, paste0(project_path, "/Tissues/", tissue, "/DVP_Ancestry.rds"))
# 
# if(sexual){
#   topDV <- topVar(fitvar, coef=8, number=10000)
#   topDV <- topDV[topDV$Adj.P.Value<0.05,]
#   saveRDS(topDV, paste0(project_path, "/Tissues/", tissue, "/DVP_Age.rds"))
#   
#   topDV <- topVar(fitvar, coef=9, number=10000)
#   topDV <- topDV[topDV$Adj.P.Value<0.05,]
#   saveRDS(topDV, paste0(project_path, "/Tissues/", tissue, "/DVP_BMI.rds"))
# } else{
#   topDV <- topVar(fitvar, coef=8, number=10000)
#   topDV <- topDV[topDV$Adj.P.Value<0.05,]
#   saveRDS(topDV, paste0(project_path, "/Tissues/", tissue, "/DVP_Sex.rds"))
#   
#   topDV <- topVar(fitvar, coef=9, number=10000)
#   topDV <- topDV[topDV$Adj.P.Value<0.05,]
#   saveRDS(topDV, paste0(project_path, "/Tissues/", tissue, "/DVP_Age.rds"))
#   
#   topDV <- topVar(fitvar, coef=10, number=10000)
#   topDV <- topDV[topDV$Adj.P.Value<0.05,]
#   saveRDS(topDV, paste0(project_path, "/Tissues/", tissue, "/DVP_BMI.rds"))
# }
