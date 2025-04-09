###### Variance partition for demographical traits ######
#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Run variance partition analysis on DNA methylation data per tissue on demographic traits
# @software version: R=4.2.2


first_dir <- "/gpfs/"
setwd(paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Methylation/"))

# -------------- #
print(Sys.time())
#-------------- #

# Parsing
library(optparse)
parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-t", "--tissue"), type="character",
                     dest="tissue",
                     help="Tissue")
options=parse_args(parser)
tissue=options$tissue
print(tissue)

first_dir <- "/gpfs/"
project_path <- paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Methylation/")

print("Reading data")
Sys.time()
data <- readRDS(paste0(project_path,'Tissues/', tissue, "/data.rds")) 
Sys.time() #12 minutes to load 15 Gb
beta <- data

print("Reading metadata")
metadata <- readRDS(paste0(project_path, 'Tissues/', tissue, "/metadata.rds"))

print("Reading Admixture results")
admixture_ancestry <- read.table('/gpfs/scratch/bsc83/bsc83535/GTEx/v8/genotype_data/admixture_inferred_ancestry.txt')
colnames(admixture_ancestry) <- c('SUBJID','AFRv1','EURv1','inferred_ancestry','AFRv2','EURv2')
metadata <- merge(metadata, admixture_ancestry[,c("SUBJID","EURv1")], by='SUBJID')

metadata$SEX <- as.factor(metadata$SEX)
if(length(levels(metadata$SEX))==1){
  print("Sexual tissue")
  metadata <- metadata[,-which(names(metadata) == "SEX")]
   individual_variables <- c("EURv1", "AGE", "BMI", "TRISCHD", "DTHHRDY")
 
} else{
   individual_variables <- c("EURv1", "SEX", "AGE", "BMI", "TRISCHD", "DTHHRDY")
  
}
metadata$DTHHRDY <- as.factor(metadata$DTHHRDY)
rownames(metadata) <- metadata$SUBJID
metadata$SUBJID <- NULL

### performing linear model #### 

data <- data[,rownames(metadata)]
probes <- rownames(data)

# Create M values
data <- sapply(data, as.numeric)
rownames(data) <- probes
M <- log2(data/(1-data)) 

print('Betas prepared')

### run model ####
library('variancePartition')
library('BiocParallel')
# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param = SnowParam(40, "SOCK", progressbar=TRUE)
register(param)

print('create model')
if (tissue == 'Testis') {
  metadata_2 <- metadata[,c("PEER1", "PEER2", "PEER3", "PEER4", individual_variables)]
} else {
  metadata_2 <- metadata[,c("PEER1", "PEER2", "PEER3", "PEER4", "PEER5", individual_variables)]
}

form <- paste0("~", paste0(colnames(metadata_2), collapse="+"))
print(form)

print('variance Partition')

### select only DMPs
print(paste0("Computing var part for ", tissue))

dma_res <- readRDS(paste0("Tissues/", tissue, "/DML_results_5_PEERs_continous.rds"))
traits <- names(dma_res)
dm_cpgs <- unique(unlist(lapply(traits, function(trait) rownames(dma_res[[trait]][dma_res[[trait]][,"adj.P.Val"] < 0.05,]) ) ))

print(paste0("Computing var part for nº", length(dm_cpgs)))
M <- M[dm_cpgs,]
print(dim(M))

vp = fitExtractVarPartModel(M, form, metadata_2)
pl <- plotVarPart( sortCols(vp))

saveRDS(vp, paste0(tissue,'_var_part.rds'))

pdf(file = paste0("Plots/",tissue,"_var_part.pdf"), w = 6, h = 3.5)
print(pl)
dev.off()