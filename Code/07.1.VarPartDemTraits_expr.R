###### Variance partition for demographical traits ######
#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Run variance partition analysis on gene expression data per tissue on demographic traits
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
data <- readRDS(paste0(project_path,'Tissues/', tissue, "/counts.rds")) 

print("Reading metadata")
metadata <- readRDS(paste0(project_path, 'Tissues/', tissue, "/metadata_expression.rds"))
metadata_m <- readRDS(paste0(project_path, 'Tissues/', tissue, "/metadata.rds"))

print("Reading Admixture results")
admixture_ancestry <- read.table('Data/admixture_inferred_ancestry.txt')
colnames(admixture_ancestry) <- c('Donor','AFRv1','EURv1','inferred_ancestry','AFRv2','EURv2')
metadata <- merge(metadata, admixture_ancestry[,c("Donor","EURv1")], by='Donor')

rownames(metadata) <- metadata$Donor
metadata <- metadata[rownames(metadata) %in% metadata_m$SUBJID,]

metadata$Sex <- as.factor(metadata$Sex)
if(length(levels(metadata$Sex))==1){
  print("Sexual tissue")
  metadata <- metadata[,-which(names(metadata) == "Sex")]
  individual_variables <- c("EURv1", "Age", "BMI", "IschemicTime", "HardyScale")

} else{
  individual_variables <- c("EURv1", "Sex", "Age", "BMI", "IschemicTime", "HardyScale")

}

metadata$HardyScale <- droplevels(metadata$HardyScale)
metadata$IschemicTime <- as.numeric(metadata$IschemicTime)

### performing linear model #### 

colnames(data) <- sapply(colnames(data), function(id) paste0(strsplit(id, "-")[[1]][-3], collapse="-"))
data <- data[,rownames(metadata)]
print(ncol(metadata))
print(nrow(data))

print('genes prepared')

### run model ####
library('variancePartition')
library('BiocParallel')
# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param = SnowParam(40, "SOCK", progressbar=TRUE)
register(param)

print('create model')

metadata_2 <- metadata[,c("PEER1", "PEER2","RIN","ExonicRate", individual_variables)]

form <- paste0("~", paste0(colnames(metadata_2), collapse="+"))
print(form)

print('variance Partition')

print(paste0("Computing var part for ", tissue))

print(paste0("Computing var part for nº", nrow(data)))
print(dim(data))


vp = fitExtractVarPartModel(data, form, metadata_2)
pl <- plotVarPart( sortCols(vp))

saveRDS(vp, paste0(tissue,'_var_part_expression.rds'))

# pdf(file = paste0("Plots/",tissue,"_var_part.pdf"), w = 6, h = 3.5)
# print(pl)
# dev.off()