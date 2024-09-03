###### Variance partition for demographical traits ######
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")
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
data <- readRDS(paste0(project_path,'Tissues/', tissue, "/counts.rds")) #From whole compressed data in 5.6G to compressed 1.4G/1.1Gb only in Lung (the highest number of samples)
Sys.time() #12 minutes to load 15 Gb
#beta <- data

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
  # individual_variables <- c("EURv1", "AGE", "TRISCHD", "DTHHRDY")
} else{
  individual_variables <- c("EURv1", "Sex", "Age", "BMI", "IschemicTime", "HardyScale")
  # individual_variables <- c("EURv1", "SEX", "AGE", "TRISCHD", "DTHHRDY")
}

metadata$HardyScale <- droplevels(metadata$HardyScale)
metadata$IschemicTime <- as.numeric(metadata$IschemicTime)

### performing linear model #### 

#data <- data[,rownames(metadata)]
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

### select only DMPs
print(paste0("Computing var part for ", tissue))

#dma_res <- readRDS(paste0("Tissues/", tissue, "/DML_results_5_PEERs_continous.rds"))
#traits <- names(dma_res)
#dm_cpgs <- unique(unlist(lapply(traits, function(trait) rownames(dma_res[[trait]][dma_res[[trait]][,"adj.P.Val"] < 0.05,]) ) ))

#cpgs_assessed <- rownames(dma_res$AGE) #Rownames of any variable would work
print(paste0("Computing var part for nº", nrow(data)))
print(dim(data))


vp = fitExtractVarPartModel(data, form, metadata_2)
pl <- plotVarPart( sortCols(vp))

saveRDS(vp, paste0(tissue,'_var_part_expression.rds'))

# pdf(file = paste0("Plots/",tissue,"_var_part.pdf"), w = 6, h = 3.5)
# print(pl)
# dev.off()