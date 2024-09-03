#Creating path variables
first_dir <- "/gpfs/" #I use first_dir in case I run something locally at some point
data_path <- paste0(first_dir, "/scratch/bsc83/bsc83535/GTEx/")
project_path <- paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Methylation/")

#Reading beta values
print("Reading whole data")
Sys.time()
data <- read.table(paste0(data_path, "v9/Oliva/GSE213478_methylation_DNAm_noob_final_BMIQ_all_tissues_987.txt.gz")) #It takes one hour
Sys.time()

#Reading metadata from Oliva et al. (the covariates metadata with PEER info will miss more than 100 donors that do not contain genome data)
metadata <- read.delim(paste0(data_path, "v9/Oliva/eGTExDNA_Pierce_Jan18.09-11-2021.tsv"))
names(metadata)[1] <- "ID"
names(metadata)[15] <- "Tissue"

#Processing data format
samples <- strsplit(data[1,2], ",")[[1]]
samples <- samples[2:length(samples)]

samples_2 <- sapply(samples, function(x) substr(x, 2, nchar(x)-1))
names(samples_2) <- NULL

tissues <- sapply(samples_2, function(x) metadata$Tissue[metadata$ID==x])

#Splitting data into different tissues
print("Splitting whole data")
n <- nrow(data)
library(parallel)
Sys.time()
splitted_data <- mclapply(2:n, function(x) strsplit(data[x,2], ",")[[1]], mc.cores=8) #3000 lasts 23 seconds, with mclapply, 5, 10000 lasts 20 seconds with mclapply
Sys.time()
names(splitted_data) <- data[2:n, 1]

splitted_data_frame <- as.data.frame(do.call(rbind, splitted_data))
splitted_data_frame <- splitted_data_frame[,2:ncol(splitted_data_frame)] #to remove the first column, which is always "", maybe we could have read the data frame in another way to solve that
colnames(splitted_data_frame) <- samples_2

#Reading file with the different tissue names
translations <- readRDS(paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Jose/00_Data/Tissue_info_whole.rds"))

for(tissue in unique(tissues)){
  print(tissue)
  if(tissue=="Breast - Mammary Tissue"){ #This is not in the tissue translations file
    tissue_name <- "BreastMammaryTissue"
  } else{
    tissue_name <- translations$tissue_ID[translations$tissue_name==tissue]
  }
  dir.create(file.path(project_path, "Tissues/", tissue_name, "/"), showWarnings = FALSE) #Creating a folder per tissues
  trues <- tissues == tissue
  data_to_save <- splitted_data_frame[,trues]
  
  #Reading metadata provided by Oliva et al. which includes the PEER factors for the donors with genotype
  oliva_metadata <- read.delim(paste0(data_path, "v9/Oliva/", tissue_name, ".covariates.txt"))
  rownames(oliva_metadata) <- oliva_metadata[,1]
  oliva_metadata <- oliva_metadata[,2:ncol(oliva_metadata)]
  oliva_metadata <- as.data.frame(t(oliva_metadata))
  oliva_metadata$SUBJID <- rownames(oliva_metadata)
  oliva_metadata$SUBJID <- sapply(oliva_metadata$SUBJID, function(x) gsub("\\.", "-", x))
  
  #Adding metadata on inferred ancestry only available for donors with genotype
  ancestry <- read.delim(paste0(data_path, "v8/metadata/inferred_ancestry_838donors.txt"))
  names(ancestry) <- c("SUBJID", "Ancestry")
  metadata_to_save <- merge(oliva_metadata, ancestry, by="SUBJID")
  
  #Adding the rest of donor's metadata
  donor_metadata <- read.delim(paste0(data_path, "v8/metadata/GTEx_Subject_Phenotypes.GRU.txt.gz"))
  donor_metadata <- donor_metadata[,c("SUBJID", "SEX", "AGE", "BMI", "TRISCHD", "DTHHRDY")]
  metadata_to_save <- merge(metadata_to_save, donor_metadata, by="SUBJID")
  saveRDS(metadata_to_save, file = paste0(project_path, "Tissues/", tissue_name, "/metadata.rds"))

  #Data to save only for all donors in the given tissue:
  colnames(data_to_save) <- sapply(colnames(data_to_save), function(x) paste0(strsplit(x, "-")[[1]][1:2], collapse = "-")) #Now the column name is the subject id
  print(dim(data_to_save))
  data_to_save <- data_to_save[,colnames(data_to_save) %in% metadata_to_save$SUBJID]
  print(dim(data_to_save))
  saveRDS(data_to_save, file = paste0(project_path, "Tissues/", tissue_name, "/data.rds"))
}
