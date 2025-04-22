#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to prepare the necessary metadata per tissue for gene expression
# @software version: R=4.2.2

#Set path

# Reading tissue information: names, abbreviations and colors
tissue_info <- read.csv("/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/data/public/tissue_abreviation.txt")
tissues <- tissue_info$SMTSD #4 tissues that we do not use
tissues <- tissues[!tissues %in% c("Cells - EBV-transformed lymphocytes", "Cells - Cultured fibroblasts")] #Non-tissues are excluded
tissues <- tissues[tissues %in% tissue_info$SMTSD[tissue_info$tissue %in% c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")]]

#Reading protected metadata
donor_metadata <- read.delim("/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/data/protected/GTEx_Subject_Phenotypes.GRU.txt.gz")
donor_metadata <- donor_metadata[, colnames(donor_metadata) %in% c("SUBJID", "DTHHRDY", "AGE", "SEX", "BMI")] #Variables of interest

sample_metadata <- read.delim("/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/data/protected/GTEx_Sample_Attributes.GRU.txt.gz")
ancestry_metadata <- read.delim("/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/data/protected/inferred_ancestry_838donors.txt")

#Reading and filtering gene information
gene_annotation <- read.delim("/gpfs/projects/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/data/public/gencode.v26.GRCh38.genes.bed")
colnames(gene_annotation) <- c("chr","start","end","strand","feature","ensembl.id","gene.name", "biotype","source") #Renaming variables
#Keeping only PC genes and lincRNAs
gene_annotation <- gene_annotation[gene_annotation$biotype %in% c("protein_coding","lincRNA"),]
#Excluding PAR genes
PAR_genes <- sapply(grep(".Y", gene_annotation$ensembl.id, value = T), function(gene)
  unlist(strsplit(gene, split = "_"))[[1]])
gene_annotation <- gene_annotation[-unlist(lapply(PAR_genes, function(gene)
  grep(gene, gene_annotation$ensembl.id))),]
#write.csv(gene_annotation, "data/public/gene_annotation.csv", row.names = F)


#Reading expression data (these files are not shared by us. They can be downloaded directly from the GTEx portal: https://gtexportal.org/home/datasets)
# data_dir <- "~/Documents/mn4/scratch/bsc83/bsc83535/GTEx/v8/" #My path to these heavy files
print("About to read expression data")
data_dir <- "/gpfs/scratch/bsc83/MN4/bsc83/bsc83535/GTEx/v8/" #My path to these heavy files
tpm <- read.delim(paste0(data_dir, "expression_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"), skip=2)
colnames(tpm) <- gsub("\\.","-", colnames(tpm))
rownames(tpm) <- tpm[,1]
tpm <- tpm[,-c(1:2)]
print("tpm read")

counts <- read.delim(paste0(data_dir, "expression_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"), skip = 2)
colnames(counts) <- gsub("\\.","-", colnames(counts))
rownames(counts) <- counts[,1]
counts <- counts[,-c(1:2)]
print("counts read")

# Subset TPM and count matrix for the genes we are interested in
tpm <- tpm[rownames(tpm) %in% gene_annotation$ensembl.id,]
counts <- counts[rownames(counts) %in% gene_annotation$ensembl.id,]

#Reading the PEER factors. The files are publicly available in the GTEx portal under "GTEx_Analysis_v8_eQTL_covariates.tar.gz"
peer_path <- paste0(data_dir, "cis_QTLs/cis_eQTLs/GTEx_Analysis_v8_eQTL_covariates/") # This is the path where we have located the extracted information



#Function to call later on 
keep_diseases <- function(table){
  if(length(table)==2 & sum(table>20)==2){
    return(TRUE)
  } else{return(FALSE)}
}

create_metadata <- function(tissue){ #Function to get the metadata per tissue
  metadata_subset <- sample_metadata[sample_metadata$SMTSD==tissue,] #Sample metadata for the tissue of interest
  metadata_subset <- metadata_subset[,colnames(metadata_subset) %in% c("SAMPID", "SMTSISCH", "SMRIN", "SMEXNCRT")] #Variables of interest
  metadata_subset$SUBJID <- sapply(metadata_subset$SAMPID, function(i) paste(unlist(strsplit(i,split = "-"))[1:2],collapse="." ) ) #Getting donor ID based on the sample ID
  
  #Adding PEER information in the metadata
  tissue_id <- tissue
  tissue_id <- gsub("\\(|\\)|-", "", tissue_id) #Replace "(", ")" and "-"
  tissue_id <- gsub("[[:space:]]+", "_", tissue_id) #Replace any amount of blank for an "_"
  if(length(grep(tissue_id, list.files(peer_path)))==0){return(NA)}
  print(tissue)
  peer_metadata <- read.delim(list.files(peer_path, full.names = T)[grep(tissue_id, list.files(peer_path))])
  peer_metadata <- peer_metadata[6:7,2:ncol(peer_metadata)] #PEER1 and PEER2
  peer_metadata <- t(peer_metadata)
  colnames(peer_metadata) <- c("PEER1", "PEER2")
  metadata_subset <- merge(metadata_subset, peer_metadata, by.x="SUBJID", by.y="row.names")
  
  #Adding donor metadata:
  metadata_subset$SUBJID <- gsub("\\.", "-", metadata_subset$SUBJID)
  metadata_subset <- merge(metadata_subset, donor_metadata, by="SUBJID")
  
  #Adding genetically inferred ancestry:
  metadata_subset <- merge(metadata_subset, ancestry_metadata, by.x="SUBJID", by.y="ID")
  metadata_subset <- metadata_subset[metadata_subset$inferred_ancestry!="ASN",]
  
  #Renaming and reordering variables
  colnames(metadata_subset) <- c("Donor", "Sample", "RIN", "IschemicTime", "ExonicRate", "PEER1", "PEER2", "Sex", "Age", "BMI", "HardyScale", "Ancestry")
  metadata_subset <- metadata_subset[,c("Donor", "Sample", "HardyScale", "IschemicTime", "RIN", "ExonicRate", "PEER1", "PEER2", "Age", "Ancestry", "Sex", "BMI")]
  metadata_subset$HardyScale <- as.factor(metadata_subset$HardyScale)
  metadata_subset$Ancestry <- as.factor(metadata_subset$Ancestry)
  metadata_subset$Sex <- as.factor(metadata_subset$Sex)
  
  metadata_subset <- na.omit(metadata_subset) #Removing samples with missing data
  
  tissue_id <- tissue_info$tissue[tissue_info$SMTSD==tissue] #Getting tissue name without space or _
  
  #Subsetting the tpm and counts per tissue
  count_tissue <- counts[,colnames(counts) %in% metadata_subset$Sample]
  tpm_tissue <- tpm[,colnames(tpm) %in% metadata_subset$Sample]
  
  #Subset of metadata for which we have expression data
  metadata_subset <- metadata_subset[metadata_subset$Sample %in% colnames(count_tissue),]
  #Renaming sample IDs
  metadata_subset$Sample <- sapply(metadata_subset$Sample, function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) ) #This shorter id will be to match the clinical annotation
  colnames(count_tissue) <- sapply(colnames(count_tissue), function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) ) #This shorter id will be to match the clinical annotation
  colnames(tpm_tissue) <- sapply(colnames(tpm_tissue), function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) ) #This shorter id will be to match the clinical annotation
  
  #Removing genes that are lowly expressed:
  #TPM>=0.1 in at least 20% of the tissue samples. 
  exprs_genes.tpm <- rownames(tpm)[apply(tpm, 1, function(x) sum(x>=0.1) ) >= 0.2*ncol(tpm)  ]
  #Count >=6 in at least 20% of the tissue samples. 
  exprs_genes.counts <- rownames(counts)[ apply(counts, 1, function(x) sum(x>=6) ) >= 0.2*ncol(counts)  ]
  exprs_genes <- intersect(exprs_genes.tpm,
                           exprs_genes.counts)
  
  # Excluding chrY genes in tissues associated to females
  if(tissue %in% c("Vagina", "Uterus", "Ovary")){
    Y_genes <- gene_annotation[gene_annotation$chr=="chrY",]$ensembl.id 
    exprs_genes <- exprs_genes[!exprs_genes %in% Y_genes]
  }
  
  #saving expression count data
  count_tissue <- count_tissue[exprs_genes,]
  tpm_tissue <- tpm_tissue[exprs_genes,]
  saveRDS(count_tissue, paste0("Tissues/", tissue_id, "/counts.rds"))
  saveRDS(tpm_tissue, paste0("Tissues/", tissue_id, "/tpm.rds"))
  saveRDS(exprs_genes, paste0("Tissues/", tissue_id, "/expressed_genes.rds"))
  
   saveRDS(metadata_subset, paste0("Tissues/", tissue_id, "/metadata_expression.rds"))
  return(metadata_subset)
}

metadata <- lapply(tissues, create_metadata)
