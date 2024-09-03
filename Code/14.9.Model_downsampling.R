#!/usr/bin/env Rscript

# Parsing
library(optparse)
library(limma)

parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-t", "--tissue"), type="character",
                     dest="tissue",
                     help="Tissue")

options=parse_args(parser)
tissue=options$tissue

Sys.time()
#Create functions:
limma_function <- function(fit, x){
  covariate <<- x #makeContrast does not read the function's environment, so I add covariate to the general environment in my session
  if(covariate %in% nonEstimable(fit$design)){
    print("Covariate non estimable")
    return(0)
  }
  contrast.matrix <- suppressWarnings(makeContrasts(covariate, levels=fit$design)) #Warnings due to change of name from (Intercept) to Intercept
  fitConstrasts <- suppressWarnings(contrasts.fit(fit, contrast.matrix)) #Warning due to Intercept name
  eb = eBayes(fitConstrasts)
  tt.smart.sv <- topTable(eb,adjust.method = "BH",number=Inf)
  return(tt.smart.sv)
}
# tissue <- "Lung"
# tissue <- "WholeBlood"


print(tissue)

# first_dir <- "~/Documents/mn4/"
first_dir <- "/gpfs/"
project_path <- paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Methylation/")
#project_path <- paste0(first_dir, "/Projects/GTEx_v8/Methylation/")


print("Reading data")
Sys.time()
data <- readRDS(paste0(project_path,'Tissues/', tissue, "/data.rds")) #From whole compressed data in 5.6G to compressed 1.4G/1.1Gb only in Lung (the highest number of samples)
Sys.time() #12 minutes to load 15 Gb
probes <- rownames(data)
data <- sapply(data, as.numeric)
rownames(data) <- probes


print("Iterating:")
# for(number in c(30,40,50,100)){ #At first I was only doing 42 to match Muscle, but muscle was not being downsampled
for(number in c(35)){ 
    print(paste("Number of samples", number))
  Sys.time()
  # for(i in 1:50){
  for(i in 1:100){
      print(paste("i =", i))
    #Reading metadata
    metadata <- paste0(project_path, "Downsampling/", number, "/", tissue, "/metadata_downsampling_", i, ".rds")
    if(!file.exists(metadata)){
      print("no file")
      break}
    
    metadata <- readRDS(metadata)
    metadata$SEX <- as.factor(metadata$SEX)
    if(length(levels(metadata$SEX))==1){
      # print("Sexual tissue")
      metadata <- metadata[,-which(names(metadata) == "SEX")]
      individual_variables <- c("EURv1", "AGE", "BMI")
    } else{
      individual_variables <- c("EURv1", "SEX", "AGE", "BMI")
    }
    metadata$DTHHRDY <- as.factor(metadata$DTHHRDY)
    rownames(metadata) <- metadata$SUBJID
    metadata$SUBJID <- NULL
    
    #Subset beta values to keep only the needed samples
    beta <- data
    beta <- beta[,rownames(metadata)]   ### Note: the order should be the same in beta and metadata

    metadata_2 <- metadata[,c("PEER1", "PEER2", "PEER3", "PEER4", "PEER5", "TRISCHD", "DTHHRDY", individual_variables)]
    
    M <- log2(beta/(1-beta)) #Computing m-values to model
    
    if(i==1){
      print(as.formula(paste0("~", paste0(colnames(metadata_2), collapse="+"))))
    }
    mod <- model.matrix( as.formula(paste0("~", paste0(colnames(metadata_2), collapse="+"))), data =  metadata_2)
    
    fit <- lmFit(M, mod)
    
    to_run <- c("EURv1", "SEX2", "AGE", "BMI")
    if(!"SEX" %in% colnames(metadata_2)){
      to_run <- to_run[!to_run=="SEX2"]
    }
    res <- lapply(to_run, function(x) limma_function(fit, x))
    names(res) <- to_run
    
    output_name <- paste0(project_path, "Downsampling/", number, "/", tissue, "/results_downsampling_", i, ".rds")
    saveRDS(res, output_name) 
  }
}


