#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to get results from downsampling, median nº of DMPs
# @software version: R=4.2.2


#Set path 
setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")


#Plot results
library(parallel)

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "KidneyCortex", "Testis", "WholeBlood", "MuscleSkeletal")

table <- matrix(nrow = length(tissues), ncol = 3, dimnames = list(tissues, c("Ancestry", "Age", "Sex")))

read_file <- function(file) {
  data <- readRDS(file)
  list(age = sum(data$AGE$adj.P.Val < 0.05),
      ancestry = sum(data$EURv1$adj.P.Val < 0.05),
      sex = sum(data$SEX2$adj.P.Val < 0.05)) #If sexual tissue this will output a 0
}

# for(number in c(30,35,40,50,100)){
for(number in c(35)){
    for(tissue in tissues){
    Sys.time()
    print(tissue)
    files <- list.files(paste0("Downsampling/", number, "/", tissue, "/"), pattern = "results", full.names = T)
    
    start <- Sys.time()
    results <- mclapply(files, read_file, mc.cores = parallel::detectCores())
    
    medians_age <- sapply(results, function(x) x$age)
    medians_ancestry <- sapply(results, function(x) x$ancestry)
    medians_sex <- sapply(results, function(x) x$sex)
    
    table[tissue,"Ancestry"] <- median(medians_ancestry)
    table[tissue,"Age"] <- median(medians_age)
    table[tissue,"Sex"] <- median(medians_sex)
    
    stop <- Sys.time()
    print(stop-start)
  }
  
  print(table)
  saveRDS(table, paste0("Downsampling/downsampling_", number, ".rds"))
  Sys.time()
}

