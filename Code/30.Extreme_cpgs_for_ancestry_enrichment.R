#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Get CpGs with extreme methylation values
# @software version: R=4.2.2


##### Check enrichment in xtreme regions ####
### first get extreme regions per tissue ####
first_dir <- "/gpfs/projects/bsc83/"

project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Age", "Ancestry", "BMI", "Sex")

tissue_info <- readRDS(paste0(project_path, "Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

get_extremes <- function(tissue,dir) {
  betas <- readRDS(paste0(project_path,'Tissues/', tissue, "/data.rds"))
  probes <- rownames(betas)
  betas <- sapply(betas, as.numeric)
  rownames(betas) <- probes
  
  if (dir=='hyper') {
    means <- rowMeans(betas)
    extr <- names(which(means >0.9))
    return(extr)
  } else {
    means <- rowMeans(betas)
    extr <- names(which(means <0.1))
    return(extr)
  }
  
}

extremes <- lapply(c('hyper','hypo'), function(dir) lapply(tissues, function(tissue) get_extremes(tissue,dir)))
names(extremes) <- c('hyper','hypo')

for (dir in c('hyper','hypo')) {
  names(extremes[[dir]]) <- tissues
}

saveRDS(extremes, '/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/extreme_cpgs_09_01.rds')
