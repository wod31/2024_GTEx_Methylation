
# library(minfi)
# library(lumi)
library(dplyr)
library(tidyr)
library(tidyverse)
## Input data #####
#library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

first_dir <- "/gpfs/"
annotation <- read.csv(paste0(first_dir, "scratch/bsc83/bsc83535/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv"))

project_path <- paste0(first_dir, "projects/bsc83/Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Ancestry", "Age", "Sex")

results_Age <- lapply(tissues, function(tis) 
  readRDS(paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Methylation/Tissues/",tis,"/DVP_Age.rds")))
names(results_Age) <- tissues

results_Ancestry <- lapply(tissues, function(tis) 
  readRDS(paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Methylation/Tissues/",tis,"/DVP_Ancestry.rds")))
names(results_Ancestry) <- tissues

non_sexual_tissues <- c("Lung", "ColonTransverse", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "WholeBlood")
results_Sex <- lapply(non_sexual_tissues, function(tis) 
  readRDS(paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Methylation/Tissues/",tis,"/DVP_Sex.rds")))
names(results_Sex) <- non_sexual_tissues

files <- list.files(paste0(first_dir, '/projects/bsc83/Projects/GTEx_v8/Methylation/Data/EpiMap/'), pattern='.bed.gz',full.names=T)

names_chrom <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "Breast", "MuscleSkeletal", "KidneyCortex", "Testis", "PBMC")
chromhmm <- lapply(names_chrom, function(tis)
  read.delim(files[grep(tis, files)], sep='\t', header=F))
names(chromhmm) <- tissues

### overlap chromhmm and cpgs tested #### 
ann_bed <- annotation[!is.na(annotation$MAPINFO) & !is.na(annotation$CHR),] %>%
  dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>% distinct()
ann_bed$chrom <- paste0('chr',ann_bed$chrom)
head(ann_bed)
ann_bed$start <- ann_bed$start-1

library(valr)
chromhmm_cpgs <- lapply(tissues, function(tis) {
  chrom_df <- chromhmm[[tis]][,c(1:4)]
  colnames(chrom_df) <- c('chrom','start','end','region')
  bed_intersect(ann_bed, chrom_df, suffix = c("_ann", "_chromhmm"))})
names(chromhmm_cpgs) <- tissues

#### enrichment -- chi-squared first ####

#From 18 states to 14

# enrichment of DSE in particular types of AS events
my_fisher <- function(type, tissue, trait, direction){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  print(tissue)
  print(trait)
  res <- get(paste0("results_", trait))[[tissue]]
  chrom_tissue <- as.data.frame(chromhmm_cpgs[[tissue]])
  chrom_tissue$region_chromhmm_new <- chrom_tissue$region_chromhmm
  chrom_tissue <- chrom_tissue[chrom_tissue$.overlap ==1,]
  # 
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TssFlnkD", "TssFlnk", "TssFlnkU","TssA")] <- "TSS"
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("EnhA2", "EnhA1","EnhWk","EnhG1", "EnhG2")] <- "Enh"
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("ReprPCWk","ReprPC")] <- "ReprPC"
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TxWk","Tx")] <- "Tx"
  
  type_df <- chrom_tissue[chrom_tissue$region_chromhmm_new == type,]
  other_type <- chrom_tissue[chrom_tissue$region_chromhmm_new != type,]
  if(direction=="Hypo"){
    type_diff <- nrow(type_df[type_df$name_ann %in% rownames(res[res$DiffLevene<0,]),])
    other_type_diff <- nrow(other_type[other_type$name_ann %in% rownames(res[res$DiffLevene<0,]),])
  }else if(direction=="Hyper"){
    type_diff <- nrow(type_df[type_df$name_ann %in% rownames(res[res$DiffLevene>0,]),])
    other_type_diff <- nrow(other_type[other_type$name_ann %in% rownames(res[res$DiffLevene>0,]),])
  }
  type_notdiff <- nrow(type_df) - type_diff
  other_type_notdiff <- nrow(other_type) - other_type_diff
  
  ### test significance
  m <- matrix(c(type_diff, type_notdiff, other_type_diff, other_type_notdiff), 2,2, byrow = T)
  # print(m)  
  
  m[is.na(m)] <- 0
  # print(m)
  f <- fisher.test(m)
  # print(f)
  return(list("region"=type, "oddsRatio"=f$estimate, "p_value"=f$p.value, "CI_down"=f$conf.int[1], "CI_up"=f$conf.int[2], "sample_size" = type_diff, "tissue"=tissue, "direction"=direction))
}
# Two-tailed Fisher test
#families <- as.vector(unique(chromhmm_cpgs$Lung$region_chromhmm))
families <- c('Enh','EnhBiv','Het','Quies','ReprPC','TSS','TssBiv','Tx','ZNF/Rpts')
fisher_results <- lapply(non_sexual_tissues, function(tis) lapply(names, function(trait) lapply(families, function(region) lapply(c("Hyper", "Hypo"), function(direction) my_fisher(region, tis, trait, direction)))))
names(fisher_results) <- non_sexual_tissues

for (name in non_sexual_tissues) {names(fisher_results[[name]]) <- names}
for (name in non_sexual_tissues) {for (trait in names(fisher_results[[name]])) {names(fisher_results[[name]][[trait]]) <- families}}

fisher_results_sexual <- lapply(c("Ovary", "Prostate", "Testis"), function(tis) lapply(names, function(trait) lapply(families, function(region) lapply(c("Hyper", "Hypo"), function(direction) my_fisher(region, tis, trait, direction)))))
names(fisher_results_sexual) <- c("Ovary", "Prostate", "Testis")

for (name in c("Ovary", "Prostate", "Testis")) {names(fisher_results_sexual[[name]]) <- names}
for (name in c("Ovary", "Prostate", "Testis")) {for (trait in names(fisher_results_sexual[[name]])) {names(fisher_results_sexual[[name]][[trait]]) <- families}}


final <- c(fisher_results, fisher_results_sexual)
saveRDS(final, paste0(first_dir, 'projects/bsc83/Projects/GTEx_v8/Methylation/Tissues/DVP_enrichment_chromhmm.rds'))


