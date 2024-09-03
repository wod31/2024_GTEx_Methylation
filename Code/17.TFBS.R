#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Enrichment on TFBS for differentially methylated positions/loci
# @software version: R=4.2.2

suppressMessages(library(dplyr))
library(parallel)
library(ggplot2)
suppressMessages(library(data.table))

#Set path 
setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")


#Variables
# tissue <- "Lng"
# tissue <- "Myo"
n_cores <- 16
# n_cores <- 1

# first_path <- "~/Documents/mn4/"
first_path <- "/gpfs/"
data_path <- paste0(first_path, "projects/bsc83/Projects/GTEx_v8/Methylation/TFBS/")

fisher_function <- function(target, dmps){
  target_positions <- unique(annotation$name_ann[annotation$TF==target])
  non_target_positions <- rownames(results)[!rownames(results) %in% target_positions]
  non_dmp <- rownames(results)[!rownames(results) %in% dmps & results$adj.P.Val<0.05]
  m <- matrix(c(sum(dmps %in% target_positions),
                sum(dmps %in% non_target_positions),
                sum(non_dmp %in% target_positions),
                sum(non_dmp %in% non_target_positions)), nrow=2)
  t <- fisher.test(m)
  return(list(t$p.value, t$estimate, t$conf.int[1], t$conf.int[2]))
}

#Splitting by chromatin state
# fisher_function <- function(target, region_cpgs, dmps){
#   bg <- rownames(results)[rownames(results) %in% region_cpgs]
#   target_positions <- unique(annotation$name_ann[annotation$TF==target & annotation$name_ann %in% bg]) #cpgs in the chromatin state where the TF binds
#   non_target_positions <- region_cpgs[!region_cpgs %in% target_positions] #cpgs in the chromatin state where the TF does not bind
#   non_dmp <- region_cpgs[!region_cpgs %in% dmps] #cpgs in the chromatin state that are not hyper/hypo. the background could be different, for instance dmps in the chromatin state, so we would be comparing hyper vs dmp
#   m <- matrix(c(sum(dmps %in% target_positions), #hyper cpgs in the given TFBS
#                 sum(dmps %in% non_target_positions),#hyper cpgs not in the given TFBS
#                 sum(non_dmp %in% target_positions), #non hyper cpgs in the given TFBS
#                 sum(non_dmp %in% non_target_positions)), nrow=2) #non hyper cpgs not in the given TFBS
#   t <- fisher.test(m)
#   return(list(t$p.value, t$estimate, t$conf.int[1], t$conf.int[2], target))
# }

Sys.time()


tissues <- list.dirs("Tissues/", full.names = F)[-1]
# hyper_enriched <- matrix(nrow=9, ncol=3, dimnames = list(tissues, c("AGE", "EURv1", "SEX2")))
# hypo_enriched <- matrix(nrow=9, ncol=3, dimnames = list(tissues, c("AGE", "EURv1", "SEX2")))

files <- list.files('/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/EpiMap/', pattern='.bed.gz',full.names=T)

names_chrom <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "Breast", "KidneyCortex", "Testis", "PBMC","MuscleSkeletal")
#names_chrom <- c("Breast", "ColonTransverse", "KidneyCortex", "Lung", "MuscleSkeletal", "Ovary", "Prostate", "Testis", "PBMC")

chromhmm <- lapply(names_chrom, function(tis)
  read.delim(files[grep(tis, files)], sep='\t', header=F))
names(chromhmm) <- tissues

annotation <- read.csv("/gpfs/scratch/bsc83/MN4/bsc83/bsc83535/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")
ann_bed <- annotation[!is.na(annotation$MAPINFO) & !is.na(annotation$CHR),] %>%
  dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>% distinct()
ann_bed$chrom <- paste0('chr',ann_bed$chrom)
head(ann_bed)
ann_bed$start <- ann_bed$start-1

# library(valr)
# chromhmm_cpgs <- lapply(tissues, function(tis) {
#   chrom_df <- chromhmm[[tis]][,c(1:4)]
#   colnames(chrom_df) <- c('chrom','start','end','region')
#   bed_intersect(ann_bed, chrom_df, suffix = c("_ann", "_chromhmm"))})
# names(chromhmm_cpgs) <- tissues
# 
# #names_table <- cbind(tissues, "short"=c("Brs", "Dig", "Kid", "Lng", "Myo", "ALL", "Prs", "ALL", "Bld"))
names_table <- cbind(tissues, "short"=c("Lng", "Dig", "ALL", "Prs", "Brs", "Kid", "ALL", "Bld", "Myo"))
# 
# families <- c('Enh','EnhBiv','Het','Quies','ReprPC','TSS','TssBiv','Tx','ZNF/Rpts')
# families <- unique(chromhmm_cpgs$Lung$region_chromhmm)

sharing <- readRDS('/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/Sharing_DMP.rds')
tissues_shared <- c('shared')
names_table <- cbind(tissues_shared, "short"=c('ALL'))

for(tissue in c('shared')){
#  for (chrom in c('ReprPC','TSS')) {
    
  print(tissue)
#  print(chrom)
  short <- names_table[names_table[,1]==tissue,2]
  
  
  print("Reading data")
  results_DML <- readRDS(paste0("Tissues/", 'Lung' , "/DML_results_5_PEERs_continous.rds"))
  
  annotation <- read.csv(paste0(data_path, "chip_seq_", short,".csv"))[,-1]
  colnames(annotation)[5] <- c("TF")
  print(paste("We will be testing", length(unique(annotation$TF)), "TFs"))
  
  to_save <- data.frame(TF="1", tissue="1",state="1", variable="1", methylation="1", direction="1", p_value="1", odds_ratio="1", CI_low="1", CI_high="1")
  to_save$p_value <- as.numeric(to_save$p_value) #These 4 lines are needed to run in previous versions of R
  to_save$odds_ratio <- as.numeric(to_save$odds_ratio)
  to_save$CI_low <- as.numeric(to_save$CI_low)
  to_save$CI_high <- as.numeric(to_save$CI_high)
  #for(variable in c("AGE", "EURv1", "SEX2")){
  for(variable in c("SEX2")){
    print(variable)
    results <- results_DML[[variable]]
    shared_cpgs <- sharing[sharing$number>=2 & sharing$trait==variable & sharing$dir %in% c(1,-1),]
    #region_cpgs <- chromhmm_cpgs[[tissue]]$name_ann[chromhmm_cpgs[[tissue]]$region_chromhmm==chrom & chromhmm_cpgs[[tissue]]$.overlap ==1] ## chromhmm
    #results <- results[rownames(results) %in% region_cpgs,]
    #results <- results[results$adj.P.Val<0.05,] ## add this to use as background the DMPs
    # results_hyper <- rownames(results[results$adj.P.Val<0.05 & results$logFC>0 & rownames(results) %in% shared_cpgs$CG,])
    # results_hypo <- rownames(results[results$adj.P.Val<0.05 & results$logFC<0 & rownames(results) %in% shared_cpgs$CG,])
    results_hyper <- shared_cpgs$CG[shared_cpgs$dir==1]
    results_hypo <- shared_cpgs$CG[shared_cpgs$dir==-1]
    
    if(length(results_hyper)>0){
      print("Computing fisher's for hypermethylated positions")
      Sys.time()
      hyper <- mclapply(unique(annotation$TF), function(target) fisher_function(target,results_hyper),
                        mc.cores=n_cores) 
      names(hyper) <- unique(annotation$TF)
      hyper <- do.call(rbind, hyper)
      hyper <- as.data.frame(hyper)
      hyper[,1] <- unlist(hyper[,1])
      hyper[,2] <- unlist(hyper[,2])
      hyper[,3] <- unlist(hyper[,3])
      hyper[,4] <- unlist(hyper[,4])
      hyper[,5] <- NULL
      colnames(hyper) <- c("p_value", "odds_ratio", "CI_low", "CI_high")
      hyper$TF <- rownames(hyper)
      #saveRDS(hyper, paste0("Tissues/", tissue , "/TFBS_hyper_", variable,".rds"))
      Sys.time()
      
      # hyper_enriched[tissue, variable] <- sum(hyper$adj_p_val<0.05 & hyper$odds_ratio>1)
      direction <- ifelse(hyper$odds_ratio>1, "Enriched", "Depleted")
      to_save <- rbind(to_save, cbind(tissue=tissue, state='shared', variable=variable, methylation="Hyper", direction=direction, hyper))
    }
  
    if(length(results_hypo)>0){
      print("Computing fisher's for hypomethylated positions")
      Sys.time()
      hypo <- mclapply(unique(annotation$TF), function(target) fisher_function(target, results_hypo),
                       mc.cores=n_cores)
      names(hypo) <- unique(annotation$TF)
      hypo <- do.call(rbind, hypo)
      hypo <- as.data.frame(hypo)
      hypo[,1] <- unlist(hypo[,1])
      hypo[,2] <- unlist(hypo[,2])
      hypo[,3] <- unlist(hypo[,3])
      hypo[,4] <- unlist(hypo[,4])
      hypo[,5] <- NULL
      colnames(hypo) <- c("p_value", "odds_ratio", "CI_low", "CI_high")
      hypo$TF <- rownames(hypo)
      
      # hypo$adj_p_val <- p.adjust(hypo$p_value, "BH")
      # saveRDS(hypo, paste0("Tissues/", tissue , "/TFBS_hypo_", variable,".rds"))
      Sys.time()
      
      # hypo_enriched[tissue, variable] <- sum(hypo$adj_p_val<0.05 & hypo$odds_ratio>1)
      direction <- ifelse(hypo$odds_ratio>1, "Enriched", "Depleted")
      to_save <- rbind(to_save, cbind(tissue=tissue, state='shared', variable=variable, methylation="Hypo", direction=direction, hypo))
    }
  }
  to_save <- to_save[-1,]
  to_save$adj_p_val <- p.adjust(to_save$p_value, "BH")
  to_save$direction <- factor(to_save$direction, levels=c("Enriched", "Depleted"))
  to_save$adj_p_val <- as.numeric(to_save$adj_p_val)
  to_save <- to_save[order(to_save$variable, to_save$methylation, to_save$direction, to_save$adj_p_val), ]
  write.csv(to_save, paste0(data_path, "tfbs_results_", tissue, "_tested_",'shared',"_dmpsbg.csv"))
  print("done")
  
}

