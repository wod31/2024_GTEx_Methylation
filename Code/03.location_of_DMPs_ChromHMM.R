#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Get fisher enrichments of DMPs per chromatin state and tissue
# @software version: R=4.2.2


###### Do enrichment of DMPs for each chromatin state per tissue and trait #####

library(dplyr)
library(tidyr)
library(tidyverse)


first_dir <- "/gpfs/projects/bsc83/"
annotation <- read.csv("/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")

project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Age", "Ancestry", "BMI", "Sex")
traits_to_use <- c('EURv1','SEX2','AGE','BMI')

tissues <- c('Ovary')

results_DML <- lapply(tissues, function(tis) 
  readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Tissues/",tis,"/DML_results_5_PEERs_continous.rds")))
names(results_DML) <- tissues

files <- list.files('/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/EpiMap/', pattern='.bed.gz',full.names=T)

names_chrom <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "Breast", "MuscleSkeletal", "KidneyCortex", "Testis", "PBMC")
names_chrom <- c('Ovary')
chromhmm <- lapply(names_chrom, function(tis)
  read.delim(files[grep(tis, files)], sep='\t', header=F))
names(chromhmm) <- tissues

### overlap chromhmm and cpgs tested #### 
# ann_granges <- annotation[!is.na(annotation$MAPINFO),] %>%
#   dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>% distinct() %>% 
#   makeGRangesFromDataFrame(keep.extra.columns=T)
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
### read sharing ####
sharing <- readRDS('/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/Sharing_DMP.rds')
shared_cpgs <- readRDS('/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/chromHMM_shared_9tissues.rds')

#From 18 states to 14

# enrichment of DSE in particular types of AS events
my_fisher <- function(type, tissue, trait){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  print(tissue)
  res <- results_DML[[tissue]][[trait]]
  chrom_tissue <- as.data.frame(chromhmm_cpgs[[tissue]])
  chrom_tissue$region_chromhmm_new <- chrom_tissue$region_chromhmm
  chrom_tissue <- chrom_tissue[chrom_tissue$.overlap ==1,]
  # 
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TssFlnkD", "TssFlnk", "TssFlnkU","TssA")] <- "TSS"
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("EnhA2", "EnhA1","EnhWk","EnhG1", "EnhG2")] <- "Enh"

  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("ReprPCWk","ReprPC")] <- "ReprPC"

  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TxWk","Tx")] <- "Tx"
  
  type_df <- chrom_tissue[chrom_tissue$region_chromhmm_new == type & chrom_tissue$name_ann %in% rownames(res),]
  other_type <- chrom_tissue[chrom_tissue$region_chromhmm_new != type & chrom_tissue$name_ann %in% rownames(res),]
  type_diff <- nrow(type_df[type_df$name_ann %in% rownames(res[res$adj.P.Val<0.05 & res$logFC<0,]),])
  type_notdiff <- nrow(type_df) - type_diff
  other_type_diff <- nrow(other_type[other_type$name_ann %in% rownames(res[res$adj.P.Val<0.05 & res$logFC<0,]),])
  other_type_notdiff <- nrow(other_type) - other_type_diff
  
  ### test significance
  m <- matrix(c(type_diff, type_notdiff, other_type_diff, other_type_notdiff), 2,2, byrow = T)
  print(m)  
  
  m[is.na(m)] <- 0
  #m <- m[c(type,paste0('No ',type)),]
  rownames(m) <- c(type, "Other")
  colnames(m) <- c("Hyper","Not Hyper")
  print(m)
  f <- fisher.test(m)
  print(f)
  return(list("f" = f, "m" = type_diff))
}
# Two-tailed Fisher test
#families <- as.vector(unique(chromhmm_cpgs$Lung$region_chromhmm))
# families <- c('Enh','EnhBiv','Het','Quies','ReprPC','TSS','TssBiv','Tx','ZNF/Rpts')
# fisher_results <- lapply(tissues, function(tis) lapply(names(results_DML[[tis]]), function(trait) lapply(families, function(region) my_fisher(region, tis,trait))))
# names(fisher_results) <- tissues
# 
# for (name in tissues) {
#   names(fisher_results[[name]]) <- names(results_DML[[name]])
# }
# 
# for (name in tissues) {
#   for (trait in names(fisher_results[[name]])) {
#     names(fisher_results[[name]][[trait]]) <- families
#   }
# }
# 
# saveRDS(fisher_results, '/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Tissues/Ovary.enrichment_chromhmm_hypo_batch_CI.simple.continous.filt.rds')

### enrichment shared positions
my_fisher <- function(type, trait){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  res <- sharing[sharing$trait == trait,]
  chrom_tissue <- shared_cpgs

  type_df <- chrom_tissue[chrom_tissue$region_chromhmm == type & chrom_tissue$name_ann %in% res$CG,]
  other_type <- chrom_tissue[chrom_tissue$region_chromhmm != type & chrom_tissue$name_ann %in% res$CG,]
  type_diff <- nrow(type_df[type_df$name_ann %in% res$CG[res$number>=2 & res$dir==1],])
  type_notdiff <- nrow(type_df) - type_diff
  other_type_diff <- nrow(other_type[other_type$name_ann %in% res$CG[res$number>=2 & res$dir==1],])
  other_type_notdiff <- nrow(other_type) - other_type_diff

  ### test significance
  m <- matrix(c(type_diff, type_notdiff, other_type_diff, other_type_notdiff), 2,2, byrow = T)
  print(m)

  m[is.na(m)] <- 0
  #m <- m[c(type,paste0('No ',type)),]
  rownames(m) <- c(type, "Other")
  colnames(m) <- c("Hyper","Not Hyper")
  print(m)
  f <- fisher.test(m)
  print(f)
  return(list("f" = f, "m" = type_diff))
}
# Two-tailed Fisher test
#families <- as.vector(unique(shared_cpgs$region_chromhmm))
families <- c('Enh','EnhBiv','Het','Quies','ReprPC','TSS','TssBiv','Tx','ZNF/Rpts')
fisher_results <- lapply(c("EURv1" ,"SEX2" , "AGE"  , "BMI"), function(trait) lapply(families, function(region) my_fisher(region,trait)))
names(fisher_results) <- c("Ancestry" ,"Sex" , "Age"  , "BMI")

for (name in c("Ancestry" ,"Sex" , "Age"  , "BMI")) {
  names(fisher_results[[name]]) <- families
}

saveRDS(fisher_results, '/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hyper_shared_CI.continous.2.rds')


#### enrichment overlap positions ancestry - Age colon transverse
# shared_positions <- read.delim('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Overlap_Ancestry_Age_COLON.txt')
# 
# first_dir <- "/gpfs/projects/bsc83/"
# annotation <- read.csv("/gpfs/scratch/bsc83/bsc83535/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")
# 
# project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")
# 
# tissues <- c("ColonTransverse")
# names <- c("Age", "Ancestry")
# traits_to_use <- c('EURv1','AGE')
# 
# results_DML <- lapply(tissues, function(tis) 
#   readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Tissues/",tis,"/DML_results_5_PEERs_continous.rds")))
# names(results_DML) <- tissues
# 
# files <- list.files('/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/EpiMap/', pattern='.bed.gz',full.names=T)
# 
# names_chrom <- c("ColonTransverse")
# chromhmm <- lapply(names_chrom, function(tis)
#   read.delim(files[grep(tis, files)], sep='\t', header=F))
# names(chromhmm) <- tissues
# 
# ### overlap chromhmm and cpgs tested #### 
# # ann_granges <- annotation[!is.na(annotation$MAPINFO),] %>%
# #   dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>% distinct() %>% 
# #   makeGRangesFromDataFrame(keep.extra.columns=T)
# ann_bed <- annotation[!is.na(annotation$MAPINFO) & !is.na(annotation$CHR),] %>%
#   dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>% distinct()
# ann_bed$chrom <- paste0('chr',ann_bed$chrom)
# head(ann_bed)
# ann_bed$start <- ann_bed$start-1
# 
# genes_age <- results_DML$ColonTransverse$AGE[results_DML$ColonTransverse$AGE$adj.P.Val<0.05,]
# genes_ancestry <- results_DML$ColonTransverse$EURv1[results_DML$ColonTransverse$EURv1$adj.P.Val<0.05,]
# 
# library(valr)
# chromhmm_cpgs <- lapply(tissues, function(tis) {
#   chrom_df <- chromhmm[[tis]][,c(1:4)]
#   colnames(chrom_df) <- c('chrom','start','end','region')
#   bed_intersect(ann_bed, chrom_df, suffix = c("_ann", "_chromhmm"))})
# names(chromhmm_cpgs) <- tissues
#   
# my_fisher <- function(type, tissue){
#   #             DS      Not DS
#   # type
#   # other_types
#   print(type)
#   res <- shared_positions[shared_positions$class==1,]
#   chrom_tissue <- as.data.frame(chromhmm_cpgs[[tissue]])
#   chrom_tissue$region_chromhmm_new <- chrom_tissue$region_chromhmm
#   chrom_tissue <- chrom_tissue[chrom_tissue$.overlap ==1,]
#   # 
#   chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TssFlnkD", "TssFlnk", "TssFlnkU","TssA")] <- "TSS"
#   chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("EnhA2", "EnhA1","EnhWk","EnhG1", "EnhG2")] <- "Enh"
#   
#   chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("ReprPCWk","ReprPC")] <- "ReprPC"
#   
#   chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TxWk","Tx")] <- "Tx"
#   
#   all_cpgs <- unique(c(unique(rownames(genes_ancestry)), unique(rownames(genes_age))))
#   
#   type_df <- chrom_tissue[chrom_tissue$region_chromhmm_new == type & chrom_tissue$name_ann %in% all_cpgs,]
#   other_type <- chrom_tissue[chrom_tissue$region_chromhmm_new != type & chrom_tissue$name_ann %in% all_cpgs,]
#   type_diff <- nrow(type_df[type_df$name_ann %in% rownames(res[res$logFC_Ancestry>0,]),])
#   type_notdiff <- nrow(type_df) - type_diff
#   other_type_diff <- nrow(other_type[other_type$name_ann %in% rownames(res[res$logFC_Ancestry>0,]),])
#   other_type_notdiff <- nrow(other_type) - other_type_diff
#   
#   ### test significance
#   m <- matrix(c(type_diff, type_notdiff, other_type_diff, other_type_notdiff), 2,2, byrow = T)
#   print(m)  
#   
#   m[is.na(m)] <- 0
#   #m <- m[c(type,paste0('No ',type)),]
#   rownames(m) <- c(type, "Other")
#   colnames(m) <- c("Hyper","Not Hyper")
#   print(m)
#   f <- fisher.test(m)
#   print(f)
#   return(list("f" = f, "m" = type_diff))
# }
# # Two-tailed Fisher test
# #families <- as.vector(unique(shared_cpgs$region_chromhmm))
# families <- c('Enh','EnhBiv','Het','Quies','ReprPC','TSS','TssBiv','Tx','ZNF/Rpts')
# fisher_results <- lapply(families, function(region) my_fisher(region,'ColonTransverse'))
# names(fisher_results) <- families
# 
# 
# saveRDS(fisher_results, '/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hyperEUR_hypoOLD_CI.continous.rds')

# Plot
# odds_ratio <- list()
# adj.P.Val <- list()
# # odds ratio
# odds_ratio <- lapply(tissues, function(tis) lapply(names(results_DML[[tis]]), function(trait) 
#   sapply(families, function(chr) fisher_results[[tis]][[trait]][[chr]][['f']]$estimate)))
# # adj.P.Val
# adj.P.Val <- lapply(tissues, function(tis) lapply(names(results_DML[[tis]]), function(trait) 
#   p.adjust(sapply(families, function(chr) fisher_results[[tis]][[trait]][[chr]][['f']]$p.value), method = "BH")))
# 
# names(odds_ratio) <- tissues
# names(adj.P.Val) <- tissues
# 
# for (name in tissues) {
#   names(odds_ratio[[name]]) <- names(results_DML[[name]])
# }
# for (name in tissues) {
#   names(adj.P.Val[[name]]) <- names(results_DML[[name]])
# }
# 
# odds_ratio_df <- as.data.frame(unlist(odds_ratio))
# odds_ratio_df$label <- rownames(odds_ratio_df)
# odds_ratio_df <- odds_ratio_df %>% separate(col=label, into=c('tissue', 'trait', 'region','type'), sep='\\.')
# colnames(odds_ratio_df) <- c('value','tissue', 'trait', 'region','type')
# 
# adj.P.Val_df <- as.data.frame(unlist(adj.P.Val))
# adj.P.Val_df$label <- rownames(adj.P.Val_df)
# adj.P.Val_df <- adj.P.Val_df %>% separate(col=label, into=c('tissue', 'trait', 'region','type'), sep='\\.')
# adj.P.Val_df$type <- 'adj.P.value'
# colnames(adj.P.Val_df) <- c('value','tissue', 'trait', 'region','type')
# 
# all_info_hypo <- rbind(odds_ratio_df,adj.P.Val_df)
# saveRDS(all_info_hypo, '/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm.rds')


###### plot ####
# all_info_hypo <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm.rds')
# # print(paste0(nrow(all_info_hypo[all_info_hypo$value<0.05 & all_info_hypo$type=='adj.P.value',])))
# 
# odds_ratio_df <- all_info_hypo[all_info_hypo$type == 'odds ratio',]
# adj.P.Val_df <- all_info_hypo[all_info_hypo$type == 'adj.P.value',]
# colnames(odds_ratio_df) <- c('OddsRatio','tissue', 'trait', 'region','type')
# colnames(adj.P.Val_df) <- c('adjPvalue','tissue', 'trait', 'region','type')
# #
# all <- merge(odds_ratio_df, adj.P.Val_df, by=c('tissue','trait','region'))
# head(all)
# all$sig <- 'not Sig'
# all$sig[all$adjPvalue<0.05] <- 'Sig'
# #
# # #### Bubble plot ####
# #
# # pdf('/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Tissues/')
# ggplot(all[all$trait=='SEX2',], aes(y=region, x=tissue, size=-log10(adjPvalue), fill=log2(OddsRatio))) +
#   geom_point(aes(color=sig), alpha=0.8, shape=21) +
#   scale_size(range = c(0, 21), name="-log10(adjPvalue)") +
#   scale_fill_gradient2(low="dark blue", high="red", name='OddsRatio', midpoint = 0)+
#   #theme_ipsum() +
#   theme(legend.position="right") +
#   ylab("") +
#   xlab("") +
#   theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
#         panel.grid = element_line(colour = 'light grey')) +
#   scale_color_manual(values = c("not Sig" = "white",
#                                 "Sig" = "black"))
# #dev.off()
# 

### Ovary and Colon new results 
# 
# first_dir <- "~/marenostrum/"
# annotation <- read.csv("~/marenostrum_scratch/MN4/bsc83/bsc83535/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")
# 
# project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")
# 
# tissues <- c("Ovary")
# names <- c("Age", "Ancestry")
# traits_to_use <- c('EURv1','AGE')
# 
# results_DML <- lapply(tissues, function(tis)
#   readRDS(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/",tis,"/DML_results_5_PEERs_continous.filt.rds")))
# names(results_DML) <- tissues
# 
# files <- list.files('/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/EpiMap/', pattern='.bed.gz',full.names=T)
# 
# names_chrom <- c("Ovary")
# chromhmm <- lapply(names_chrom, function(tis)
#   read.delim(files[grep(tis, files)], sep='\t', header=F))
# names(chromhmm) <- tissues
# 
# ### overlap chromhmm and cpgs tested ####
# # ann_granges <- annotation[!is.na(annotation$MAPINFO),] %>%
# #   dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>% distinct() %>%
# #   makeGRangesFromDataFrame(keep.extra.columns=T)
# ann_bed <- annotation[!is.na(annotation$MAPINFO) & !is.na(annotation$CHR),] %>%
#   dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>% distinct()
# ann_bed$chrom <- paste0('chr',ann_bed$chrom)
# head(ann_bed)
# ann_bed$start <- ann_bed$start-1
# 
# 
# library(valr)
# chromhmm_cpgs <- lapply(tissues, function(tis) {
#   chrom_df <- chromhmm[[tis]][,c(1:4)]
#   colnames(chrom_df) <- c('chrom','start','end','region')
#   bed_intersect(ann_bed, chrom_df, suffix = c("_ann", "_chromhmm"))})
# names(chromhmm_cpgs) <- tissues
# 
# my_fisher <- function(type, tissue){
#   #             DS      Not DS
#   # type
#   # other_types
#   print(type)
#   res <- results_DML[[tissue]][['AGE']]
#   chrom_tissue <- as.data.frame(chromhmm_cpgs[[tissue]])
#   chrom_tissue$region_chromhmm_new <- chrom_tissue$region_chromhmm
#   chrom_tissue <- chrom_tissue[chrom_tissue$.overlap ==1,]
#   #
#   chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TssFlnkD", "TssFlnk", "TssFlnkU","TssA")] <- "TSS"
#   chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("EnhA2", "EnhA1","EnhWk","EnhG1", "EnhG2")] <- "Enh"
# 
#   chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("ReprPCWk","ReprPC")] <- "ReprPC"
# 
#   chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TxWk","Tx")] <- "Tx"
# 
#   type_df <- chrom_tissue[chrom_tissue$region_chromhmm_new == type & chrom_tissue$name_ann %in% rownames(res),]
#   other_type <- chrom_tissue[chrom_tissue$region_chromhmm_new != type & chrom_tissue$name_ann %in% rownames(res),]
#   type_diff <- nrow(type_df[type_df$name_ann %in% rownames(res[res$adj.P.Val<0.05 & res$logFC>0,]),])
#   type_notdiff <- nrow(type_df) - type_diff
#   other_type_diff <- nrow(other_type[other_type$name_ann %in% rownames(res[res$adj.P.Val<0.05 & res$logFC>0,]),])
#   other_type_notdiff <- nrow(other_type) - other_type_diff
# 
#   ### test significance
#   m <- matrix(c(type_diff, type_notdiff, other_type_diff, other_type_notdiff), 2,2, byrow = T)
#   print(m)
# 
#   m[is.na(m)] <- 0
#   #m <- m[c(type,paste0('No ',type)),]
#   rownames(m) <- c(type, "Other")
#   colnames(m) <- c("Hyper","Not Hyper")
#   print(m)
#   f <- fisher.test(m)
#   print(f)
#   return(list("f" = f, "m" = type_diff))
# }
# # Two-tailed Fisher test
# #families <- as.vector(unique(shared_cpgs$region_chromhmm))
# families <- c('Enh','EnhBiv','Het','Quies','ReprPC','TSS','TssBiv','Tx','ZNF/Rpts')
# fisher_results <- lapply(families, function(region) my_fisher(region,'Ovary'))
# names(fisher_results) <- families
# 
# 
# saveRDS(fisher_results, '~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/Ovary.filt.enrichment_chromhmm_hyper.continous.rds')
# 
