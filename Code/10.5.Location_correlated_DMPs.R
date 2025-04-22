#!/usr/bin/env Rscript
# @Author: Winona Oliveros
# @E-mail: winona.oliveros@bsc.es
# @Description: Code to perform enrichment correlations EPIC annotation
# @software version: R=4.2.2

### Location ####
# ---- Data ----- ####
# Demographic traits ----
traits_cols <- c("Ancestry" = "#E69F00",
                 "Age" = "#56B4E9",
                 "Sex" =  "#009E73",
                 "BMI" = "#CC79A7")
traits <- names(traits_cols)

# Tissues ----
### read methylation results ####
first_dir <- "~/marenostrum/"
project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Age", "Ancestry", "BMI", "Sex")

tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Jose/00_Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

sex_tissues <- c('Ovary','Prostate','Testis')

tissues <- tissue_info$tissue_ID
# tissues cols
tissues_cols <- tissue_info$colcodes 
names(tissues_cols) <- tissue_info$tissue_abbrv

# Metadata ----
metadata <- lapply(tissues, function(tissue) {
  metadata_ind <- readRDS(paste0(project_path, "Tissues/",tissue, "/metadata.rds"))
  print("Reading Admixture results")
  admixture_ancestry <- read.table('~/marenostrum_scratch/GTEx/v8/genotype_data/admixture_inferred_ancestry.txt')
  colnames(admixture_ancestry) <- c('SUBJID','AFRv1','EURv1','inferred_ancestry','AFRv2','EURv2')
  metadata_ind <- merge(metadata_ind, admixture_ancestry[,c("SUBJID","EURv1")], by='SUBJID')
  metadata_ind
})
names(metadata) <- tissues

for(tissue in sex_tissues[c(1)]){
  metadata[[tissue]]$SEX <- "2"
  metadata[[tissue]]$SEX <- as.factor(metadata[[tissue]]$SEX)
  metadata[[tissue]] <- metadata[[tissue]][, colnames(metadata$MuscleSkeletal)[c(1:8,10:21)]]
}
for(tissue in sex_tissues[c(2,3)]){
  metadata[[tissue]]$SEX <- "1"
  metadata[[tissue]]$SEX <- as.factor(metadata[[tissue]]$SEX)
  metadata[[tissue]] <- metadata[[tissue]][, colnames(metadata$MuscleSkeletal)[c(1:8,10:21)]]
}

get_corr <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "SEX2"){
    NA
  }else{
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/",trait,'_Correlations_probes_genes_DEG_DMP.rds'))
    #rownames(model[[trait]][model[[trait]]$adj.P.Val<0.05,])
    model[model$p.adj<0.05,]
  }
}
DMPs_cor <- lapply(c('EURv1','SEX2','AGE','BMI'), function(trait) lapply(tissues, function(tissue) get_corr(tissue, trait)))
names(DMPs_cor) <-c('EURv1','SEX2','AGE','BMI')
for(trait in c('EURv1','SEX2','AGE','BMI')){names(DMPs_cor[[trait]]) <- tissues}

results_DML <- lapply(tissues, function(tis) 
  readRDS(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/",tis,"/DML_results_5_PEERs_continous.rds")))
names(results_DML) <- tissues

fisher_function <- function(tissue, trait, variable, direction) {
  print(variable)
  
  type_df <- DMPs_cor[[trait]][[tissue]]$probe[DMPs_cor[[trait]][[tissue]]$class==variable]
  other_type <- DMPs_cor[[trait]][[tissue]]$probe[DMPs_cor[[trait]][[tissue]]$class!=variable]
  
  if(direction=="hypo"){
    type_diff <- length(type_df[type_df %in% rownames(results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$logFC<0,]) & type_df %in% DMPs_cor[[trait]][[tissue]]$probe[DMPs_cor[[trait]][[tissue]]$p.adj<0.05]])
    type_notdiff <- length(type_df) - type_diff
    other_type_diff <- length(other_type[other_type %in% rownames(results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$logFC<0,]) & other_type %in% DMPs_cor[[trait]][[tissue]]$probe[DMPs_cor[[trait]][[tissue]]$p.adj<0.05]])
    other_type_notdiff <- length(other_type) - other_type_diff
  } else if(direction=="hyper"){
    type_diff <- length(type_df[type_df %in% rownames(results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$logFC>0,]) & type_df %in% DMPs_cor[[trait]][[tissue]]$probe[DMPs_cor[[trait]][[tissue]]$p.adj<0.05]])
    type_notdiff <- length(type_df) - type_diff
    other_type_diff <- length(other_type[other_type %in% rownames(results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$logFC>0,]) & other_type %in% DMPs_cor[[trait]][[tissue]]$probe[DMPs_cor[[trait]][[tissue]]$p.adj<0.05]])
    other_type_notdiff <- length(other_type) - other_type_diff
  }
  
  ### test significance
  m <- matrix(c(type_diff, type_notdiff, other_type_diff, other_type_notdiff), 2,2, byrow = T)
  f <- fisher.test(m)
  return(list("f" = f, "m" = type_diff))
}

families <- c("promoter", "enhancer", "gene_body")
fisher_results <- lapply(tissues, function(tis) lapply(names(results_DML[[tis]]), function(trait) lapply(families, function(region) fisher_function(tis,trait,region, "hypo"))))
names(fisher_results) <- tissues

for (name in tissues) {
  names(fisher_results[[name]]) <- names(results_DML[[name]])
}

for (name in tissues) {
  for (trait in names(fisher_results[[name]])) {
    names(fisher_results[[name]][[trait]]) <- families
  }
}

saveRDS(fisher_results, paste0(project_path,'Data/enrichment_anno_hypo_correlated.rds'))



fisher_results <- lapply(tissues, function(tis) lapply(names(results_DML[[tis]]), function(trait) lapply(families, function(region) fisher_function(tis,trait,region, "hyper"))))
names(fisher_results) <- tissues

for (name in tissues) {
  names(fisher_results[[name]]) <- names(results_DML[[name]])
}

for (name in tissues) {
  for (trait in names(fisher_results[[name]])) {
    names(fisher_results[[name]][[trait]]) <- families
  }
}

saveRDS(fisher_results, paste0(project_path,'Data/enrichment_anno_hyper_correlated.rds'))


### plot ####
