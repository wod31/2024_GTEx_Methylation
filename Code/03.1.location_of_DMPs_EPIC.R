#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Enrichment on genomic location for differentially methylated positions/loci
# @software version: R=4.2.2

library(dplyr)
library(tidyr)
library(tidyverse)

#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

Sys.time()
data_path <- "~/marenostrum_scratch/MN4/bsc83/bsc83535/GTEx/v9/Oliva/"

annotation <- read.csv(paste0(data_path, "GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv"))
Sys.time()

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Age", "Ancestry", "BMI", "Sex")
traits_to_use <- c('EURv1','SEX2','AGE','BMI')

results_DML <- lapply(tissues, function(tis) 
  readRDS(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/",tis,"/DML_results_5_PEERs_continous.rds")))
names(results_DML) <- tissues


fisher_function <- function(variable, direction, tissue, trait){
  print(variable)
  if (variable=="island"){
    type_df <- annotation[annotation$Relation_to_UCSC_CpG_Island=="Island","IlmnID"]
    other_type <- annotation[!annotation$Relation_to_UCSC_CpG_Island=="Island","IlmnID"]
  } else if(variable=="shelf"){
    type_df <- annotation[annotation$Relation_to_UCSC_CpG_Island %in% c("S_Shelf", "N_Shelf"),"IlmnID"]
    other_type <- annotation[!annotation$Relation_to_UCSC_CpG_Island %in% c("S_Shelf", "N_Shelf"),"IlmnID"]
  } else if(variable=="shore"){
    type_df <- annotation[annotation$Relation_to_UCSC_CpG_Island %in% c("S_Shore", "N_Shore"),"IlmnID"]
    other_type <- annotation[!annotation$Relation_to_UCSC_CpG_Island %in% c("S_Shore", "N_Shore"),"IlmnID"]
  } else if(variable=="open_sea"){
    type_df <- annotation[!annotation$Relation_to_UCSC_CpG_Island %in% c("Island", "S_Shelf", "N_Shelf", "S_Shore", "N_Shore"),"IlmnID"]
    other_type <- annotation[annotation$Relation_to_UCSC_CpG_Island %in% c("Island", "S_Shelf", "N_Shelf", "S_Shore", "N_Shore"),"IlmnID"]
  } 
  
  res <- results_DML[[tissue]][[trait]]
  
  if(direction=="hypo"){
    type_diff <- length(type_df[type_df %in% rownames(res[res$adj.P.Val<0.05 & res$logFC<0,])])
    type_notdiff <- length(type_df) - type_diff
    other_type_diff <- length(other_type[other_type %in% rownames(res[res$adj.P.Val<0.05 & res$logFC<0,])])
    other_type_notdiff <- length(other_type) - other_type_diff
  } else if(direction=="hyper"){
    type_diff <- length(type_df[type_df %in% rownames(res[res$adj.P.Val<0.05 & res$logFC>0,])])
    type_notdiff <- length(type_df) - type_diff
    other_type_diff <- length(other_type[other_type %in% rownames(res[res$adj.P.Val<0.05 & res$logFC>0,])])
    other_type_notdiff <- length(other_type) - other_type_diff
  }
  
  ### test significance
  m <- matrix(c(type_diff, type_notdiff, other_type_diff, other_type_notdiff), 2,2, byrow = T)
  f <- fisher.test(m)
  return(list("f" = f, "m" = type_diff))
}

families <- c("island", "shelf", "shore", "open_sea")
fisher_results <- lapply(tissues, function(tis) lapply(names(results_DML[[tis]]), function(trait) lapply(families, function(region) fisher_function(region, "hypo", tis, trait))))
names(fisher_results) <- tissues

for (name in tissues) {
  names(fisher_results[[name]]) <- names(results_DML[[name]])
}

for (name in tissues) {
  for (trait in names(fisher_results[[name]])) {
    names(fisher_results[[name]][[trait]]) <- families
  }
}

saveRDS(fisher_results, paste0('Tissues//final_enrichment_hypo_epic','.rds'))

fisher_results <- lapply(tissues, function(tis) lapply(names(results_DML[[tis]]), function(trait)lapply(families, function(region) fisher_function(region, "hyper", tis, trait))))
names(fisher_results) <- tissues

for (name in tissues) {
  names(fisher_results[[name]]) <- names(results_DML[[name]])
}

for (name in tissues) {
  for (trait in names(fisher_results[[name]])) {
    names(fisher_results[[name]][[trait]]) <- families
  }
}


saveRDS(fisher_results, paste0('Tissues/final_enrichment_hyper_epic','.rds'))

read_data <- function(variables, data, type, trait){ #Function to prepare data to plot and compute adjusted p value
  

  odds_ratio <- lapply(variables, function(tissue) data[[tissue]][[trait]][[type]][['f']]$estimate)
  adj.P.Val <- p.adjust(sapply(variables, function(tissue) data[[tissue]][[trait]][[type]][['f']]$p.value), method = "BH")
  CI_down <- lapply(variables, function(tissue) data[[tissue]][[trait]][[type]][['f']]$conf.int[1])
  CI_up <- lapply(variables, function(tissue) data[[tissue]][[trait]][[type]][['f']]$conf.int[2])
  sample_size <- lapply(variables, function(tissue) data[[tissue]][[trait]][[type]][['m']])
  
  
  names(odds_ratio) <- variables
  names(adj.P.Val) <- variables
  names(CI_down) <- variables
  names(CI_up) <- variables
  names(sample_size) <- variables
  
  odds_ratio_df <- as.data.frame(unlist(odds_ratio))
  odds_ratio_df$label <- variables
  odds_ratio_df$type <- deparse(substitute(data)) #Either hypo or hyper
  colnames(odds_ratio_df) <- c('oddsRatio', 'tissue','type')
  
  adj.P.Val_df <- as.data.frame(unlist(adj.P.Val))
  adj.P.Val_df$label <- variables
  adj.P.Val_df$type <- deparse(substitute(data))
  colnames(adj.P.Val_df) <- c('adjPvalue','tissue','type')
  
  CI_down_df <- as.data.frame(unlist(CI_down))
  CI_down_df$label <- variables
  CI_down_df$type <- deparse(substitute(data))
  colnames(CI_down_df) <- c('CI_down','tissue','type')
  
  CI_up_df <- as.data.frame(unlist(CI_up))
  CI_up_df$label <- variables
  CI_up_df$type <- deparse(substitute(data))
  colnames(CI_up_df) <- c('CI_up','tissue','type')
  
  sample_size_df <- as.data.frame(unlist(sample_size))
  sample_size_df$label <- variables
  sample_size_df$type <- deparse(substitute(data))
  colnames(sample_size_df) <- c('sample_size','tissue','type')
  
  all <- Reduce(function(x, y) merge(x, y, all=TRUE), list(odds_ratio_df, adj.P.Val_df, CI_down_df, CI_up_df, sample_size_df))
  head(all)
  all$sig <- 'not Sig'
  all$sig[all$adjPvalue<0.05] <- 'Sig'
  all <- all[,c("tissue","oddsRatio","adjPvalue","CI_down","CI_up","sig","type", "sample_size")]
  return(all)
}


hypo <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/final_enrichment_hypo_epic.rds')
hyper <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/final_enrichment_hyper_epic.rds')

colors_traits <- list('AGE'=c('#3D7CD0','#B4D6F6'),
                      'SEX2'=c('#3B734E','#89AA94'),
                      'EURv1'=c('#F0AE21','#F9DE8B'))

tissue_info <- readRDS(paste0('~/marenostrum/', "Projects/GTEx_v8/Methylation/Data/Tissue_info_whole.rds"))

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "KidneyCortex", "Testis", "WholeBlood","MuscleSkeletal")

library(ggpubr)
sex_tissues <- c('Ovary','Prostate','Testis')
for (type in c("island", "shelf", "shore", "open_sea")) {
  for (trait in c('AGE','SEX2','EURv1')) {
    if (trait == 'SEX2') {
      tissues_to_test <- tissues[!tissues %in% sex_tissues]
    } 
    if (trait == 'AGE') {
      tissues_to_test <- tissues[tissues != 'MuscleSkeletal']
    }
    if (trait == 'EURv1') {
      tissues_to_test <- tissues
    }
    hypo_d <- read_data(tissues_to_test, hypo, type, trait)
    hyper_d <- read_data(tissues_to_test, hyper, type, trait)
    hyper_hypo <- rbind(hypo_d, hyper_d)
    hyper_hypo$sig[hyper_hypo$sig =="Sig"] <- "FDR < 0.05"
    hyper_hypo$sig[hyper_hypo$sig =="not Sig"] <- "FDR >= 0.05"
    hyper_hypo$sig <- factor(hyper_hypo$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
    hyper_hypo$tissue <- factor(hyper_hypo$tissue, levels=tissues_to_test)
    hyper_hypo$type[hyper_hypo$type =="hypo"] <- "Hypomethylation"
    hyper_hypo$type[hyper_hypo$type =="hyper"] <- "Hypermethylation"
    g <- ggplot(hyper_hypo, aes(x=log2(oddsRatio), y=tissue, colour=type, alpha=sig)) +
      geom_errorbar(aes(xmin=log2(CI_down), xmax=log2(CI_up)), width=.3) +
      geom_vline(xintercept = 0) +
      #xlim(0,20) + #Only for Lung to show the 0
      geom_point(size=3) + ylab('') + theme_bw() +
      scale_colour_manual(values=colors_traits[[trait]]) +
      xlab("log2(Odds ratio)") +
      scale_alpha_discrete(range = c(0.4, 1), drop = FALSE) +
      theme(legend.title = element_blank(),
            axis.text.x = element_text(colour="black", size=13),
            axis.text.y = element_text(colour="black", size=14),
            legend.text = element_text(colour="black", size=13),
            axis.title.x = element_text(size=16),
            legend.spacing.y = unit(-0.05, "cm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", linewidth=1)) #+

    
    #Plot sample sizes:
    
    g2 <- ggplot(hyper_hypo) + geom_col(aes(sample_size, tissue, fill=type), width = 0.6) +
      theme_classic() + xlab("Number of DMPs") + ylab("") +
      scale_fill_manual(values=colors_traits[[trait]]) +
      theme(legend.position = "none",
            axis.text.x = element_text(colour="black", size=13),
            axis.text.y=element_blank(),  #remove y axis labels,
            axis.title.x = element_text(size=16)) +
      scale_x_continuous(n.breaks=3)

    p <- ggarrange(g, g2, labels = c("A", "B"),
                   common.legend = TRUE, legend = "right", widths = c(0.8,0.3))
    pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/enrichment", gsub('\\/','_',type),'_',trait,".v2.EPIC.pdf"), w = 8, h = 4)
    print(p)
    dev.off()
  }
}


