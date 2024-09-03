#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to correlate gene expression and DNA methylation
# @software version: R=4.2.2

#Set path 
# setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")

# Tissues ----
### read methylation results ####
first_dir <- "/gpfs/projects/bsc83/"
project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Age", "Ancestry", "BMI", "Sex")

tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Jose/00_Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

sex_tissues <- c('Ovary','Prostate','Testis')

tissues <- tissue_info$tissue_ID

### read Correlations
get_corr <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "SEX2"){
    NA
  }else{
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/",trait,'_Correlations_probes_genes_DEG_DMP.rds'))
    #rownames(model[[trait]][model[[trait]]$adj.P.Val<0.05,])
    model
  }
}
DMPs_Res <- lapply(c('EURv1','SEX2','AGE','BMI'), function(trait) lapply(tissues, function(tissue) get_corr(tissue, trait)))
names(DMPs_Res) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(DMPs_Res[[trait]]) <- tissues}
# cor <- readRDS("tissues/ColonTransverse/Correlations.rds")
option <- 1
# option <- 3

#Option 1: Correlate all DMPs to all genes
#Option 2: Correlate all probes to all genes
#Option 3: Correlate all DMPs to all DEGs:

get_sig <- function(tissue, trait){
  print(paste0(trait,':',tissue))
  if(tissue %in% sex_tissues & trait == "Sex"){
    NA
  }else{
    DMPs_Res[[trait]][[tissue]] <- DMPs_Res[[trait]][[tissue]][!is.na(DMPs_Res[[trait]][[tissue]]$gene),]
    table(DMPs_Res[[trait]][[tissue]]$p.adj<0.05)
  }
}

lapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait) lapply(tissues, function(tissue) get_sig(tissue, trait)))

get_signif <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "Sex"){
    NA
  }else{
    DMPs_Res[[trait]][[tissue]][DMPs_Res[[trait]][[tissue]]$p.adj<0.05,]
  }
}
signif <- lapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait) lapply(tissues, function(tissue) get_signif(tissue, trait)))
names(signif) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(signif[[trait]]) <- tissues}

#Plot rho values:
plot <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "Sex"){
    NA
  }else{
    ggplot(signif[[trait]][[tissue]], aes(class, abs(cor))) + geom_violin() + 
      geom_boxplot(outlier.shape = NA,
                   notch = T,
                   width = 0.25) + theme_bw() + xlab("") + ylab("abs(rho)") +
  ggtitle(paste0(tissue,':',trait))
  }
}

library(ggplot2)
lapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait) lapply(tissues, function(tissue) plot(tissue, trait)))

### read DEG ####
# GTEx_v8 <- readRDS('/gpfs/projects/bsc83/Projects/ribosomal_proteins/Winona/2022_Ribosomal_analysis/Data/Data_set_1.rds')
# head(GTEx_v8)

to_plot_2 <- c(1:6)
to_plot_1 <- c(1:6)
for (tissue in tissues) {
  for (trait in c("Ancestry", "Sex", "Age", "BMI")) {
    if(tissue %in% sex_tissues & trait == "Sex"){
      next
    } else {
      if(option==2){
        #Percentage of DEGs significantly correlated with a DMP
        # degs <- GTEx_v8[[tissue]][[trait]][GTEx_v8[[tissue]][[trait]][['adj.P.Val']]<0.05,'gene_name']
        # 
        # model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.rds"))
        # if (trait == 'Sex') {
        #   trait2 <- 'SEX2'
        #   dmps <- rownames(model[[trait2]][model[[trait2]]$adj.P.Val<0.05,])
        # } else if (trait == 'Age') {
        #   trait2 <- 'AGE'
        #   dmps <- rownames(model[[trait2]][model[[trait2]]$adj.P.Val<0.05,])
        # } else if (trait == 'Ancestry') {
        #   trait2 <- 'EURv1'
        #   dmps <- rownames(model[[trait2]][model[[trait2]]$adj.P.Val<0.05,])
        # } else {
        #   dmps <- rownames(model[[trait]][model[[trait]]$adj.P.Val<0.05,])
        # }
        
        #Promoters
        pro_all <- DMPs_Res[[trait]][[tissue]][!is.na(DMPs_Res[[trait]][[tissue]]$gene) & DMPs_Res[[trait]][[tissue]]$class=="promoter",]
        #pro_deg <- pro_all[pro_all$gene %in% degs,] #DEG-probe pairs
        bg <- length(unique(pro_all$gene)) # DEGs that have at least one probe in promoters
        corr <- pro_all[pro_all$p.adj<0.05,]
        # to_plot_2 <- rbind(c(length(unique(corr[corr$cor>0,"gene"]))/bg, "Promoter", "Positive"),
        #                    c(length(unique(corr[corr$cor<0,"gene"]))/bg, "Promoter", "Negative"))
        #How many of these probes are dmp?
        to_plot_2 <- rbind(to_plot_2,c(length(unique(corr[corr$cor>0 & corr$probe %in% dmps,"gene"]))/bg, "Promoter", "Positive", "DMP", tissue, trait,length(unique(corr[corr$cor>0 & corr$probe %in% dmps,"gene"]))),
                           c(length(unique(corr[corr$cor>0 & !corr$probe %in% dmps,"gene"]))/bg, "Promoter", "Positive", "not", tissue, trait,length(unique(corr[corr$cor>0 & !corr$probe %in% dmps,"gene"]))),
                           c(length(unique(corr[corr$cor<0 & corr$probe %in% dmps,"gene"]))/bg, "Promoter", "Negative", "DMP", tissue, trait,length(unique(corr[corr$cor<0 & corr$probe %in% dmps,"gene"]))),
                           c(length(unique(corr[corr$cor<0 & !corr$probe %in% dmps,"gene"]))/bg, "Promoter", "Negative", "not", tissue, trait,length(unique(corr[corr$cor<0 & !corr$probe %in% dmps,"gene"]))))
        
        pro_dmp <- pro_all[pro_all$probe %in% dmps,] #DMP-gene pairs
        bg <- length(unique(pro_dmp$probe)) # DMPs that are associated to a gene
        corr <- pro_dmp[pro_dmp$p.adj<0.05,]
        # to_plot_1 <- rbind(c(length(unique(corr[corr$cor>0,"probe"]))/bg, "Promoter", "Positive"), c(length(unique(corr[corr$cor<0,"probe"]))/bg, "Promoter", "Negative"))
        to_plot_1 <- rbind(to_plot_1,c(length(unique(corr[corr$cor>0 & corr$gene %in% degs,"probe"]))/bg, "Promoter", "Positive", "DEG", tissue, trait,length(unique(corr[corr$cor>0 & corr$gene %in% degs,"probe"]))),
                           c(length(unique(corr[corr$cor>0 & !corr$gene %in% degs,"probe"]))/bg, "Promoter", "Positive", "not", tissue, trait,length(unique(corr[corr$cor>0 & !corr$gene %in% degs,"probe"]))),
                           c(length(unique(corr[corr$cor<0 & corr$gene %in% degs,"probe"]))/bg, "Promoter", "Negative", "DEG", tissue, trait,length(unique(corr[corr$cor<0 & corr$gene %in% degs,"probe"]))),
                           c(length(unique(corr[corr$cor<0 & !corr$gene %in% degs,"probe"]))/bg, "Promoter", "Negative", "not", tissue, trait,length(unique(corr[corr$cor<0 & !corr$gene %in% degs,"probe"]))))
        
        #Enhancers
        enh_all <- DMPs_Res[[trait]][[tissue]][DMPs_Res[[trait]][[tissue]]$class=="enhancer",]
        enh_deg <- enh_all[enh_all$gene %in% degs,]
        bg <- length(unique(enh_deg$gene)) # genes that have at least one probe in promoters
        corr <- enh_deg[enh_deg$p.adj<0.05,]
        # to_plot_2 <- rbind(to_plot_2, c(length(unique(corr[corr$cor>0,"gene"]))/bg, "Enhancer", "Positive"), c(length(unique(corr[corr$cor<0,"gene"]))/bg, "Enhancer", "Negative"))
        to_plot_2 <- rbind(to_plot_2,c(length(unique(corr[corr$cor>0 & corr$probe %in% dmps,"gene"]))/bg, "Enhancer", "Positive", "DMP", tissue, trait,length(unique(corr[corr$cor>0 & corr$probe %in% dmps,"gene"]))),
                           c(length(unique(corr[corr$cor>0 & !corr$probe %in% dmps,"gene"]))/bg, "Enhancer", "Positive", "not", tissue, trait,length(unique(corr[corr$cor>0 & !corr$probe %in% dmps,"gene"]))),
                           c(length(unique(corr[corr$cor<0 & corr$probe %in% dmps,"gene"]))/bg, "Enhancer", "Negative", "DMP", tissue, trait,length(unique(corr[corr$cor<0 & corr$probe %in% dmps,"gene"]))),
                           c(length(unique(corr[corr$cor<0 & !corr$probe %in% dmps,"gene"]))/bg, "Enhancer", "Negative", "not", tissue, trait,length(unique(corr[corr$cor<0 & !corr$probe %in% dmps,"gene"]))))
        
        enh_dmp <- enh_all[enh_all$probe %in% dmps,] #DEG-probe pairs
        bg <- length(unique(enh_dmp$probe)) # DMPs that are associated to a gene
        corr <- enh_dmp[enh_dmp$p.adj<0.05,]
        # to_plot_1 <- rbind(to_plot_1, c(length(unique(corr[corr$cor>0,"probe"]))/bg, "Enhancer", "Positive"), c(length(unique(corr[corr$cor<0,"probe"]))/bg, "Enhancer", "Negative"))
        to_plot_1 <- rbind(to_plot_1,c(length(unique(corr[corr$cor>0 & corr$gene %in% degs,"probe"]))/bg, "Enhancer", "Positive", "DEG", tissue, trait,length(unique(corr[corr$cor>0 & corr$gene %in% degs,"probe"]))),
                           c(length(unique(corr[corr$cor>0 & !corr$gene %in% degs,"probe"]))/bg, "Enhancer", "Positive", "not", tissue, trait,length(unique(corr[corr$cor>0 & !corr$gene %in% degs,"probe"]))),
                           c(length(unique(corr[corr$cor<0 & corr$gene %in% degs,"probe"]))/bg, "Enhancer", "Negative", "DEG", tissue, trait,length(unique(corr[corr$cor<0 & corr$gene %in% degs,"probe"]))),
                           c(length(unique(corr[corr$cor<0 & !corr$gene %in% degs,"probe"]))/bg, "Enhancer", "Negative", "not", tissue, trait,length(unique(corr[corr$cor<0 & !corr$gene %in% degs,"probe"]))))
        
        #Gene body
        body_all <- DMPs_Res[[trait]][[tissue]][DMPs_Res[[trait]][[tissue]]$class=="gene_body",]
        body_deg <- body_all[body_all$gene %in% degs,]
        bg <- length(unique(body_deg$gene)) # genes that have at least one probe in gene body
        corr <- body_deg[body_deg$p.adj<0.05,]
        # to_plot_2 <- rbind(to_plot_2, c(length(unique(corr[corr$cor>0,"gene"]))/bg, "Gene body", "Positive"), c(length(unique(corr[corr$cor<0,"gene"]))/bg, "Gene body", "Negative"))
        to_plot_2 <- rbind(to_plot_2,c(length(unique(corr[corr$cor>0 & corr$probe %in% dmps,"gene"]))/bg, "Gene body", "Positive", "DMP", tissue, trait,length(unique(corr[corr$cor>0 & corr$probe %in% dmps,"gene"]))),
                           c(length(unique(corr[corr$cor>0 & !corr$probe %in% dmps,"gene"]))/bg, "Gene body", "Positive", "not", tissue, trait,length(unique(corr[corr$cor>0 & !corr$probe %in% dmps,"gene"]))),
                           c(length(unique(corr[corr$cor<0 & corr$probe %in% dmps,"gene"]))/bg, "Gene body", "Negative", "DMP", tissue, trait,length(unique(corr[corr$cor<0 & corr$probe %in% dmps,"gene"]))),
                           c(length(unique(corr[corr$cor<0 & !corr$probe %in% dmps,"gene"]))/bg, "Gene body", "Negative", "not", tissue, trait,length(unique(corr[corr$cor<0 & !corr$probe %in% dmps,"gene"]))))
        
        body_dmp <- body_all[body_all$probe %in% dmps,] #DEG-probe pairs
        bg <- length(unique(body_dmp$probe)) # DMPs that are associated to a gene
        corr <- body_dmp[body_dmp$p.adj<0.05,]
        # to_plot_1 <- rbind(to_plot_1, c(length(unique(corr[corr$cor>0,"probe"]))/bg, "Gene body", "Positive"), c(length(unique(corr[corr$cor<0,"probe"]))/bg, "Gene body", "Negative"))
        to_plot_1 <- rbind(to_plot_1,c(length(unique(corr[corr$cor>0 & corr$gene %in% degs,"probe"]))/bg, "Gene body", "Positive", "DEG", tissue, trait,length(unique(corr[corr$cor>0 & corr$gene %in% degs,"probe"]))),
                           c(length(unique(corr[corr$cor>0 & !corr$gene %in% degs,"probe"]))/bg, "Gene body", "Positive", "not", tissue, trait,length(unique(corr[corr$cor>0 & !corr$gene %in% degs,"probe"]))),
                           c(length(unique(corr[corr$cor<0 & corr$gene %in% degs,"probe"]))/bg, "Gene body", "Negative", "DEG", tissue, trait,length(unique(corr[corr$cor<0 & corr$gene %in% degs,"probe"]))),
                           c(length(unique(corr[corr$cor<0 & !corr$gene %in% degs,"probe"]))/bg, "Gene body", "Negative", "not", tissue, trait,length(unique(corr[corr$cor<0 & !corr$gene %in% degs,"probe"]))))
        
       
      } else {
        #Promoters
        pro_all <- DMPs_Res[[trait]][[tissue]][!is.na(DMPs_Res[[trait]][[tissue]]$gene) & DMPs_Res[[trait]][[tissue]]$class=="promoter",]
        #pro_deg <- pro_all[pro_all$gene %in% degs,] #DEG-probe pairs
        bg <- length(unique(pro_all$gene)) # DEGs that have at least one probe in promoters
        corr <- pro_all[pro_all$p.adj<0.05,]
        
        ### gene numbers 
        to_plot_2 <- rbind(to_plot_2,c(length(unique(corr[corr$cor>0,"gene"]))/bg, "Promoter", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"gene"]))),
                           c(length(unique(corr[corr$cor<0,"gene"]))/bg, "Promoter", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"gene"]))))

        bg <- length(unique(pro_all$probe)) # DMPs that are associated to a gene
        ### probe numbers 
        to_plot_1 <- rbind(to_plot_1,c(length(unique(corr[corr$cor>0,"probe"]))/bg, "Promoter", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"probe"]))),
                           c(length(unique(corr[corr$cor<0,"probe"]))/bg, "Promoter", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"probe"]))))
        
        #Enhancers
        pro_all <- DMPs_Res[[trait]][[tissue]][!is.na(DMPs_Res[[trait]][[tissue]]$gene) & DMPs_Res[[trait]][[tissue]]$class=="enhancer",]
        #pro_deg <- pro_all[pro_all$gene %in% degs,] #DEG-probe pairs
        bg <- length(unique(pro_all$gene)) # DEGs that have at least one probe in promoters
        corr <- pro_all[pro_all$p.adj<0.05,]
        
        ### gene numbers 
        to_plot_2 <- rbind(to_plot_2,c(length(unique(corr[corr$cor>0,"gene"]))/bg, "Enhancer", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"gene"]))),
                           c(length(unique(corr[corr$cor<0,"gene"]))/bg, "Enhancer", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"gene"]))))
        
        bg <- length(unique(pro_all$probe)) # DMPs that are associated to a gene
        ### probe numbers 
        to_plot_1 <- rbind(to_plot_1,c(length(unique(corr[corr$cor>0,"probe"]))/bg, "Enhancer", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"probe"]))),
                           c(length(unique(corr[corr$cor<0,"probe"]))/bg, "Enhancer", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"probe"]))))
        
        #Gene Body
        pro_all <- DMPs_Res[[trait]][[tissue]][!is.na(DMPs_Res[[trait]][[tissue]]$gene) & DMPs_Res[[trait]][[tissue]]$class=="gene_body",]
        #pro_deg <- pro_all[pro_all$gene %in% degs,] #DEG-probe pairs
        bg <- length(unique(pro_all$gene)) # DEGs that have at least one probe in promoters
        corr <- pro_all[pro_all$p.adj<0.05,]
        
        ### gene numbers 
        to_plot_2 <- rbind(to_plot_2,c(length(unique(corr[corr$cor>0,"gene"]))/bg, "Gene Body", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"gene"]))),
                           c(length(unique(corr[corr$cor<0,"gene"]))/bg, "Gene Body", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"gene"]))))
        
        bg <- length(unique(pro_all$probe)) # DMPs that are associated to a gene
        ### probe numbers 
        to_plot_1 <- rbind(to_plot_1,c(length(unique(corr[corr$cor>0,"probe"]))/bg, "Gene Body", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"probe"]))),
                           c(length(unique(corr[corr$cor<0,"probe"]))/bg, "Gene Body", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"probe"]))))
        
      }
    }
  }
}



saveRDS(to_plot_2, '/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/correlations_to_plot_2.new.rds')
saveRDS(to_plot_1, '/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/correlations_to_plot_1.new.rds')


