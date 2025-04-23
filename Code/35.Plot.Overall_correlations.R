#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Plot numbers of correlated DMPs
# @software version: R=4.2.2

library(ComplexHeatmap)
library(RColorBrewer)
suppressPackageStartupMessages(library(circlize))

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

tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Methylation/Data/Tissue_info_whole.rds"))

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
  admixture_ancestry <- read.table('~/marenostrum_scratch/MN4/bsc83/bsc83535/GTEx/v8/genotype_data/admixture_inferred_ancestry.txt')
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

to_plot_2 <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/correlations_to_plot_2.new.rds')
to_plot_1 <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/correlations_to_plot_1.new.rds')


#% of DMPs correlated
### read Correlations
get_corr <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "SEX2"){
    NA
  }else{
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/",trait,'_Correlations_probes_genes_DEG_DMP.rds'))
    model
  }
}
DMPs_Res <- lapply(c('EURv1','SEX2','AGE','BMI'), function(trait) lapply(tissues, function(tissue) get_corr(tissue, trait)))
names(DMPs_Res) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(DMPs_Res[[trait]]) <- tissues}
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


to_plot_1 <- c(1:5)
for (tissue in tissues) {
  for (trait in c("Ancestry", "Sex", "Age", "BMI")) {
    if(tissue %in% sex_tissues & trait == "Sex"){
      next
    } else {
      if(option==2){
        next
      } else {
        # all
        pro_all <- DMPs_Res[[trait]][[tissue]][!is.na(DMPs_Res[[trait]][[tissue]]$gene),]
        #pro_deg <- pro_all[pro_all$gene %in% degs,] #DEG-probe pairs
        bg <- length(unique(pro_all$gene)) # DEGs that have at least one probe in promoters
        corr <- pro_all[pro_all$p.adj<0.05,]
        
        ### gene numbers 
        # to_plot_2 <- rbind(to_plot_2,c(length(unique(corr[corr$cor>0,"gene"]))/bg, "Positive", tissue, trait,length(unique(corr[corr$cor>0,"gene"]))),
        #                    c(length(unique(corr[corr$cor<0,"gene"]))/bg, "Negative", tissue, trait,length(unique(corr[corr$cor<0,"gene"]))))
        
        bg <- length(unique(pro_all$probe)) # DMPs that are associated to a gene
        ### probe numbers 
        to_plot_1 <- rbind(to_plot_1,c(length(unique(corr[corr$cor>0,"probe"]))/bg, "Positive", tissue, trait,length(unique(corr[corr$cor>0,"probe"]))),
                           c(length(unique(corr[corr$cor<0,"probe"]))/bg, "Negative", tissue, trait,length(unique(corr[corr$cor<0,"probe"]))))
        
        #Promoters
        # pro_all <- DMPs_Res[[trait]][[tissue]][!is.na(DMPs_Res[[trait]][[tissue]]$gene) & DMPs_Res[[trait]][[tissue]]$class=="promoter",]
        # #pro_deg <- pro_all[pro_all$gene %in% degs,] #DEG-probe pairs
        # bg <- length(unique(pro_all$gene)) # DEGs that have at least one probe in promoters
        # corr <- pro_all[pro_all$p.adj<0.05,]
        # 
        # ### gene numbers 
        # to_plot_2 <- rbind(to_plot_2,c(length(unique(corr[corr$cor>0,"gene"]))/bg, "Promoter", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"gene"]))),
        #                    c(length(unique(corr[corr$cor<0,"gene"]))/bg, "Promoter", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"gene"]))))
        # 
        # bg <- length(unique(pro_all$probe)) # DMPs that are associated to a gene
        # ### probe numbers 
        # to_plot_1 <- rbind(to_plot_1,c(length(unique(corr[corr$cor>0,"probe"]))/bg, "Promoter", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"probe"]))),
        #                    c(length(unique(corr[corr$cor<0,"probe"]))/bg, "Promoter", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"probe"]))))
        # 
        # #Enhancers
        # pro_all <- DMPs_Res[[trait]][[tissue]][!is.na(DMPs_Res[[trait]][[tissue]]$gene) & DMPs_Res[[trait]][[tissue]]$class=="enhancer",]
        # #pro_deg <- pro_all[pro_all$gene %in% degs,] #DEG-probe pairs
        # bg <- length(unique(pro_all$gene)) # DEGs that have at least one probe in promoters
        # corr <- pro_all[pro_all$p.adj<0.05,]
        # 
        # ### gene numbers 
        # to_plot_2 <- rbind(to_plot_2,c(length(unique(corr[corr$cor>0,"gene"]))/bg, "Enhancer", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"gene"]))),
        #                    c(length(unique(corr[corr$cor<0,"gene"]))/bg, "Enhancer", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"gene"]))))
        # 
        # bg <- length(unique(pro_all$probe)) # DMPs that are associated to a gene
        # ### probe numbers 
        # to_plot_1 <- rbind(to_plot_1,c(length(unique(corr[corr$cor>0,"probe"]))/bg, "Enhancer", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"probe"]))),
        #                    c(length(unique(corr[corr$cor<0,"probe"]))/bg, "Enhancer", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"probe"]))))
        # 
        # #Gene Body
        # pro_all <- DMPs_Res[[trait]][[tissue]][!is.na(DMPs_Res[[trait]][[tissue]]$gene) & DMPs_Res[[trait]][[tissue]]$class=="gene_body",]
        # #pro_deg <- pro_all[pro_all$gene %in% degs,] #DEG-probe pairs
        # bg <- length(unique(pro_all$gene)) # DEGs that have at least one probe in promoters
        # corr <- pro_all[pro_all$p.adj<0.05,]
        # 
        # ### gene numbers 
        # to_plot_2 <- rbind(to_plot_2,c(length(unique(corr[corr$cor>0,"gene"]))/bg, "Gene Body", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"gene"]))),
        #                    c(length(unique(corr[corr$cor<0,"gene"]))/bg, "Gene Body", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"gene"]))))
        # 
        # bg <- length(unique(pro_all$probe)) # DMPs that are associated to a gene
        # ### probe numbers 
        # to_plot_1 <- rbind(to_plot_1,c(length(unique(corr[corr$cor>0,"probe"]))/bg, "Gene Body", "Positive", tissue, trait,length(unique(corr[corr$cor>0,"probe"]))),
        #                    c(length(unique(corr[corr$cor<0,"probe"]))/bg, "Gene Body", "Negative", tissue, trait,length(unique(corr[corr$cor<0,"probe"]))))
        # 
      }
    }
  }
}



saveRDS(to_plot_1, '/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/correlations_to_plot_1.new.all.rds')

to_plot_1 <- as.data.frame(to_plot_1[c(2:nrow(to_plot_1)),])
to_plot_1$V1 <- 100*as.numeric(to_plot_1$V1)
names(to_plot_1) <- c("N", "Correlation", 'Tissue','Trait',"Number")
#to_plot_1$type <- factor(to_plot_1$type, levels = c("Promoter", "Enhancer", "Gene Body"))


### plot numbers
library(ggh4x)
traits_cols <- c('#C49122','#4B8C61','#70A0DF','#A76595')
names(traits_cols) <- c("Ancestry", "Sex", "Age", "BMI")
strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols[1:3]))
to_plot_1$Trait <- factor(to_plot_1$Trait, levels = c("Ancestry", "Sex", "Age", "BMI"))
colors_types <- c('#274c77','#6096ba','#a3cef1')
colors_types <- c('#320A28','#511730','#8E443D')
colors_types <- c('#223843','#C6D7C1','#86AC8F')
colors_types <- c('#a3cef1','#8E443D')
pdf('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/corr_all_n_direction.pdf', width = 6, height = 2.5)
ggplot(to_plot_1[!is.na(to_plot_1$Tissue) & to_plot_1$Tissue %in% c('Testis','Ovary','Prostate','Lung','ColonTransverse','BreastMammaryTissue') & to_plot_1$Trait %in% c("Ancestry", "Sex", "Age"),], aes(y = Tissue, x=N, fill=Correlation)) +
  geom_bar(position = position_stack(reverse = TRUE), alpha=0.8, stat = 'identity') + 
  geom_text(aes(label=Number, x = as.numeric(N)), position = position_stack(reverse = TRUE), stat='identity', size=3,hjust=1.6, angle=20) +
  theme_bw() + ylab('') + xlab('Proportion correlated DMPs') + scale_fill_manual(values = colors_types)+
  facet_wrap2(~ Trait, strip = strip, nrow = 1) + theme_classic() + scale_x_continuous(breaks=seq(0, 100, 25))
dev.off()

### plot direction 
# Do we have more negative than expected?
binom <- list()
for (tissue in tissues[!tissues %in% c('KidneyCortex','MuscleSkeletal','WholeBlood')]) {
  corr <- all_cor_df[all_cor_df$Trait=='Ancestry' & all_cor_df$Tissue==tissue,]
  neg <- corr[corr$cor<0,]
  binom[[tissue]] <- binom.test(nrow(neg), nrow(corr), p = 0.5,
                                alternative = c("greater"),
                                conf.level = 0.95)
}
names(binom) <- tissues[!tissues %in% c('KidneyCortex','MuscleSkeletal','WholeBlood')]
