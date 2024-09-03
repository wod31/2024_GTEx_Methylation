#!/usr/bin/env Rscript

# Libraries ----
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(ggrepel)

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

# Differential expression analyses: results tables ----
# for(tissue in tissues){
#   if(!file.exists(paste0("~/GTEx_v8/Raquel/03_DEA/01.DEA/Tissues/",tissue,"/", tissue,".voom_limma.covariates_and_traits.results.rds"))){print(tissue)}
# }
get_DMPs <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "SEX2"){
    NA
  }else{
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.rds"))[[trait]]
    #rownames(model[[trait]][model[[trait]]$adj.P.Val<0.05,])
    model
  }
}
DMPs_Res <- lapply(c('EURv1','SEX2','AGE','BMI'), function(trait) lapply(tissues, function(tissue) get_DMPs(tissue, trait)))
names(DMPs_Res) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(DMPs_Res[[trait]]) <- tissues}


# Lists of DEGs ----
get_DMPs <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "Sex"){
    NA
  }else{
    rownames(DMPs_Res[[trait]][[tissue]][DMPs_Res[[trait]][[tissue]]$adj.P.Val<0.05,])
  }
}
DMPs <- lapply(tissues, function(tissue) lapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait) get_DMPs(tissue, trait)))
names(DMPs) <- tissues
for(trait in tissues){names(DMPs[[trait]]) <- c("Ancestry", "Sex", "Age", "BMI")}


# ---- Analysis ----- ####

# Pairwise combination of traits --
pw_traits <- c("Ancestry:Age",
               "Sex:Age")

# List of DEGs per tissue
deg_pw_traits <- lapply(pw_traits, function(pw)
  lapply(tissues, function(tissue) 
    intersect(DMPs[[tissue]][[unlist(strsplit(pw, split = ":"))[[1]]]], DMPs[[tissue]][[unlist(strsplit(pw, split = ":"))[[2]]]])
  )
)
names(deg_pw_traits) <- pw_traits
for(pw in pw_traits){
  names(deg_pw_traits[[pw]]) <- tissues
}


# Ancestry:sex ----
tissue <- "Lung" 
trait1 <- "Sex"
trait2 <- "Age"
pw <- "Sex:Age"

d1 <- DMPs_Res[[trait1]][[tissue]][ deg_pw_traits[[pw]][[tissue]], ]
colnames(d1) <- paste0(colnames(d1), "_", trait1)
d2 <- DMPs_Res[[trait2]][[tissue]][ deg_pw_traits[[pw]][[tissue]], ]
colnames(d2) <- paste0(colnames(d2), "_", trait2)
d <- cbind.data.frame(d1, d2)

d$class <- NA
d$class[d$logFC_Sex < 0 & d$logFC_Age < 0] <- "1"
d$class[d$logFC_Sex > 0 & d$logFC_Age > 0] <- "1"
d$class[is.na(d$class)] <- "0"
d$class <- factor(d$class)
table(d$class)
class_cols <- c(brewer.pal(11, "RdBu")[c(2)], "grey")
#class_cols <- brewer.pal(11, "PuOr")[c(3, 9)]
names(class_cols) <- c("1", "0")
d$gene_name <- rownames(d)

p1 <- ggplot(d) +
  geom_point(aes( x = logFC_Sex,
                  y = logFC_Age,
                  col = class)) +
  theme_bw() + 
  geom_hline(yintercept = 0, lty= 2) +
  geom_vline(xintercept = 0, lty= 2) +
  xlab(paste0(trait1, " (logFC)")) +
  ylab(paste0(trait2, " (logFC)")) +
  #xlim(-2,2) + 
  labs( title = paste0(trait1, " - ", trait2, " DEGs (",  tissue_info[tissue_info$tissue_ID==tissue, "tissue_abbrv"],")")) +
  scale_color_manual(values = class_cols) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))  +
  geom_text_repel(data=d[d$logFC_Sex > 0.5 &
                           d$logFC_Age > 0.6, ], 
                  aes(logFC_Sex, logFC_Age, label=gene_name), 
                  max.overlaps = 10,
                  cex = 3) +
  geom_text_repel(data=d[d$logFC_Sex < -0.6 &
                           d$logFC_Age < -1, ], 
                  aes(logFC_Sex, logFC_Age, label=gene_name), 
                  max.overlaps = 10,
                  cex = 3) 

pdf("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/Additive_sex_age_lung.pdf",
    width = 5, height = 5)
p1 
dev.off()

saveRDS(d, '~/marenostrum/Projects/GTEx_v8/Methylation/Data/Additive_Sex_Age_lung.rds')

# Ancestry:Age ----
tissue <- "Ovary"
trait1 <- "Ancestry"
trait2 <- "Age"
pw <- "Ancestry:Age"

d1 <- DMPs_Res[[trait1]][[tissue]][ deg_pw_traits[[pw]][[tissue]], ]
colnames(d1) <- paste0(colnames(d1), "_", trait1)
d2 <- DMPs_Res[[trait2]][[tissue]][ deg_pw_traits[[pw]][[tissue]], ]
colnames(d2) <- paste0(colnames(d2), "_", trait2)
d <- cbind.data.frame(d1, d2)

d$class <- NA
d$class[d$logFC_Ancestry>0 & d$logFC_Age > 0] <- "1"
d$class[d$logFC_Ancestry<0 & d$logFC_Age < 0] <- "1"
d$class[is.na(d$class)] <- "0"
d$class <- factor(d$class)
table(d$class)
class_cols <- c(brewer.pal(11, "RdBu")[c(2)], "grey")
#class_cols <- brewer.pal(11, "PuOr")[c(3, 9)]
names(class_cols) <- c("1", "0")
d$gene_name <- rownames(d)

p3 <- ggplot(d) +
  geom_point(aes( x = logFC_Ancestry,
                  y = logFC_Age,
                  col = class)) +
  theme_bw() + 
  geom_hline(yintercept = 0, lty= 2) +
  geom_vline(xintercept = 0, lty= 2) +
  xlab(paste0(trait1, " (logFC)")) +
  ylab(paste0(trait2, " (logFC)")) +
  xlim(-2,2) + 
  labs( title = paste0(trait1, " - ", trait2, " DEGs (",  tissue_info[tissue_info$tissue_ID==tissue, "tissue_abbrv"],")")) +
  scale_color_manual(values = class_cols) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))  +
  geom_text_repel(data=d[d$logFC_Ancestry < -1 &
                           d$logFC_Age > 0.025, ], 
                  aes(logFC_Ancestry, logFC_Age, label=gene_name), max.overlaps = 10,
                  cex = 3) +
  geom_text_repel(data=d[d$logFC_Ancestry > 0.6 &
                           d$logFC_Age < -0.007,], 
                  aes(logFC_Ancestry, logFC_Age, label=gene_name), max.overlaps = 10,
                  cex = 3) 

pdf('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/scatter_ancestry_age_ovary.pdf', width = 4, height = 4)
p3
dev.off()

saveRDS(d, '~/marenostrum/Projects/GTEx_v8/Methylation/Data/Additive_Ancestry_Age_ovary.rds')

# library(ggpubr)
# pdf("~/GTEx_v8/Raquel/Draft/Analysis/expression/additive_effects_are_widespread_and_tissue_specific/figures/Figure_3E.artery_tibial_and_adipose_scatterplots.pdf",
#     width = 5, height = 10)
# ggarrange(p3, p1,nrow = 2)  
# dev.off()
# 
# 
# pdf("~/GTEx_v8/Raquel/Draft/Analysis/expression/figures/figure_3_scatterplot_artery_tibial_age_sex.pdf",
#     width = 5, height = 5)
# p3
# dev.off()

