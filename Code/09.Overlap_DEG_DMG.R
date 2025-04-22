#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Overlap between Differentially expressed genes and differentially methylated genes
# @software version: R=4.2.2

#########################################
#### Overlap DMG and DEG #####
#########################################
### read DEG ####
GTEx_v8 <- readRDS('~/marenostrum/Projects/ribosomal_proteins/Winona/2022_Ribosomal_analysis/Data/Data_set_1.rds')
head(GTEx_v8)

### read methylation results ####
first_dir <- "~/marenostrum/"
project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Age", "Ancestry", "BMI", "Sex")

tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Jose/00_Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

sex_tissues <- c('Ovary','Prostate','Testis')
get_DMPs <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "SEX2"){
    NA
  }else{
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.rds"))
    rownames(model[[trait]][model[[trait]]$adj.P.Val<0.05,])
  }
}
DMPs <- lapply(c('EURv1','SEX2','AGE','BMI'), function(trait) lapply(tissues, function(tissue) get_DMPs(tissue, trait)))
names(DMPs) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(DMPs[[trait]]) <- tissues}

#########################################
### Get genes tested in array and DM ####
#########################################
library(tidyverse)
library(tidyr)
library(dplyr)
annotation <- read.csv("~/marenostrum_scratch/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv") ### procesed file
annotation <- annotation %>% separate_rows(UCSC_RefGene_Name, sep = ';')
annotation <- annotation %>% distinct()
## I will only consider genes associated to Enhancer/promoter and/or with a gene annotated on the manifest file
promoter_associated_probes <- annotation[annotation$Regulatory_Feature_Group=='Promoter_Associated',c("IlmnID","UCSC_RefGene_Name","Regulatory_Feature_Group", "Phantom5_Enhancers", "UCSC_RefGene_Group")]
enhancer_associated_probes <- annotation[annotation$Phantom5_Enhancers!='',c("IlmnID","UCSC_RefGene_Name","Regulatory_Feature_Group", "Phantom5_Enhancers", "UCSC_RefGene_Group")]
gene_body_associated_probes <- annotation[grep('Body',annotation$UCSC_RefGene_Group),c("IlmnID","UCSC_RefGene_Name","Regulatory_Feature_Group", "Phantom5_Enhancers", "UCSC_RefGene_Group")]
exon_body_associated_probes <- annotation[grep('1stExon',annotation$UCSC_RefGene_Group),c("IlmnID","UCSC_RefGene_Name","Regulatory_Feature_Group", "Phantom5_Enhancers", "UCSC_RefGene_Group")]
exon_gene <- rbind(gene_body_associated_probes, exon_body_associated_probes)
  
exon_gene <- exon_gene[!exon_gene$IlmnID %in% enhancer_associated_probes$IlmnID,]
exon_gene <- exon_gene[!exon_gene$IlmnID %in% promoter_associated_probes$IlmnID,]
enhancer_associated_probes <- enhancer_associated_probes[!enhancer_associated_probes$IlmnID %in% promoter_associated_probes$IlmnID,]
exon_gene$Type <- 'Gene_Associated'
enhancer_associated_probes$Type <- 'Enhancer_Associated'
promoter_associated_probes$Type <- 'Promoter_Associated'
gene_prom_enh_cpgs <- rbind(exon_gene, enhancer_associated_probes, promoter_associated_probes)
write.table(gene_prom_enh_cpgs, '~/marenostrum/Projects/GTEx_v8/Methylation/Data/Methylation_Epic_gene_promoter_enhancer_processed.txt', sep = '\t', 
            col.names = T, row.names = F, quote = F)

### get tested CpGs ####
get_DMPs_tested <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "SEX2"){
    NA
  }else{
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.rds"))
    rownames(model[[trait]])
  }
}
tested_cpgs <- lapply(c('EURv1','SEX2','AGE','BMI'), function(trait) lapply(tissues, function(tissue) get_DMPs_tested(tissue, trait)))
names(tested_cpgs) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(tested_cpgs[[trait]]) <- tissues}

#### get tested genes and DMG ######
get_DMGs <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "SEX2"){
    NA
  }else{
    unique(gene_prom_enh_cpgs$UCSC_RefGene_Name[gene_prom_enh_cpgs$IlmnID %in% DMPs[[trait]][[tissue]]])
  }
}
DMGs <- lapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait) lapply(tissues, function(tissue) get_DMGs(tissue, trait)))
names(DMGs) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(DMGs[[trait]]) <- tissues}

get_MGs <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "SEX2"){
    NA
  }else{
    unique(gene_prom_enh_cpgs$UCSC_RefGene_Name[gene_prom_enh_cpgs$IlmnID %in% tested_cpgs[[trait]][[tissue]]])
  }
}
MGs <- lapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait) lapply(tissues, function(tissue) get_MGs(tissue, trait)))
names(MGs) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(MGs[[trait]]) <- tissues}

#########################################
#### subset DEG #####
#########################################
get_DEGs <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "Sex"){
    NA
  }else{
    GTEx_v8[[tissue]][[trait]][GTEx_v8[[tissue]][[trait]][['adj.P.Val']]<0.05,'gene_name']
  }
}
tissues <- tissues[tissues !='KidneyCortex']
DEGs <- lapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait) lapply(tissues, function(tissue) get_DEGs(tissue, trait)))
names(DEGs) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(DEGs[[trait]]) <- tissues}

get_tested_DEG <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "Sex"){
    NA
  }else{
    GTEx_v8[[tissue]][[trait]][,'gene_name']
  }
}
tissues <- tissues[tissues !='KidneyCortex']
tested_DEGs <- lapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait) lapply(tissues, function(tissue) get_tested_DEG(tissue, trait)))
names(tested_DEGs) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(tested_DEGs[[trait]]) <- tissues}

#########################################
#### calculate overlap #####
#########################################

n.deg.1 <- (sapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait)
  sapply(tissues, function(tissue) {
    genes_common <- MGs[[trait]][[tissue]][MGs[[trait]][[tissue]] %in% tested_DEGs[[trait]][[tissue]]]
    DM <- DMGs[[trait]][[tissue]][DMGs[[trait]][[tissue]] %in% DEGs[[trait]][[tissue]]]
    return(length(DM))
  }
  )
))

n.deg.2 <- (sapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait)
  sapply(tissues, function(tissue) {
    genes_common <- MGs[[trait]][[tissue]][MGs[[trait]][[tissue]] %in% tested_DEGs[[trait]][[tissue]]]
    DM <- DMGs[[trait]][[tissue]][DMGs[[trait]][[tissue]] %in% DEGs[[trait]][[tissue]]]
    return(length(DM)/length(DMGs[[trait]][[tissue]][DMGs[[trait]][[tissue]] %in% genes_common]))
  }
  )
))

total_shared <- sapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait) length(unique(unlist(DMGs[[trait]]))[(unique(unlist(DMGs[[trait]]))) %in% unique(unlist(DEGs[[trait]]))]))

#########################################
#### plot overlap #####
#########################################
library(ComplexHeatmap)
library(RColorBrewer)
suppressPackageStartupMessages(library(circlize))

n_samples <- c()
for(tissue in tissues){ 
  print(tissue)
  metadata <- readRDS(paste0(project_path, "Tissues/",tissue, "/metadata.rds"))
  n_samples <- c(n_samples, nrow(metadata))
}
names(n_samples) <- tissues

my_pretty_num_function <- function(n){
  if(n==""){
    return(n)
  } else{
    prettyNum(n, big.mark = ",")
  }
}

tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Jose/00_Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

create_heatmap <- function(data, tissue_info, size=12){ #It takes as input the whole data frame, the whole information on tissues, the subset of tissues and  diseases to be plotted, and the font size
  #data, tissue_info, tissues_plot = rownames(data), diseases_plot = colnames(data), size=12
  without_NA <- replace(data, is.na(data), "")
  tissue_info <- tissue_info[tissues,]
  tissues_cols <- tissue_info[, 3]
  names(tissues_cols) <- tissue_info[, 1]
  tissues_cols <- tissues_cols[tissues]
  
  traits_cols <- c('#C49122','#4B8C61','#70A0DF','#A76595')
  names(traits_cols) <- c("Ancestry","Sex","Age","BMI")
  
  row_ha_left <- HeatmapAnnotation("Samples" = anno_barplot(n_samples,  #positives if we want to show the number of positives instead of n
                                                            gp = gpar(fill = tissues_cols,
                                                                      col = tissues_cols),
                                                            border=F, width = unit(1.5, "cm")),
                                   gap = unit(0.3,"cm"),
                                   show_legend = F,
                                   show_annotation_name = T,
                                   annotation_name_rot = 90,
                                   annotation_name_gp = gpar(fontsize = size),
                                   which = "row")
  column_ha_top <- HeatmapAnnotation("total unique number of DEGs" = anno_barplot(total_shared,
                                                                                  border = F,
                                                                                  gp = gpar(fill = traits_cols,
                                                                                            col = traits_cols)),
                                     show_annotation_name = T,
                                     annotation_name_side = "left",
                                     annotation_name_rot = 90,
                                     annotation_name_gp = gpar(fontsize = 6.5),
                                     height = unit(3, "cm"))
  
  
  Heatmap(data,
          heatmap_legend_param = list(legend_height = unit(5, "cm"),
                                      grid_width = unit(1, "cm"),
                                      labels_gp=gpar(fontsize=size),
                                      title_gp=gpar(fontsize=size, fontface=2),
                                      heatmap_legend_side = "bottom"),
          col= colorRamp2( c(0,1,max(data[!is.na(data)])/4,max(data[!is.na(data)])/2,max(data[!is.na(data)])),
                           brewer.pal(8, "BuPu")[c(1,2,4,5,7)]),
          # col=colorRamp2( c(0,1,3000),
          #                 c("white","#F5F8FC","#1266b5") ),
          na_col = "white",
          cluster_rows = F,
          cluster_columns = F,
          # name = "DE signal",
          name = "#DEG",
          row_names_side = "left",
          column_names_side = "bottom",
          column_names_rot =  90,
          column_names_gp = gpar(fontsize = size),
          column_names_max_height= unit(9, "cm"),
          row_names_gp = gpar(fontsize = size),
          left_annotation = row_ha_left,
          top_annotation = column_ha_top,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(my_pretty_num_function(without_NA[i, j]), x, y, gp = gpar(fontsize = size))}
          
  )
}

create_heatmap(n.deg.1, tissue_info, size=12)

###prop
n.deg.3 <- (sapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait)
  sapply(tissues, function(tissue) {
    genes_common <- tested_DEGs[[trait]][[tissue]][tested_DEGs[[trait]][[tissue]] %in% MGs[[trait]][[tissue]]]
    DM <- DEGs[[trait]][[tissue]][DEGs[[trait]][[tissue]] %in% DMGs[[trait]][[tissue]]]
    return(length(DM)/length(DEGs[[trait]][[tissue]][DEGs[[trait]][[tissue]] %in% genes_common]))
  }
  )
))

data <- cbind.data.frame("Tissue" = rep(tissues, 4),
                         "Trait" = unlist(lapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait) rep(trait, length(tissues)))),
                         "proportion_of_DMGs" = unlist(as.data.frame(n.deg.2)),
                         "proportion_of_DEGs" = unlist(as.data.frame(n.deg.3))
)
data$Tissue <- factor(data$Tissue, levels = rev(tissues), order = T)
data$Trait <- factor(data$Trait, levels = c("Ancestry", "Sex", "Age", "BMI"), order = T)

traits_cols <- c('#C49122','#4B8C61','#70A0DF','#A76595')
names(traits_cols) <- c("Ancestry","Sex","Age","BMI")

pdf(paste0(plot_path, "Figure_1C.gene_expression_variation_explained.pdf"),
    width = 3, height = 9)
ggplot(data,
       aes(x = Tissue,
           y = proportion_of_DMGs,
           col = Trait,
           size = proportion_of_DEGs)) +
  geom_point(alpha = 0.5) +
  coord_flip() +
  theme_bw() +
  scale_color_manual(values = traits_cols) +
  ylab("Nº of genes DM and DE (%)") +
  xlab("") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12))
dev.off()

pdf(paste0(plot_path, "Figure_1C.gene_expression_variation_explained.legend.pdf"),
    width = 5, height = 9)
ggplot(data,
       aes(x = Tissue,
           y = proportion_of_DMGs,
           col = Trait,
           size = proportion_of_DEGs)) +
  geom_point(alpha = 0.5) +
  coord_flip() +
  theme_bw() +
  scale_color_manual(values = traits_cols) +
  ylab("Nº of genes DM and DE (%)") +
  xlab("") +
  theme(#axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        #legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12))
dev.off()

### compare DM and DE prop overlap 
tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Jose/00_Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

tissue_info <- tissue_info[tissues,]
tissues_cols <- tissue_info[, 3]
names(tissues_cols) <- tissue_info[, 1]
tissues_cols <- tissues_cols[tissues]

ggplot(data,
       aes(x = proportion_of_DEGs,
           y = proportion_of_DMGs,
           col = Tissue,
           group = Trait)) +
  geom_point(aes(shape=Trait),alpha = 0.6, size=4) +
  coord_flip() +
  theme_bw() +
  scale_color_manual(values = tissues_cols) +
  ylab("Nº of genes DM (%)") +
  xlab("Nº of genes DE (%)") +
  geom_abline(intercept = 0, slope = 1)+
  theme(#axis.text.y = element_blank(),
    #axis.ticks.y = element_blank(),
    #legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 12))

#########################################
#### get overlapped genes to make enrichments #####
#########################################
get_overlap_names <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "Sex"){
    NA
  }else{
    DMGs[[trait]][[tissue]][DMGs[[trait]][[tissue]] %in% DEGs[[trait]][[tissue]]]
  }
}
tissues <- tissues[tissues !='KidneyCortex']
DM_names <- lapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait) lapply(tissues, function(tissue) get_overlap_names(tissue, trait)))
names(DM_names) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(DM_names[[trait]]) <- tissues}

get_overlap_names <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "Sex"){
    NA
  }else{
    tested_DEGs[[trait]][[tissue]][tested_DEGs[[trait]][[tissue]] %in% MGs[[trait]][[tissue]]]
  }
}
tissues <- tissues[tissues !='KidneyCortex']
common_names <- lapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait) lapply(tissues, function(tissue) get_overlap_names(tissue, trait)))
names(common_names) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(common_names[[trait]]) <- tissues}

saveRDS(common_names, '~/marenostrum/Projects/GTEx_v8/Methylation/Data/common_genes_tested_meth_expr.rds')
saveRDS(DM_names, '~/marenostrum/Projects/GTEx_v8/Methylation/Data/common_genes_diff_meth_expr.rds')

