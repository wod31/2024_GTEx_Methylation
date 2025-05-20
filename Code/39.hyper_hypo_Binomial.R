#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Plot DMPs and get binomial for their directionality
# @software version: R=4.2.2


rm(list=ls())
suppressPackageStartupMessages(library(ComplexHeatmap))
library(RColorBrewer)
suppressPackageStartupMessages(library(circlize))

first_dir <- "~/marenostrum/"

project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Sex")
data <- matrix(nrow = length(tissues), ncol=1, dimnames = list(tissues, names))

tissue_info <- readRDS(paste0(project_path, "Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

n_samples <- c()
for(tissue in tissues){ 
  print(tissue)
  # model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results.rds"))
  # model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_batch.rds"))
  # model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_batch.rds"))
  model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.rds")) #Jose used the wrong name to the file DMR instead of DML
  if(TRUE %in% grepl("SEX",names(model))){
    data[tissue, "Sex"] <- sum(model$SEX2$adj.P.Val<0.05)# & model$SEX2$logFC<0)
  } else{data[tissue, "Sex"] <- NA}
  
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

data_up <- matrix(nrow = length(tissues), ncol=1, dimnames = list(tissues, names))
for(tissue in tissues){ 
  print(tissue)
  # model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results.rds"))
  # model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_batch.rds"))
  # model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_batch.rds"))
  model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.rds")) #Jose used the wrong name to the file DMR instead of DML
  
  if(TRUE %in% grepl("SEX",names(model))){
    data_up[tissue, "Sex"] <- sum(model$SEX2$adj.P.Val<0.05 & model$SEX2$logFC>0)
  } else{data_up[tissue, "Sex"] <- NA} 
  
}

data_down <- matrix(nrow = length(tissues), ncol=1, dimnames = list(tissues, names))
for(tissue in tissues){ 
  print(tissue)
  # model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results.rds"))
  # model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_batch.rds"))
  # model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_batch.rds"))
  model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.rds")) #Jose used the wrong name to the file DMR instead of DML
  
  if(TRUE %in% grepl("SEX",names(model))){
    data_down[tissue, "Sex"] <- sum(model$SEX2$adj.P.Val<0.05 & model$SEX2$logFC<0)
  } else{data_down[tissue, "Sex"] <- NA} 
  
}

data <- data[,c("Ancestry", "Sex", "Age", "BMI")]
create_heatmap <- function(data, tissue_info, size=12){ #It takes as input the whole data frame, the whole information on tissues, the subset of tissues and  diseases to be plotted, and the font size
  #data, tissue_info, tissues_plot = rownames(data), diseases_plot = colnames(data), size=12
  without_NA <- replace(data, is.na(data), "")
  
  tissues_cols <- tissue_info[, 3]
  names(tissues_cols) <- tissue_info[, 1]
  tissues_cols <- tissues_cols[tissues]
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
  
  Heatmap(data,
          heatmap_legend_param = list(legend_height = unit(5, "cm"),
                                      grid_width = unit(1, "cm"),
                                      labels_gp=gpar(fontsize=size),
                                      title_gp=gpar(fontsize=size, fontface=2)),
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
          column_names_side = "top",
          column_names_rot =  60,
          column_names_gp = gpar(fontsize = size),
          column_names_max_height= unit(9, "cm"),
          row_names_gp = gpar(fontsize = size),
          left_annotation = row_ha_left,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(my_pretty_num_function(without_NA[i, j]), x, y, gp = gpar(fontsize = size))}
          
  )
}

create_heatmap(data, tissue_info, size=12)

#### number hyper hypo barplot ####
sex_tissues <- c('Ovary','Prostate','Testis')
names <- c("Ancestry", "Sex", "Age", "BMI")

get_tissue_proportion <- function(tissue){
  print(tissue)
  #dea <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_batch_5_PEERs.rds"))
  if(tissue %in% sex_tissues){
    tissue_expression_variation_explained <- list(Hypo=sapply(names[-2], function(trait) data_down[tissue,trait]/sum(data_up[tissue,trait],data_down[tissue,trait])),
                                                  Hyper=sapply(names[-2], function(trait) data_up[tissue,trait]/sum(data_up[tissue,trait],data_down[tissue,trait])))
    tissue_expression_variation_explained <- list(Hypo=c(tissue_expression_variation_explained$Hypo[c(1)], 0, tissue_expression_variation_explained$Hypo[c(2,3)]),
                                                  Hyper=c(tissue_expression_variation_explained$Hyper[c(1)], 0, tissue_expression_variation_explained$Hyper[c(2,3)]))
    names(tissue_expression_variation_explained$Hypo)[2] <- "Sex"
    names(tissue_expression_variation_explained$Hyper)[2] <- "Sex"
  }else{
    tissue_expression_variation_explained <- list(Hypo=sapply(names, function(trait) data_down[tissue,trait]/sum(data_up[tissue,trait],data_down[tissue,trait])),
                                                  Hyper=sapply(names, function(trait) data_up[tissue,trait]/sum(data_up[tissue,trait],data_down[tissue,trait])))
  }
  return(tissue_expression_variation_explained)   
}
tissue_expression_variation_explained <- lapply(tissues, function(tissue) get_tissue_proportion(tissue))
names(tissue_expression_variation_explained) <- tissues

tissue_expression_variation_explained_df <- as.data.frame(unlist(tissue_expression_variation_explained))
tissue_expression_variation_explained_df$label <- rownames(tissue_expression_variation_explained_df)
tissue_expression_variation_explained_df <- tissue_expression_variation_explained_df %>% separate(label, c("Tissue", "Direction","Trait"))
colnames(tissue_expression_variation_explained_df) <- c('Prop',"Tissue", "Direction","Trait")

## binomial
binom_test <- function(trait,tissue) {
  print(tissue)
  print(trait)
  #dea <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_batch_5_PEERs.rds"))
  if(tissue %in% sex_tissues){
    if (trait == 'Sex') {
      return(NA)
    } else {
      up <- data_up[tissue,trait]
      dea <- sum(data_up[tissue,trait],data_down[tissue,trait])
    }
    
  }else{
    up <- data_up[tissue,trait]
    dea <- sum(data_up[tissue,trait],data_down[tissue,trait])
  }
  if (dea < 1) {
    return(NA)
  }
  binom <- binom.test(up, dea, 0.5)
  return(binom)   
}

binom_res <- lapply(traits, function(trait) lapply(tissues, function(tissue) binom_test(trait,tissue)))
names(binom_res) <- traits

for (trait in traits) {
  names(binom_res[[trait]]) <- tissues
}

read_data <- function(variables, data, trait){ #Function to prepare data to plot and compute adjusted p value
  
  
  odds_ratio <- lapply(variables, function(tissue) data[[trait]][[tissue]]$estimate)
  adj.P.Val <- p.adjust(sapply(variables, function(tissue) data[[trait]][[tissue]]$p.value), method = "BH")
  CI_down <- lapply(variables, function(tissue) data[[trait]][[tissue]]$conf.int[1])
  CI_up <- lapply(variables, function(tissue) data[[trait]][[tissue]]$conf.int[2])
  #sample_size <- lapply(variables, function(tissue) data[[trait]][[tissue]][['m']])
  
  
  names(odds_ratio) <- variables
  names(adj.P.Val) <- variables
  names(CI_down) <- variables
  names(CI_up) <- variables
  #names(sample_size) <- variables
  
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
  
  # sample_size_df <- as.data.frame(unlist(sample_size))
  # sample_size_df$label <- variables
  # sample_size_df$type <- deparse(substitute(data))
  # colnames(sample_size_df) <- c('sample_size','tissue','type')
  
  all <- Reduce(function(x, y) merge(x, y, all=TRUE), list(odds_ratio_df, adj.P.Val_df, CI_down_df, CI_up_df))
  head(all)
  all$sig <- 'not Sig'
  all$sig[all$adjPvalue<0.05] <- 'Sig'
  all <- all[,c("tissue","oddsRatio","adjPvalue","CI_down","CI_up","sig","type")]
  return(all)
}

binom_all <- list()

for (trait in c('Sex','Age','Ancestry')) {
  if (trait == 'Sex') {
    tissues_to_test <- tissues[!tissues %in% sex_tissues]
  } 
  if (trait == 'Age') {
    tissues_to_test <- tissues[tissues != 'MuscleSkeletal']
  }
  if (trait == 'Ancestry') {
    tissues_to_test <- tissues
  }
  binom_all[[trait]] <- read_data(tissues_to_test, binom_res, trait)
}

saveRDS(binom_all, '~/marenostrum/Projects/GTEx_v8/Methylation/Data/binomial_DMPs.rds')

# bar plot
library(ggh4x)
traits_cols <- c('#C49122','#4B8C61','#70A0DF')
names(traits_cols) <- c('Ancestry','Sex','Age')
strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols))

tissue_expression_variation_explained_df <- tissue_expression_variation_explained_df[tissue_expression_variation_explained_df$Trait!='BMI',]
tissue_expression_variation_explained_df$Trait <- factor(tissue_expression_variation_explained_df$Trait, levels = c('Ancestry','Sex','Age'))
tissue_expression_variation_explained_df$Tissue <- factor(tissue_expression_variation_explained_df$Tissue, 
                                                          levels = rev(tissues))
#tissue_expression_variation_explained_df$Trait <- droplevels(tissue_expression_variation_explained_df$Trait)

tissue_expression_variation_explained_df$label <- paste0(tissue_expression_variation_explained_df$Direction,
                                                         ':',tissue_expression_variation_explained_df$Trait)

## new colors 
colors_traits <- list('Hyper:Age'=c('#3D7CD0'),'Hypo:Age'=c('#B4D6F6'),
                      'Hyper:Sex'=c('#3B734E'),'Hypo:Sex'=c('#89AA94'),
                      'Hyper:Ancestry'=c('#F0AE21'),'Hypo:Ancestry'=c('#F9DE8B'))

pdf('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/DMP.barplot.hyper.hypo.new.colors.pdf', height = 4, width = 6)
ggplot(tissue_expression_variation_explained_df, aes(fill=label, y=Tissue, x=Prop)) + 
  geom_bar(position="stack", stat="identity", alpha=0.8) + ylab('')+
  scale_fill_manual(values = colors_traits, labels = c("Old", "EA","Female",'Young','AA','Male')) + xlab('Proportion DMPs')+
  facet_wrap2(~ Trait, strip = strip) + theme_bw() + scale_x_continuous(breaks=seq(0, 1, 0.5))
dev.off()
