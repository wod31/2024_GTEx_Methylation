##### analyze sex by age ####
rm(list=ls())
suppressPackageStartupMessages(library(ComplexHeatmap))
library(RColorBrewer)
suppressPackageStartupMessages(library(circlize))

first_dir <- "~/marenostrum/"

project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Age")
data <- matrix(nrow = length(tissues), ncol=1, dimnames = list(tissues, names))

tissue_info <- readRDS(paste0(project_path, "Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

n_samples <- c()
# sex-specific
for(tissue in tissues){ 
  print(tissue)
  if (tissue %in% c('Prostate','Testis')) {
    #model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.rds"))
    data[tissue, "Age"] <- NA
    
  } else if (tissue != 'Ovary') {
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.peer.AgeFemale.rds"))
    data[tissue, "Age"] <- sum(model$AGE$adj.P.Val<0.05)# & model$SEX2$logFC<0)
    
  } else if (tissue== 'Ovary') {
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.rds"))
    data[tissue, "Age"] <- sum(model$AGE$adj.P.Val<0.05)# & model$SEX2$logFC<0)
    
  }
  
  metadata <- readRDS(paste0(project_path, "Tissues/",tissue, "/metadata.rds"))
  n_samples <- c(n_samples, nrow(metadata[metadata$SEX==2,]))
}

names(n_samples) <- tissues

n_samples <- c()
# sex-specific
for(tissue in tissues){ 
  print(tissue)
  if (tissue %in% c('Ovary')) {
    #model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.rds"))
    data[tissue, "Age"] <- NA
    
  } else if (!tissue %in% c('Prostate','Testis')) {
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.peer.AgeMale.rds"))
    data[tissue, "Age"] <- sum(model$AGE$adj.P.Val<0.05)# & model$SEX2$logFC<0)
    
  } else if (tissue %in% c('Prostate','Testis')) {
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.rds"))
    data[tissue, "Age"] <- sum(model$AGE$adj.P.Val<0.05)# & model$SEX2$logFC<0)
    
  }
  
  metadata <- readRDS(paste0(project_path, "Tissues/",tissue, "/metadata.rds"))
  n_samples <- c(n_samples, nrow(metadata[metadata$SEX==2,]))
}

for (tissue in tissues) {
  metadata <- readRDS(paste0(project_path, "Tissues/",tissue, "/metadata.rds"))
  n_samples <- c(n_samples, nrow(metadata[metadata$SEX==1,]))
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
for(tissue in tissues[!tissues %in% c('Ovary')]){ 
  print(tissue)
  if (tissue %in% c('Prostate','Testis')) {
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.rds"))
  } else {
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.peer.AgeMale.rds"))
  }#Jose used the wrong name to the file DMR instead of DML 
  data_up[tissue, "Age"] <- sum(model$AGE$adj.P.Val<0.05 & model$AGE$logFC>0)
  
}

data_down <- matrix(nrow = length(tissues), ncol=1, dimnames = list(tissues, names))
for(tissue in tissues[!tissues %in% c('Ovary')]){ 
  print(tissue)
  if (tissue %in% c('Prostate','Testis')) {
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.rds"))
  } else {
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.peer.AgeMale.rds"))
  }#Jose used the wrong name to the file DMR instead of DML   
  data_down[tissue, "Age"] <- sum(model$AGE$adj.P.Val<0.05 & model$AGE$logFC<0)
 
}

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
          name = "#DMP",
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

pdf('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/DMP.heatmap.age.male.v2.pdf', height = 5, width = 4.5)
create_heatmap(data, tissue_info, size=12)
dev.off()

#### number hyper hypo barplot ####
names <- c("Age")

get_tissue_proportion <- function(tissue){
  print(tissue)

  tissue_expression_variation_explained <- list(Hypo=sapply(names, function(trait) data_down[tissue,trait]/sum(data_up[tissue,trait],data_down[tissue,trait])),
                                                  Hyper=sapply(names, function(trait) data_up[tissue,trait]/sum(data_up[tissue,trait],data_down[tissue,trait])))
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
 
  up <- data_up[tissue,trait]
  dea <- sum(data_up[tissue,trait],data_down[tissue,trait])
  if (is.na(dea)) {
    return(NA)
    } else if (dea < 1) {
    return(NA)
  } else {
  binom <- binom.test(up, dea, 0.5)
  return(binom)   
  }
}

traits <- 'Age'
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

tissues_to_test <- tissues[!tissues %in% c('MuscleSkeletal','BreastMammaryTissue','KidneyCortex','Ovary')]

binom_all[[trait]] <- read_data(tissues_to_test, binom_res, trait)

saveRDS(binom_all, '~/marenostrum/Projects/GTEx_v8/Methylation/Data/binomial_DMPs.rds')

# bar plot
library(ggh4x)
traits_cols <- c('#70A0DF')
names(traits_cols) <- c('Age')
strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols))

tissue_expression_variation_explained_df$Trait <- factor(tissue_expression_variation_explained_df$Trait, levels = c('Age'))
tissue_expression_variation_explained_df$Tissue <- factor(tissue_expression_variation_explained_df$Tissue, 
                                                          levels = rev(tissues))
#tissue_expression_variation_explained_df$Trait <- droplevels(tissue_expression_variation_explained_df$Trait)

tissue_expression_variation_explained_df$label <- paste0(tissue_expression_variation_explained_df$Direction,
                                                         ':',tissue_expression_variation_explained_df$Trait)

## new colors 
colors_traits <- list('Hyper:Age'=c('#3D7CD0'),'Hypo:Age'=c('#B4D6F6'))

pdf('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/DMP.barplot.hyper.hypo.age_male.Downsample.pdf', height = 4, width = 4)
ggplot(tissue_expression_variation_explained_df, aes(fill=label, y=Tissue, x=Prop)) + 
  geom_bar(position="stack", stat="identity", alpha=0.8) + ylab('')+
  scale_fill_manual(values = colors_traits, labels = c("Old",'Young')) + xlab('Proportion DMPs')+
  facet_wrap2(~ Trait, strip = strip) + theme_bw() + scale_x_continuous(breaks=seq(0, 1, 0.5))
dev.off()


## female
data <- matrix(nrow = length(tissues), ncol=1, dimnames = list(tissues, names))
for(tissue in tissues){ 
  print(tissue)
  if (tissue %in% c('Ovary')) {
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.rds"))
    data[tissue, "Age"] <- sum(model$AGE$adj.P.Val<0.05)# & model$SEX2$logFC<0)
    
  } else if (!tissue %in% c('Prostate','Testis')) {
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.peer.AgeFemale.Downsample.rds"))
    data[tissue, "Age"] <- sum(model$AGE$adj.P.Val<0.05)# & model$SEX2$logFC<0)
    
  }
}

data <- matrix(nrow = length(tissues), ncol=1, dimnames = list(tissues, names))
for(tissue in tissues){ 
  print(tissue)
  if (tissue %in% c('Prostate','Testis')) {
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.rds"))
    data[tissue, "Age"] <- sum(model$AGE$adj.P.Val<0.05)# & model$SEX2$logFC<0)
    
  } else if (!tissue %in% c('Ovary')) {
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.peer.AgeMale.Downsample.rds"))
    data[tissue, "Age"] <- sum(model$AGE$adj.P.Val<0.05)# & model$SEX2$logFC<0)
    
  }
}

my_pretty_num_function <- function(n){
  if(n==""){
    return(n)
  } else{
    prettyNum(n, big.mark = ",")
  }
}

data_up <- matrix(nrow = length(tissues), ncol=1, dimnames = list(tissues, names))
for(tissue in tissues[!tissues %in% c('Prostate','Testis')]){ 
  print(tissue)
  if (tissue %in% c('Ovary')) {
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.rds"))
  } else {
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.peer.AgeFemale.Downsample.rds"))
  }#Jose used the wrong name to the file DMR instead of DML 
  data_up[tissue, "Age"] <- sum(model$AGE$adj.P.Val<0.05 & model$AGE$logFC>0)
  
}

data_down <- matrix(nrow = length(tissues), ncol=1, dimnames = list(tissues, names))
for(tissue in tissues[!tissues %in% c('Prostate','Testis')]){ 
  print(tissue)
  if (tissue %in% c('Ovary')) {
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.rds"))
  } else {
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.peer.AgeFemale.Downsample.rds"))
  }#Jose used the wrong name to the file DMR instead of DML   
  data_down[tissue, "Age"] <- sum(model$AGE$adj.P.Val<0.05 & model$AGE$logFC<0)
  
}

create_heatmap <- function(data, tissue_info, size=12){ #It takes as input the whole data frame, the whole information on tissues, the subset of tissues and  diseases to be plotted, and the font size
  #data, tissue_info, tissues_plot = rownames(data), diseases_plot = colnames(data), size=12
  without_NA <- replace(data, is.na(data), "")
  
  tissues_cols <- tissue_info[, 3]
  names(tissues_cols) <- tissue_info[, 1]
  tissues_cols <- tissues_cols[tissues]
  # row_ha_left <- HeatmapAnnotation("Samples" = anno_barplot(n_samples,  #positives if we want to show the number of positives instead of n
  #                                                           gp = gpar(fill = tissues_cols,
  #                                                                     col = tissues_cols),
  #                                                           border=F, width = unit(1.5, "cm")),
  #                                  gap = unit(0.3,"cm"),
  #                                  show_legend = F,
  #                                  show_annotation_name = T,
  #                                  annotation_name_rot = 90,
  #                                  annotation_name_gp = gpar(fontsize = size),
  #                                  which = "row")
  
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
          name = "#DMP",
          row_names_side = "left",
          column_names_side = "top",
          column_names_rot =  60,
          column_names_gp = gpar(fontsize = size),
          column_names_max_height= unit(9, "cm"),
          row_names_gp = gpar(fontsize = size),
          #left_annotation = row_ha_left,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(my_pretty_num_function(without_NA[i, j]), x, y, gp = gpar(fontsize = size))}
          
  )
}

pdf('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/DMP.heatmap.age.Downsample.male.pdf', height = 5, width = 4.5)
create_heatmap(data, tissue_info, size=12)
dev.off()

#### number hyper hypo barplot ####
names <- c("Age")

get_tissue_proportion <- function(tissue){
  print(tissue)
  
  tissue_expression_variation_explained <- list(Hypo=sapply(names, function(trait) data_down[tissue,trait]/sum(data_up[tissue,trait],data_down[tissue,trait])),
                                                Hyper=sapply(names, function(trait) data_up[tissue,trait]/sum(data_up[tissue,trait],data_down[tissue,trait])))
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
  
  up <- data_up[tissue,trait]
  dea <- sum(data_up[tissue,trait],data_down[tissue,trait])
  if (is.na(dea)) {
    return(NA)
  } else if (dea < 1) {
    return(NA)
  } else {
    binom <- binom.test(up, dea, 0.5)
    return(binom)   
  }
}

traits <- 'Age'
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

tissues_to_test <- tissues[!tissues %in% c('MuscleSkeletal','BreastMammaryTissue','KidneyCortex','Ovary')]
tissues_to_test <- tissues[!tissues %in% c('MuscleSkeletal','Prostate','KidneyCortex','Testis','WholeBlood')]

binom_all[[trait]] <- read_data(tissues_to_test, binom_res, trait)

saveRDS(binom_all, '~/marenostrum/Projects/GTEx_v8/Methylation/Data/binomial_DMPs.rds')

# bar plot
library(ggh4x)
traits_cols <- c('#70A0DF')
names(traits_cols) <- c('Age')
strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols))

tissue_expression_variation_explained_df$Trait <- factor(tissue_expression_variation_explained_df$Trait, levels = c('Age'))
tissue_expression_variation_explained_df$Tissue <- factor(tissue_expression_variation_explained_df$Tissue, 
                                                          levels = rev(tissues))
#tissue_expression_variation_explained_df$Trait <- droplevels(tissue_expression_variation_explained_df$Trait)

tissue_expression_variation_explained_df$label <- paste0(tissue_expression_variation_explained_df$Direction,
                                                         ':',tissue_expression_variation_explained_df$Trait)

## new colors 
colors_traits <- list('Hyper:Age'=c('#3D7CD0'),'Hypo:Age'=c('#B4D6F6'))

pdf('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/DMP.barplot.hyper.hypo.age_female.Downsample.pdf', height = 4, width = 4)
ggplot(tissue_expression_variation_explained_df, aes(fill=label, y=Tissue, x=Prop)) + 
  geom_bar(position="stack", stat="identity", alpha=0.8) + ylab('')+
  scale_fill_manual(values = colors_traits, labels = c("Old",'Young')) + xlab('Proportion DMPs')+
  facet_wrap2(~ Trait, strip = strip) + theme_bw() + scale_x_continuous(breaks=seq(0, 1, 0.5))
dev.off()


#####Tissue Sharing######
final_table <- data.frame(cg="chr1",logFC=1, adj.P.Val=1, trait=1, tissue=1)
for(trait in c("AGE", "EURv1", "BMI", "SEX2")){
  print(trait)
  probes <- data.frame(cg="chr1",logFC=1, adj.P.Val=1,trait=1, tissue=1)
  for(tissue in tissues){ 
    print(tissue)
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.rds"))
    if (!trait %in% names(model)) {
      next}
    if(is.na((model[[trait]][[1]][1]))){
      next
    }
    res <- as.data.frame(model[[trait]][model[[trait]]$adj.P.Val<0.05,])
    if (nrow(res)<1) {
      next}
    res$tissue <- tissue
    res$trait <- trait
    res$sign <- sign(res$logFC)
    res$cg <- rownames(res)
    probes <- rbind(probes, res[,c('cg','logFC','adj.P.Val','trait','tissue')])
  }
  probes <- probes[-1,]
  final_table <- rbind(final_table, probes)
}
final_table <- final_table[-1,]

library(plyr) #Counting
#final_table$name <- paste0(final_table$seqnames,':',final_table$start,';',final_table$end,';',final_table$overlapping.genes)
counts <- ddply(final_table, .(cg, trait), nrow)
names(counts) <- c("CG", "trait", "number")
final_table$sign <- as.numeric(sign(final_table$logFC))
logfc <- ddply(final_table, .(cg, trait), summarise, sign=length(unique(sign))*sign(unique(sign)[1]))
logfc[abs(logfc$sign)>1,'sign'] <- abs(logfc$sign)[abs(logfc$sign)>1]
names(logfc) <- c("CG", "trait", "dir")

to_plot <- merge(counts, logfc, by=c("CG", "trait"))
library(ggplot2)
to_plot$dir <- as.factor(to_plot$dir)
ggplot(to_plot) + geom_jitter(aes(trait, number, col=dir), alpha=0.5) +
  theme_bw() + geom_violin(aes(x = trait, y = number, fill=trait),alpha=0.8) +
  scale_fill_manual(values = c('#70A0DF','#C49122','#4B8C61'))

table(to_plot$trait, to_plot$dir)
table(to_plot$trait, to_plot$number)
table(to_plot$trait, to_plot$dir, to_plot$number)

View(to_plot[to_plot$number>=4,])

## binomial shared 
# sex 
binom <- binom.test(nrow(to_plot[to_plot$number>=2 & to_plot$trait=='SEX2' & to_plot$dir==1,]), nrow(to_plot[to_plot$number>=2 & to_plot$trait=='SEX2' & to_plot$dir %in% c(1,-1),]), 0.5)
## age
binom <- binom.test(nrow(to_plot[to_plot$number>=3 & to_plot$trait=='AGE' & to_plot$dir==1,]), nrow(to_plot[to_plot$number>=3 & to_plot$trait=='AGE' & to_plot$dir %in% c(1,-1),]), 0.5)


##### Final plot ####
head(to_plot)
to_plot$number <- as.factor(to_plot$number)

traits_cols <- c('#C49122','#4B8C61','#70A0DF','#A76595')
names(traits_cols) <- c('EURv1','SEX2','AGE','BMI')
strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols))
to_plot$trait <- factor(to_plot$trait, levels = c('EURv1','SEX2','AGE','BMI'))

ggplot(to_plot, aes(y = number, fill=dir)) +
  geom_bar(position = 'fill', alpha=0.8) + 
  geom_text(aes(label=after_stat(count), x = after_stat(count+1.5)), stat='count', position='fill', size=3,hjust=.8, angle=45) +
  theme_bw() + ylab('Nº of Tissues') + xlab('Proportion shared CpG') +
  facet_wrap2(~ trait, strip = strip, nrow = 1) + theme_bw() + scale_x_continuous(breaks=seq(0, 1, 0.5))

### proportion of shared per variable ####
to_plot <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Sharing_DMP.rds')
total_DEGs <- sapply(c('EURv1','SEX2','AGE','BMI'), function(trait) nrow(to_plot[to_plot$trait==trait,]))

barplot(total_DEGs, col = traits_cols, names.arg=c("Ancestry","Sex","Age","BMI"), ylab = 'total unique nº of DMP')
plot.new()
legend("center",
       c("Ancestry","Sex","Age","BMI"),
       col = traits_cols,
       pch = 15,
       bty = 'n', ncol = 4)

d <- rbind.data.frame(sapply(c('EURv1','SEX2','AGE','BMI'), function(trait) 
  100*(sum(to_plot$trait == trait & to_plot$number == '1')/total_DEGs[trait])
),
sapply(c('EURv1','SEX2','AGE','BMI'), function(trait) 
  100*(sum(to_plot$trait == trait & to_plot$number %in% c('2','3'))/total_DEGs[trait])
),
sapply(c('EURv1','SEX2','AGE','BMI'), function(trait) 
  100*(sum(to_plot$trait == trait & to_plot$number %in% c('4','5','6','7','8'))/total_DEGs[trait])
))
colnames(d) <- c("Ancestry","Sex","Age","BMI")
rownames(d) <- c("Tissue-specific", "Moderate sharing", "High sharing")

library(reshape)
library(ggh4x)
library(RColorBrewer)
df2 <- melt(d)
df2$type <- rep(c("Tissue-specific", "Moderate sharing", "High sharing"), 4)
df2$type <- factor(df2$type, levels = rev(c("Tissue-specific", "Low sharing", "Moderate sharing", "High sharing")), order = T) 
cols <- brewer.pal(4, "Greys")[c(2:4)]
names(cols) <- c("Tissue-specific","Moderate sharing", "High sharing")

traits_cols <- c('#C49122','#4B8C61','#70A0DF','#A76595')
names(traits_cols) <- c("Ancestry","Sex","Age","BMI")

strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols[1]))
to_plot$trait <- factor(to_plot$trait, levels = c('EURv1','SEX2','AGE','BMI'))

p2 <- ggplot(df2[df2$variable == 'Ancestry',], aes(x = 1,
                                                   y = value,
                                                   fill = type)) +
  geom_bar(stat= "identity") + 
  coord_flip() +
  theme_bw() +
  facet_wrap2(~ variable, strip = strip, nrow = 1) + 
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  xlab("") + ylab("DMPs (%)")
p2

### Plot jitter for trait specific results #####
to_plot$type <- 'Hyper'
to_plot$type[to_plot$dir==-1] <- 'Hypo'

colors_traits <- list('AGE'=c('#3D7CD0','#B4D6F6'),
                      'SEX2'=c('#3B734E','#89AA94'),
                      'EURv1'=c('#F0AE21','#F9DE8B'))

to_plot$label <- 'No'
to_plot$label[to_plot$number >= 5] <- 'Yes'

library(ggrepel)
#pdf('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/tissue_sharing_sex_violin.pdf', width = 3, height = 3)
#png('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/tissue_sharing_ancestry_violin.png', width = 900, height = 700, res = 300)
p <- (ggplot(to_plot[to_plot$trait=='SEX2',]) + geom_jitter(aes(type, number, col=type), alpha=0.5) +
        theme_bw() + #geom_violin(aes(x = type, y = number, fill=type),alpha=0.8) +
        geom_boxplot(aes(type, number),col = "black",
                     fill = "white",
                     outlier.shape = NA,
                     notch = T,
                     width = 0.25) +
        geom_violin(aes(type, number,fill=type),col = "black")+
        scale_fill_manual(values = colors_traits[['SEX2']]) +
        scale_color_manual(values = colors_traits[['SEX2']]) + xlab('') + ylab('Nº Tissues') +
        ggrepel::geom_text_repel(aes(type, number,label=ifelse(label=='Yes',as.character(UCSC_RefGene_Name),'')), max.overlaps = Inf, )+
        theme(legend.position = "none",
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 14)))
#dev.off()
ggsave('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/tissue_sharing_sex_violin.png', p, dpi = 300, width = 2, height = 2, units = 'in')


# gene examples of highly shared genes


# Age --
cpgname <- "cg16867657"
betas <- lapply(tissues, function(tissue) readRDS(paste0(project_path,"Tissues/", tissue, "/data.rds"))[cpgname,])
names(betas) <- tissues

young_value <- sapply(tissues, function(tissue)
  median(as.numeric(betas[[tissue]][,metadata[[tissue]][as.numeric(metadata[[tissue]]$AGE) < 45, "SUBJID"]]))
)
old_value <- sapply(tissues, function(tissue)
  median(as.numeric(betas[[tissue]][,metadata[[tissue]][as.numeric(metadata[[tissue]]$AGE) >= 45, "SUBJID"]]))
)

data <- cbind.data.frame("Age" = c(rep("[20-45)",length(young_value)),
                                   rep("[45-70]",length(old_value))),
                         "value" = c(young_value, old_value),
                         "Tissue" = rep(tissues, 2))
data$Age <- factor(data$Age, levels = c("[20-45)", "[45-70]"), order = T)
data$Tissue <- factor(data$Tissue, levels = tissues,
                      order = T)
tissue_cols <- tissue_info[tissues, "colcodes"]
names(tissue_cols) <- as.character(levels(data$Tissue))
data$cpg <- rep(cpgname, nrow(data))
data$cpg <- factor(data$cpg)

p4 <- ggplot(data = data,
             aes(x = Age,
                 y = (value),
                 col = Tissue),
) +
  geom_violin(col = "black") +
  geom_boxplot(col = "black",
               fill = "white",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  geom_point(aes(col = Tissue),
             size = 0.8) +
  geom_line(aes(group=Tissue)) +
  xlab("") +
  ylab("Beta (median)") +
  scale_color_manual(values = tissue_cols) +
  labs(title="") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill=traits_cols["Age"]),
        legend.position = "none") +
  facet_grid(~cpg)

p4
library(ggpubr)
pdf(paste0(plot_path, "Figure_2E.tissue_sharing_examples.pdf"),
    width = 2, height = 8)
ggarrange(p3, p4, ncol = 2)
dev.off()


#### Enrichments ####
library(stringr)
genes <- to_plot_genes$UCSC_RefGene_Name[to_plot_genes$trait=='AGE' & to_plot_genes$number>=4]
#genes_non_conc <- unique(unlist(str_split(genes, ";")))[2:640]

write.table(genes_non_conc, '~/marenostrum/Projects/GTEx_v8/Methylation/genes_non_concordant_sharing_age.txt', sep = '\n', quote = F, 
            col.names = F, row.names = F)

library(missMethyl)

sharing <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Sharing_DMP.rds')
betas <- read.csv("~/marenostrum_scratch/MN4/bsc83/bsc83535/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")$Name
to_plot <- sharing
to_plot$number <- as.numeric(to_plot$number)
res <- gometh(to_plot$CG[to_plot$trait=='SEX2' & to_plot$number>=2 & to_plot$dir==1], all.cpg=(betas),
              collection="GO", array.type="EPIC")
res <- res[res$ONTOLOGY=="BP",]
print(table(res$FDR<0.05))

### in TSS ##
shared_cpgs <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/chromHMM_shared_9tissues.rds')
tss <- shared_cpgs[shared_cpgs$region_chromhmm=='TSS',]
res <- gometh(to_plot$CG[to_plot$trait=='SEX2' & to_plot$number>=2 & to_plot$dir==1],
              all.cpg=to_plot$CG[to_plot$trait=='SEX2'], collection="GO", array.type="EPIC", sig.genes = TRUE,
              genomic.features = c("TSS1500","TSS200","1stExon"))

write.table(res[res$FDR<0.05,], '~/marenostrum/Projects/GTEx_v8/Methylation/Data/SupplementaryTable4_enrichment_shared_age.hyper.txt', 
            sep = '\t', quote = F, 
            col.names = T, row.names = T)

topgo <- topGSA(res, n=20)
(ggplot(data = topgo[topgo$FDR<0.05,], aes(x=DE/N, y = factor(TERM, levels=rev(TERM)),
                                           color = -log10(FDR), size = DE)) +
    geom_point() + scale_color_gradient(low = "red", high = "blue") +
    theme_bw() +  ylab("") +  xlab("Gene Ratio") +
    ggtitle(paste0('AGE Shared'," GO BP sig")))

