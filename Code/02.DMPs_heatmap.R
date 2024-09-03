#This code creates a heatmap to compare tissues and demographic traits
rm(list=ls())
suppressPackageStartupMessages(library(ComplexHeatmap))
library(RColorBrewer)
suppressPackageStartupMessages(library(circlize))

first_dir <- "~/Documents/mn4/"
first_dir <- "~/marenostrum/"

project_path <- paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Methylation/")
project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "KidneyCortex", "Testis", "WholeBlood", "MuscleSkeletal")
names <- c("Age", "Ancestry", "BMI", "Sex")
# names <- c("Age", "Ancestry", "Sex")
data <- matrix(nrow = length(tissues), ncol=4, dimnames = list(tissues, names))
# data <- matrix(nrow = length(tissues), ncol=3, dimnames = list(tissues, names))

tissue_info <- readRDS(paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Jose/00_Data/Tissue_info_whole.rds"))
tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Methylation/Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]
tissue_info$name <- gsub("- ", "", tissue_info$tissue_name)

n_samples <- c()
for(tissue in tissues){ 
  print(tissue)
  # model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results.rds"))
  model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.rds"))
  # model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_batch.rds"))
  # model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_batch_5_PEERs.rds"))
  # model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_no_BMI.rds"))
  # model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DMR_results_5_PEERs.rds")) #Jose used the wrong name to the file DMR instead of DML
  data[tissue, "Age"] <- sum(model$AGE$adj.P.Val<0.05)
  if(TRUE %in% grepl("EURv1",names(model))){
    data[tissue, "Ancestry"] <- sum(model$EURv1$adj.P.Val<0.05)
  } else{data[tissue, "Ancestry"] <- NA}
  data[tissue, "BMI"] <- sum(model$BMI$adj.P.Val<0.05)
  if(TRUE %in% grepl("SEX",names(model))){
    data[tissue, "Sex"] <- sum(model$SEX2$adj.P.Val<0.05)
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

data <- data[,c("Ancestry", "Sex", "Age", "BMI")]
# data <- data[,c("Ancestry", "Sex", "Age")]
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
          na_col = "white",
          cluster_rows = F,
          cluster_columns = F,
          name = "#DVPs",
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
rownames(data) <- tissue_info$name[match(rownames(data), tissue_info$tissue_ID)]
pdf(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/Heatmap.pdf"), height = 4.5, width = 6.5)
create_heatmap(data, tissue_info, size=12)
dev.off()

#Tissue Sharing
final_table <- data.frame(probe="1", tissue="1", logFC=1, trait="1")
for(trait in c("AGE", "AncestryEUR", "BMI", "SEX2")){
  print(trait)
  probes <- data.frame(probe="1", tissue="1", logFC=1, trait="1")
  for(tissue in tissues){ 
    print(tissue)
    model <- readRDS(paste0(project_path, tissue, "/DML_results_no_smoking.rds"))
    if(is.null(nrow(model[[trait]]))){
      next
    }
    subset <- model[[trait]][model[[trait]]$adj.P.Val<0.05,]
    if(nrow(subset)==0){
      next
    }
    data_frame <- cbind("probe"=rownames(subset), "tissue"=tissue, "logFC"=subset[,1], "trait"=trait)
    probes <- rbind(probes, data_frame)
  }
  probes <- probes[-1,]
  final_table <- rbind(final_table, probes)
}
final_table <- final_table[-1,]

library(plyr) #Counting
counts <- ddply(final_table, .(probe, trait), nrow)
names(counts) <- c("probe", "trait", "number")
final_table$logFC <- as.numeric(final_table$logFC)
logfc <- ddply(final_table, .(probe, trait), summarise, logFC=mean(logFC))


nrow(counts[counts$number>1,])
nrow(counts[counts$number>2,])
nrow(counts[counts$number>3,])
nrow(counts[counts$number>4,])
nrow(counts[counts$number>5,])
nrow(counts[counts$number>6,])

counts <- counts[counts$number>=2,]

to_plot <- merge(counts, logfc, by=c("probe", "trait"))
library(ggplot2)
ggplot(to_plot) + geom_jitter(aes(trait, number, col=logFC), alpha=0.5) +
  scale_color_gradient2(low="red", mid="gray", high="blue") + theme_bw()

nrow(to_plot[to_plot$trait=="AncestryEUR" & to_plot$logFC>0,])
nrow(to_plot[to_plot$trait=="AncestryEUR" & to_plot$logFC<0,])

nrow(to_plot[to_plot$trait=="SEX2" & to_plot$logFC>0,])
nrow(to_plot[to_plot$trait=="SEX2" & to_plot$logFC<0,])

nrow(to_plot[to_plot$trait=="AGE" & to_plot$logFC>0,])
nrow(to_plot[to_plot$trait=="AGE" & to_plot$logFC<0,])
