#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Enrichment on TFBS for differentially methylated positions/loci
# @software version: R=4.2.2

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

tissues <- list.dirs("Tissues/", full.names = F)[-1]
tissues <- tissues[-grep('Old', tissues)]
hypo_e <- as.data.frame(matrix(nrow = length(tissues), ncol = 3, dimnames = list(tissues, c("AGE", "EURv1", "SEX2"))))
hyper_e <- as.data.frame(matrix(nrow = length(tissues), ncol = 3, dimnames = list(tissues, c("AGE", "EURv1", "SEX2"))))
hypo_d <- as.data.frame(matrix(nrow = length(tissues), ncol = 3, dimnames = list(tissues, c("AGE", "EURv1", "SEX2"))))
hyper_d <- as.data.frame(matrix(nrow = length(tissues), ncol = 3, dimnames = list(tissues, c("AGE", "EURv1", "SEX2"))))

n_samples <- c()
for(tissue in tissues){
  print(tissue)
  data <- read.csv(paste0("TFBS/tfbs_results_", tissue, "_tested.csv"))
  data$variable <- factor(data$variable, levels=c("AGE", "EURv1", "SEX2")) #For cases with 0 matches, to get a 0 in the output
  signif <- data[data$adj_p_val<0.05 & data$direction=="Enriched",]
  hyper_e[tissue,] <- tryCatch({table(signif$variable, signif$methylation)[,1]}, error = function(e) {return(c(0,0,0))})
  hypo_e[tissue,] <- tryCatch({table(signif$variable, signif$methylation)[,2]}, error = function(e) {return(c(0,0,0))})
  signif <- data[data$adj_p_val<0.05 & data$direction=="Depleted",]
  hyper_d[tissue,] <- tryCatch({table(signif$variable, signif$methylation)[,1]}, error = function(e) {return(c(0,0,0))})
  hypo_d[tissue,] <- tryCatch({table(signif$variable, signif$methylation)[,2]}, error = function(e) {return(c(0,0,0))})

  metadata <- readRDS(paste0("Tissues/",tissue, "/metadata.rds"))
  n_samples <- c(n_samples, nrow(metadata))  
}
names(n_samples) <- tissues


#Plot heatmaps
tissue_info <- readRDS("../Jose/00_Data/Tissue_info_whole.rds")
tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]
tissues_cols <- tissue_info[, 3]
names(tissues_cols) <- tissue_info[, 1]
order <- rev(sort(n_samples))
tissues_cols <- tissues_cols[names(order)]

row_ha_left <- HeatmapAnnotation("Samples" = anno_barplot(order,  #positives if we want to show the number of positives instead of n
                                                          gp = gpar(fill = tissues_cols,
                                                                    col = tissues_cols),
                                                          border=F, width = unit(1, "cm")),
                                 gap = unit(0.3,"cm"),
                                 show_legend = F,
                                 show_annotation_name = T,
                                 annotation_name_rot = 90,
                                 annotation_name_gp = gpar(fontsize = 10),
                                 which = "row")
create_heatmap <- function(data){
  colnames(data) <- c("Age", "Ancestry", "Sex")
  data <- data[names(order), c("Ancestry", "Sex", "Age")]
  h <- Heatmap(as.matrix(data), 
                col = colorRamp2( c(0, 10, 50, 100, 150),
                                  c("white", brewer.pal(9, "BuPu")[c(1, 4, 5, 7)])),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(as.matrix(data)[i, j], x, y, gp = gpar(fontsize = 9))},
                na_col = "white",
                cluster_rows = F,
                column_names_side = "top",
               column_names_rot =  60,
                cluster_columns = F,
                name = "Nº TFBS",
                # column_title_side = "bottom",
                left_annotation = row_ha_left,
                row_names_gp = gpar(fontsize = 10),
                column_names_gp = gpar(fontsize = 10),
                row_names_side = "left")
  return(h)
}
# pdf("Plots/TFBS/Hyper_enriched.pdf", height = 2.5, width = 4)
# create_heatmap(hyper_e)
# dev.off()
# pdf("Plots/TFBS/Hypo_enriched.pdf", height = 2.5, width = 4)
# create_heatmap(hypo_e)
# dev.off()
# pdf("Plots/TFBS/Hyper_depleted.pdf", height = 2.5, width = 4)
# create_heatmap(hyper_d)
# dev.off()
# pdf("Plots/TFBS/Hypo_depleted.pdf", height = 2.5, width = 4)
# create_heatmap(hypo_d)
# dev.off()


#Final format:
total <- hyper_e + hypo_e
pdf("Plots/TFBS/Total_enriched_tested.pdf", height = 3, width = 4)
create_heatmap(total)
dev.off()

total_d <- hyper_d + hypo_d
pdf("Plots/TFBS/Total_depleted_tested.pdf", height = 3, width = 4)
create_heatmap(total_d)
dev.off()

#Barplot separating hyper and hypo proportions
barplot <- function(hyper_e, hypo_e, total){
  library(ggh4x)
  traits_cols <- c('#C49122','#4B8C61','#70A0DF')
  names(traits_cols) <- c('Ancestry','Sex','Age')
  strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols))
  
  library(tidyr)
  library(dplyr)
  long_format_table_hyper <- replace(hyper_e/total, is.na(hyper_e/total), 0) %>%
    tibble::rownames_to_column(var = "Tissue") %>%
    gather(key = "Variable", value = "Value", -Tissue)
  long_format_table_hyper$Direction <- "Hyper"
  
  long_format_table_hypo <- replace(hypo_e/total, is.na(hypo_e/total), 0) %>%
    tibble::rownames_to_column(var = "Tissue") %>%
    gather(key = "Variable", value = "Value", -Tissue)
  long_format_table_hypo$Direction <- "Hypo"
  
  long_format_table <- rbind(long_format_table_hyper, long_format_table_hypo)
  long_format_table$Tissue <- factor(long_format_table$Tissue, levels=names(rev(order)))
  long_format_table$Variable <- ifelse(long_format_table$Variable=="AGE", "Age", ifelse(long_format_table$Variable=="SEX2", "Sex", "Ancestry"))
  long_format_table$Variable <- factor(long_format_table$Variable, levels = c('Ancestry','Sex','Age'))
  
  g <- ggplot(long_format_table, aes(fill=Direction, y=Tissue, x=Value)) + 
    geom_bar(position="stack", stat="identity", alpha=0.8) + ylab('')+
    scale_fill_manual(values = c('#b23c29','#1f6db6')) + xlab("Proportion") +
    facet_wrap2(~ Variable, strip = strip) + theme_bw() + scale_x_continuous(breaks=c(0, 0.5, 1), labels=c(0, 0.5, 1))
  return(g)
}
g <- barplot(hyper_e, hypo_e, total)
ggsave("Plots/TFBS/Total_enriched_direction_tested.pdf", g, height = 2.5, width = 4.5)

g <- barplot(hyper_d, hypo_d, total_d)
ggsave("Plots/TFBS/Total_depleted_direction_tested.pdf", g, height = 2.5, width = 4.5)


#Tissue sharing
final_table <- data.frame(TF="1", methylation=1, adj_p_val=1, variable=1, tissue=1)
for(trait in c("AGE", "EURv1", "SEX2")){
  print(trait)
  probes <- data.frame(TF="1",methylation=1, adj_p_val=1,variable=1, tissue=1)
  for(tissue in tissues){ 
    print(tissue)
    model <- read.csv(paste0("TFBS/tfbs_results_", tissue, "_tested.csv"))
    if (!trait %in% model$variable) {
      next}
    res <- model[model$variable==trait & model$adj_p_val<0.05 & model$direction=="Enriched",]
    if (nrow(res)<1) {
      next}
    probes <- rbind(probes, res[,c('TF','methylation','adj_p_val','variable','tissue')])
  }
  probes <- probes[-1,]
  final_table <- rbind(final_table, probes)
}
final_table <- final_table[-1,]

library(plyr) #Counting
counts <- ddply(final_table, .(TF, variable, methylation), nrow)
names(counts) <- c("TF", "Trait", "Methylation", "Number")

##### Final plot ####
traits_cols <- c('#C49122','#4B8C61','#70A0DF')
names(traits_cols) <- c('EURv1','SEX2','AGE')
strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols))
counts$Trait <- factor(counts$Trait, levels = c('EURv1','SEX2','AGE'))

g2 <- ggplot(counts, aes(y = Number, fill=Methylation)) +
  geom_bar(position = 'fill', alpha=0.8) + 
  geom_text(aes(label=after_stat(count), x = after_stat(count+1.5)), stat='count', position='fill', size=3,hjust=.8, angle=45) +
  theme_bw() + ylab('Nº of Tissues') + xlab('Proportion shared TFBS') +
  facet_wrap2(~ Trait, strip = strip, nrow = 1) + theme_bw() + 
  scale_x_continuous(breaks=seq(0, 1, 0.5), labels = c(0, 0.5, 1))

ggsave("Plots/TFBS/Tissue_sharing_tested.pdf", g2, height = 2, width = 4)

#Plot top shared
shared <- counts[counts$Number>1,]

### create table
tf_res <- list()
for(tissue in tissues){
  print(tissue)
  data <- read.csv(paste0("TFBS/tfbs_results_", tissue, "_tested.csv"))
  data$variable <- factor(data$variable, levels=c("AGE", "EURv1", "SEX2")) #For cases with 0 matches, to get a 0 in the output
  signif <- data[data$adj_p_val<0.05 & data$direction=="Enriched",] 
  tf_res[[tissue]] <- signif[signif$variable=='AGE',c(2:ncol(signif))]
  }
tf_res <- do.call('rbind.data.frame', tf_res)
write.table(tf_res[tf_res$TF!='Epitope',], '~/marenostrum/Projects/GTEx_v8/Methylation/Data/SupplementaryTable9_TFBS_enrichment_age.txt', sep = '\t', 
            col.names = T, row.names = F, quote = F)

#Plot top enriched TFBS in ggarrange
parsing_dummy_variable <- function(value) {
  return(gsub("_.*", "", value)) #Replace everything after the first _
}

for(tissue in tissues){
  print(tissue)
  data <- read.csv(paste0("TFBS/tfbs_results_", tissue, ".csv"))
  signif <- data[data$adj_p_val<0.05 & data$direction=="Enriched",]
  hyper_plot <- signif[signif$methylation=="Hyper",]
  hypo_plot <- signif[signif$methylation=="Hypo",]
  # plot <- signif[order(signif$adj_p_val),]
  plot <- rbind(hyper_plot, hypo_plot) #plot together the top 10 hyper and top 10 hypo?
  if(nrow(plot)==0){next}
  plot_top_10 <- plot %>%
    group_by(variable, methylation) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  #I create a dummy variable because the same TF can be in different variables and I want a different order in each facet
  plot_top_10$TF_dummy <- paste0(plot_top_10$TF, "_", plot_top_10$variable, "_", plot_top_10$methylation)
  plot_top_10$TF_dummy <- factor(plot_top_10$TF_dummy, levels = rev(plot_top_10$TF_dummy))

  plot_top_10$variable <- factor(plot_top_10$variable, levels = c('EURv1','SEX2','AGE'))
  
  trues <- names(traits_cols) %in% plot_top_10$variable
  strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols[trues]))
  
  g_plot <- ggplot(plot_top_10, aes(x=log(odds_ratio), y=TF_dummy, col=methylation)) +
    geom_errorbar(aes(xmin=log(CI_low), xmax=log(CI_high)), width=.4) +
    geom_vline(xintercept = 0) +
    geom_point(size=2.5) + ylab('') + theme_bw() +
    scale_colour_manual(values=c("#b23c29", "#1f6db6")) +
    xlab("Log(Odds ratio)") +
    theme(#legend.position = "none",
          axis.text.x = element_text(colour="black", size=11),
          axis.text.y = element_text(colour="black", size=13),
          legend.text = element_text(colour="black", size=12),
          axis.title.x = element_text(size=13),
          legend.spacing.y = unit(-0.05, "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", linewidth=1)) +
    facet_wrap2(~variable, scales = "free_y", strip = strip, nrow = 1) + #Function facet_wrap2 is different from facet_wrap only to color
    scale_y_discrete(labels = parsing_dummy_variable)
  ggsave(paste0("Plots/TFBS/Top_", tissue, ".pdf"), g_plot, width = 4 + sum(trues), height = 4)
}

### plot per trait
colors_traits <- list('AGE'=c('#3D7CD0','#B4D6F6'),
                      'SEX2'=c('#3B734E','#89AA94'),
                      'EURv1'=c('#F0AE21','#F9DE8B'))

shared <- counts[counts$Number>1,]
shared <- shared[shared$TF != 'Epitope',]

for(trait in c('AGE')){
  print(trait)
  data <- lapply(tissues, function(tissue) read.csv(paste0("TFBS/tfbs_results_", tissue, "_tested.csv")))
  data <- do.call('rbind.data.frame',data)
  data <- data[data$variable == trait,]
  signif <- data[data$adj_p_val<0.05 & data$direction=="Enriched",]
  hyper_plot <- signif[signif$methylation=="Hyper",]
  hypo_plot <- signif[signif$methylation=="Hypo",]
  #plot <- signif[order(signif$adj_p_val),]
  plot <- rbind(hyper_plot, hypo_plot) #plot together the top 10 hyper and top 10 hypo?
  plot <- signif[order(signif$adj_p_val),]
  plot <- plot[plot$TF!='Epitope',]
  if(nrow(plot)==0){next}
  plot_top_10 <- plot %>%
    group_by(tissue, methylation) %>%
    slice_head(n = 10) %>%
    ungroup()
  #plot_top_10 <- plot[plot$TF %in% shared$TF[shared$Trait==trait],]
  
  #I create a dummy variable because the same TF can be in different variables and I want a different order in each facet
  plot_top_10$TF_dummy <- paste0(plot_top_10$TF, "_", plot_top_10$tissue, "_", plot_top_10$methylation)
  plot_top_10$TF_dummy <- factor(plot_top_10$TF_dummy, levels = rev(plot_top_10$TF_dummy))
  
  #plot_top_10$tissue <- factor(plot_top_10$variable, levels = c('EURv1','SEX2','AGE'))
  
  trues <- names(traits_cols) %in% plot_top_10$tissue
  strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols[trues]))
  
  g_plot <- ggplot(plot_top_10, aes(x=log2(odds_ratio), y=TF_dummy, col=methylation)) +
    geom_errorbar(aes(xmin=log2(CI_low), xmax=log2(CI_high)), width=.4) +
    geom_vline(xintercept = 0) +
    geom_point(size=2.5) + ylab('') + theme_bw() +
    scale_colour_manual(values=colors_traits[[trait]]) +
    xlab("Log(Odds ratio)") +
    theme(#legend.position = "none",
      axis.text.x = element_text(colour="black", size=11),
      axis.text.y = element_text(colour="black", size=13),
      legend.text = element_text(colour="black", size=12),
      axis.title.x = element_text(size=13),
      legend.spacing.y = unit(-0.05, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", linewidth=1)) +
    facet_wrap2(~tissue, scales = "free_y", strip = strip, nrow = 1) + #Function facet_wrap2 is different from facet_wrap only to color
    scale_y_discrete(labels = parsing_dummy_variable)
  ggsave(paste0("Plots/TFBS/Top_", trait, "_tested.pdf"), g_plot, width = 9 + sum(trues), height = 4)
  
  library(reshape2)
  plot_top_10$logodds_ratio <- log2(plot_top_10$odds_ratio)
  wide <- dcast(plot_top_10[plot_top_10$methylation=='Hyper',c("tissue", "TF","logodds_ratio")], tissue ~ TF)
 # wide <- log2(wide)
  pal <- rev(c('#0061FF', '#1872FF', '#3183FF', '#4994FF', '#61A5FF', '#7AB7FF', '#92C8FF', '#AAD9FF', '#C3EAFF', '#DBFBFF'))
  rownames(wide) <- wide$tissue
  wide$tissue <- NULL
  pdf(paste0("Plots/TFBS/Top_", trait, "_tested_shared.Heatmap.pdf"), height = 3, width = 4)
  Heatmap(as.matrix(wide), 
               col = colorRamp2( c(0, 0.1, 0.6, 1.5, 2.5),
                                 c("white", pal[c(4, 6, 7, 10)])),
               # cell_fun = function(j, i, x, y, width, height, fill) {
               #   grid.text(as.matrix(data)[i, j], x, y, gp = gpar(fontsize = 9))},
               na_col = "white",
               cluster_rows = F,
               column_names_side = "top",
               column_names_rot =  60,
               cluster_columns = F,
               name = "log(OR)",
               # column_title_side = "bottom",
               #left_annotation = row_ha_left,
               row_names_gp = gpar(fontsize = 10),
               column_names_gp = gpar(fontsize = 10),
               row_names_side = "left")
  dev.off()
  
}

#Plot top enriched TFBS in ggarrange
parsing_dummy_variable <- function(value) {
  return(gsub(".*_", "", value)) #Replace everything after the first _
}


for(trait in c('SEX2')){
  print(trait)
  data <- lapply(tissues, function(tissue) read.csv(paste0("TFBS/tfbs_results_", tissue, "_tested.csv")))
  data <- do.call('rbind.data.frame',data)
  data <- data[data$variable == trait,]
  signif <- data[data$adj_p_val<0.05 & data$direction=="Enriched",]
  hyper_plot <- signif[signif$methylation=="Hyper",]
  hypo_plot <- signif[signif$methylation=="Hypo",]
  #plot <- signif[order(signif$adj_p_val),]
  plot <- rbind(hyper_plot, hypo_plot) #plot together the top 10 hyper and top 10 hypo?
  plot <- signif[order(signif$adj_p_val),]
  plot <- plot[plot$TF!='Epitope',]
  if(nrow(plot)==0){next}
  # plot_top_10 <- plot %>%
  #   group_by(tissue, methylation) %>%
  #   slice_head(n = 10) %>%
  #   ungroup()
  #plot_top_10 <- plot[plot$TF %in% shared$TF[shared$Trait==trait],]
  
  #I create a dummy variable because the same TF can be in different variables and I want a different order in each facet
  plot$TF_dummy <- paste0(plot$TF, "_", plot$tissue, "_", plot$methylation)
  plot$TF_dummy <- paste0(plot$tissue, "_", plot$methylation, "_", plot$TF)
  plot$TF_dummy <- factor(plot$TF_dummy, levels = plot$TF_dummy[order(plot$TF_dummy)])
  
  #plot_top_10$tissue <- factor(plot_top_10$variable, levels = c('EURv1','SEX2','AGE'))
  
  trues <- names(traits_cols) %in% plot$tissue
  strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols[trues]))
  
  g_plot <- ggplot(plot, aes(x=log2(odds_ratio), y=(TF_dummy), col=methylation)) +
    geom_errorbar(aes(xmin=log2(CI_low), xmax=log2(CI_high)), width=.4) +
    geom_vline(xintercept = 0) +
    geom_point(size=2.5) + ylab('') + theme_bw() +
    scale_colour_manual(values=colors_traits[[trait]]) +
    xlab("Log(Odds ratio)") +
    theme(#legend.position = "none",
      axis.text.x = element_text(colour="black", size=11),
      axis.text.y = element_text(colour="black", size=13),
      legend.text = element_text(colour="black", size=12),
      axis.title.x = element_text(size=13),
      legend.spacing.y = unit(-0.05, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", linewidth=1)) +
    facet_wrap2(~tissue, scales = "free_y", strip = strip, nrow = 1) + #Function facet_wrap2 is different from facet_wrap only to color
    scale_y_discrete(labels = parsing_dummy_variable)
  ggsave(paste0("Plots/TFBS/Top_", trait, "_tested.pdf"), g_plot, width = 6 + sum(trues), height = 4)

}


#Polycomb repressive complex genes:
PcG <- c("EZH2", "PCG1", "PCG2", "RYBP", "SUZ12", "RNF2", "CBX8", "KDM2B", "RING1B", "JARID2", "RYBP", "BMI1", "PHC1", "PHC2", "CBX", "RING1", "EED", "REST", "YY1")
hyper_plot$TF[hyper_plot$TF %in% PcG & hyper_plot$variable=="AGE"] #7 in lung



#Functional enrichment
library(WebGestaltR)
library(clusterProfiler)
library(org.Hs.eg.db)

for(tissue in tissues){
  print(tissue)
  data <- read.csv(paste0("TFBS/tfbs_results_", tissue, "_tested.csv"))
  bg <- unique(data$TF)
  signif <- data[data$adj_p_val<0.05 & data$direction=="Enriched",]
  if(nrow(signif)==0){next}
  for(variable in c("AGE", "EURv1", "SEX2")){
    signif_v <- signif[signif$variable==variable,]
    hyper <- signif_v$TF[signif_v$methylation=="Hyper"]
    hypo <- signif_v$TF[signif_v$methylation=="Hypo"]
    
    if(length(hyper)!=0){
      test <- enrichGO(gene= hyper, universe = bg,
                       keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                       ont = "BP")
      if(nrow(test@result[test@result$p.adjust<0.05,])>0){
        print(dotplot(test, title=paste(tissue, "-", variable,"- hypermethylation")))
      }
    }
    if(length(hypo)!=0){
      test <- enrichGO(gene= hypo, universe = bg,
                       keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                       ont = "BP")
      if(nrow(test@result[test@result$p.adjust<0.05,])>0){
        print(dotplot(test, title=paste(tissue, "-", variable,"- hypomethylation")))
      }
    }
    test <- enrichGO(gene= signif_v$TF, universe = bg,
                     keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                     ont = "BP")
    if(nrow(test@result[test@result$p.adjust<0.05,])>0){
      print(dotplot(test, title=paste(tissue, "-", variable,"- DMPs")))
    }
  }
}

## Quiescent
#Tissue sharing
final_table <- data.frame(TF="1", methylation=1, adj_p_val=1, variable=1, tissue=1)
for(trait in c("EURv1")){
  print(trait)
  probes <- data.frame(TF="1",methylation=1, adj_p_val=1,variable=1, tissue=1)
  for(tissue in tissues){ 
    print(tissue)
    model <- read.csv(paste0("TFBS/tfbs_results_", tissue, "_tested_Quies.csv"))
    if (!trait %in% model$variable) {
      next}
    res <- model[model$variable==trait & model$adj_p_val<0.05 & model$direction=="Enriched",]
    if (nrow(res)<1) {
      next}
    probes <- rbind(probes, res[,c('TF','methylation','adj_p_val','variable','tissue')])
  }
  probes <- probes[-1,]
  final_table <- rbind(final_table, probes)
}
final_table <- final_table[-1,]

library(plyr) #Counting
counts <- ddply(final_table, .(TF, variable, methylation), nrow)
names(counts) <- c("TF", "Trait", "Methylation", "Number")

##### Final plot ####
traits_cols <- c('#C49122','#4B8C61','#70A0DF')
names(traits_cols) <- c('EURv1','SEX2','AGE')
strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols))
counts$Trait <- factor(counts$Trait, levels = c('EURv1','SEX2','AGE'))

g2 <- ggplot(counts, aes(y = Number, fill=Methylation)) +
  geom_bar(position = 'fill', alpha=0.8) + 
  geom_text(aes(label=after_stat(count), x = after_stat(count+1.5)), stat='count', position='fill', size=3,hjust=.8, angle=45) +
  theme_bw() + ylab('Nº of Tissues') + xlab('Proportion shared TFBS') +
  facet_wrap2(~ Trait, strip = strip, nrow = 1) + theme_bw() + 
  scale_x_continuous(breaks=seq(0, 1, 0.5), labels = c(0, 0.5, 1))

ggsave("Plots/TFBS/Tissue_sharing_tested_quiescent.pdf", g2, height = 2, width = 4)

### enrichment quiescent tfs ####

for(tissue in tissues){
  print(tissue)
  data <- read.csv(paste0("TFBS/tfbs_results_", tissue, "_tested_Quies.csv"))
  bg <- unique(data$TF)
  signif <- data[data$adj_p_val<0.05 & data$direction=="Enriched",]
  if(nrow(signif)==0){next}
  for(variable in c("EURv1")){
    signif_v <- signif[signif$variable==variable,]
    hyper <- signif_v$TF[signif_v$methylation=="Hyper"]
    hypo <- signif_v$TF[signif_v$methylation=="Hypo"]
    
    if(length(hyper)!=0){
      test <- enrichGO(gene= hyper, universe = bg,
                       keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                       ont = "BP")
      if(nrow(test@result[test@result$p.adjust<0.05,])>0){
        print(dotplot(test, title=paste(tissue, "-", variable,"- hypermethylation")))
      }
    }
    if(length(hypo)!=0){
      test <- enrichGO(gene= hypo, universe = bg,
                       keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                       ont = "BP")
      if(nrow(test@result[test@result$p.adjust<0.05,])>0){
        print(dotplot(test, title=paste(tissue, "-", variable,"- hypomethylation")))
      }
    }
    test <- enrichGO(gene= signif_v$TF, universe = bg,
                     keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                     ont = "BP")
    if(nrow(test@result[test@result$p.adjust<0.05,])>0){
      print(dotplot(test, title=paste(tissue, "-", variable,"- DMPs")))
    }
  }
  
  ### enrichment reprpc tfs ####
  #Tissue sharing
  final_table <- data.frame(TF="1", methylation=1, adj_p_val=1, variable=1, tissue=1)
  for(trait in c("SEX2")){
    print(trait)
    probes <- data.frame(TF="1",methylation=1, adj_p_val=1,variable=1, tissue=1)
    for(tissue in tissues){ 
      print(tissue)
      model <- read.csv(paste0("TFBS/tfbs_results_", tissue, "_tested_ReprPC.csv"))
      if (!trait %in% model$variable) {
        next}
      res <- model[model$variable==trait & model$adj_p_val<0.05 & model$direction=="Enriched",]
      if (nrow(res)<1) {
        next}
      probes <- rbind(probes, res[,c('TF','methylation','adj_p_val','variable','tissue')])
    }
    probes <- probes[-1,]
    final_table <- rbind(final_table, probes)
  }
  final_table <- final_table[-1,]
  
  library(plyr) #Counting
  counts <- ddply(final_table, .(TF, variable, methylation), nrow)
  names(counts) <- c("TF", "Trait", "Methylation", "Number")
  
  ##### Final plot ####
  traits_cols <- c('#C49122','#4B8C61','#70A0DF')
  names(traits_cols) <- c('EURv1','SEX2','AGE')
  strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols))
  counts$Trait <- factor(counts$Trait, levels = c('EURv1','SEX2','AGE'))
  
  g2 <- ggplot(counts, aes(y = Number, fill=Methylation)) +
    geom_bar(position = 'fill', alpha=0.8) + 
    geom_text(aes(label=after_stat(count), x = after_stat(count+1.5)), stat='count', position='fill', size=3,hjust=.8, angle=45) +
    theme_bw() + ylab('Nº of Tissues') + xlab('Proportion shared TFBS') +
    facet_wrap2(~ Trait, strip = strip, nrow = 1) + theme_bw() + 
    scale_x_continuous(breaks=seq(0, 1, 0.5), labels = c(0, 0.5, 1))
  
  ggsave("Plots/TFBS/Tissue_sharing_tested_quiescent.pdf", g2, height = 2, width = 4)
  
  for(tissue in tissues){
    print(tissue)
    data <- read.csv(paste0("TFBS/tfbs_results_", tissue, "_tested_ReprPC.csv"))
    bg <- unique(data$TF)
    signif <- data[data$adj_p_val<0.05 & data$direction=="Enriched",]
    if(nrow(signif)==0){next}
    for(variable in c("SEX2")){
      signif_v <- signif[signif$variable==variable,]
      hyper <- signif_v$TF[signif_v$methylation=="Hyper"]
      hypo <- signif_v$TF[signif_v$methylation=="Hypo"]
      
      if(length(hyper)!=0){
        test <- enrichGO(gene= hyper, universe = bg,
                         keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                         ont = "BP")
        if(nrow(test@result[test@result$p.adjust<0.05,])>0){
          print(dotplot(test, title=paste(tissue, "-", variable,"- hypermethylation")))
        }
      }
      if(length(hypo)!=0){
        test <- enrichGO(gene= hypo, universe = bg,
                         keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                         ont = "BP")
        if(nrow(test@result[test@result$p.adjust<0.05,])>0){
          print(dotplot(test, title=paste(tissue, "-", variable,"- hypomethylation")))
        }
      }
      test <- enrichGO(gene= signif_v$TF, universe = bg,
                       keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                       ont = "BP")
      if(nrow(test@result[test@result$p.adjust<0.05,])>0){
        print(dotplot(test, title=paste(tissue, "-", variable,"- DMPs")))
      }
    }
  }
  
}

### TSS
#Tissue sharing
final_table <- data.frame(TF="1", methylation=1, adj_p_val=1, variable=1, tissue=1)
for(trait in c("SEX2")){
  print(trait)
  probes <- data.frame(TF="1",methylation=1, adj_p_val=1,variable=1, tissue=1)
  for(tissue in tissues){ 
    print(tissue)
    model <- read.csv(paste0("TFBS/tfbs_results_", tissue, "_tested_TSS.csv"))
    if (!trait %in% model$variable) {
      next}
    res <- model[model$variable==trait & model$adj_p_val<0.05 & model$direction=="Enriched",]
    if (nrow(res)<1) {
      next}
    probes <- rbind(probes, res[,c('TF','methylation','adj_p_val','variable','tissue')])
  }
  probes <- probes[-1,]
  final_table <- rbind(final_table, probes)
}
final_table <- final_table[-1,]

library(plyr) #Counting
counts <- ddply(final_table, .(TF, variable, methylation), nrow)
names(counts) <- c("TF", "Trait", "Methylation", "Number")

##### Final plot ####
traits_cols <- c('#C49122','#4B8C61','#70A0DF')
names(traits_cols) <- c('EURv1','SEX2','AGE')
strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols))
counts$Trait <- factor(counts$Trait, levels = c('EURv1','SEX2','AGE'))

g2 <- ggplot(counts, aes(y = Number, fill=Methylation)) +
  geom_bar(position = 'fill', alpha=0.8) + 
  geom_text(aes(label=after_stat(count), x = after_stat(count+1.5)), stat='count', position='fill', size=3,hjust=.8, angle=45) +
  theme_bw() + ylab('Nº of Tissues') + xlab('Proportion shared TFBS') +
  facet_wrap2(~ Trait, strip = strip, nrow = 1) + theme_bw() + 
  scale_x_continuous(breaks=seq(0, 1, 0.5), labels = c(0, 0.5, 1))

ggsave("Plots/TFBS/Tissue_sharing_tested_quiescent.pdf", g2, height = 2, width = 4)

for(tissue in tissues){
  print(tissue)
  data <- read.csv(paste0("TFBS/tfbs_results_", tissue, "_tested_ReprPC.csv"))
  bg <- unique(data$TF)
  signif <- data[data$adj_p_val<0.05 & data$direction=="Enriched",]
  if(nrow(signif)==0){next}
  for(variable in c("SEX2")){
    signif_v <- signif[signif$variable==variable,]
    hyper <- signif_v$TF[signif_v$methylation=="Hyper"]
    hypo <- signif_v$TF[signif_v$methylation=="Hypo"]
    
    if(length(hyper)!=0){
      test <- enrichGO(gene= hyper, universe = bg,
                       keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                       ont = "BP")
      if(nrow(test@result[test@result$p.adjust<0.05,])>0){
        print(dotplot(test, title=paste(tissue, "-", variable,"- hypermethylation")))
      }
    }
    if(length(hypo)!=0){
      test <- enrichGO(gene= hypo, universe = bg,
                       keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                       ont = "BP")
      if(nrow(test@result[test@result$p.adjust<0.05,])>0){
        print(dotplot(test, title=paste(tissue, "-", variable,"- hypomethylation")))
      }
    }
    test <- enrichGO(gene= signif_v$TF, universe = bg,
                     keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                     ont = "BP")
    if(nrow(test@result[test@result$p.adjust<0.05,])>0){
      print(dotplot(test, title=paste(tissue, "-", variable,"- DMPs")))
    }
  }
}



## shared 
#Tissue sharing
final_table <- data.frame(TF="1", methylation=1, adj_p_val=1, variable=1, tissue=1)
for(trait in c("SEX2")){
  print(trait)
  probes <- data.frame(TF="1",methylation=1, adj_p_val=1,variable=1, tissue=1)
  for(tissue in tissues){ 
    print(tissue)
    model <- read.csv(paste0("TFBS/tfbs_results_", tissue, "_tested_shared.csv"))
    if (!trait %in% model$variable) {
      next}
    res <- model[model$variable==trait & model$adj_p_val<0.05 & model$direction=="Enriched",]
    if (nrow(res)<1) {
      next}
    probes <- rbind(probes, res[,c('TF','methylation','adj_p_val','variable','tissue')])
  }
  probes <- probes[-1,]
  final_table <- rbind(final_table, probes)
}
final_table <- final_table[-1,]

library(plyr) #Counting
counts <- ddply(final_table, .(TF, variable, methylation), nrow)
names(counts) <- c("TF", "Trait", "Methylation", "Number")

##### Final plot ####
traits_cols <- c('#C49122','#4B8C61','#70A0DF')
names(traits_cols) <- c('EURv1','SEX2','AGE')
strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols[2]))
counts$Trait <- factor(counts$Trait, levels = c('EURv1','SEX2','AGE'))

g2 <- ggplot(counts, aes(y = Number, fill=Methylation)) +
  geom_bar(position = 'fill', alpha=0.8) + 
  geom_text(aes(label=after_stat(count), x = after_stat(count+1.5)), stat='count', position='fill', size=3,hjust=.8, angle=45) +
  theme_bw() + ylab('Nº of Tissues') + xlab('Proportion shared TFBS') +
  facet_wrap2(~ Trait, strip = strip, nrow = 1) + theme_bw() + 
  scale_x_continuous(breaks=seq(0, 1, 0.5), labels = c(0, 0.5, 1))

ggsave("Plots/TFBS/Tissue_sharing_tested_TSS.pdf", g2, height = 2, width = 4)

for(tissue in tissues){
  print(tissue)
  data <- read.csv(paste0("TFBS/tfbs_results_", tissue, "_tested_shared.csv"))
  bg <- unique(data$TF)
  signif <- data[data$adj_p_val<0.05 & data$direction=="Enriched",]
  if(nrow(signif)==0){next}
  for(variable in c("SEX2")){
    signif_v <- signif[signif$variable==variable,]
    hyper <- signif_v$TF[signif_v$methylation=="Hyper"]
    hypo <- signif_v$TF[signif_v$methylation=="Hypo"]
    
    if(length(hyper)!=0){
      test <- enrichGO(gene= hyper, universe = bg,
                       keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                       ont = "BP")
      if(nrow(test@result[test@result$p.adjust<0.05,])>0){
        print(dotplot(test, title=paste(tissue, "-", variable,"- hypermethylation")))
      }
    }
    if(length(hypo)!=0){
      test <- enrichGO(gene= hypo, universe = bg,
                       keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                       ont = "BP")
      if(nrow(test@result[test@result$p.adjust<0.05,])>0){
        print(dotplot(test, title=paste(tissue, "-", variable,"- hypomethylation")))
      }
    }
    test <- enrichGO(gene= signif_v$TF, universe = bg,
                     keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                     ont = "BP")
    if(nrow(test@result[test@result$p.adjust<0.05,])>0){
      print(dotplot(test, title=paste(tissue, "-", variable,"- DMPs")))
    }
  }
}

### hyper_tfbs 
lung_tfbs <- read.csv(paste0("TFBS/tfbs_results_", 'Lung', "_tested_shared.csv"))
colon_tfbs <- read.csv(paste0("TFBS/tfbs_results_", 'ColonTransverse', "_tested_shared.csv"))

bg <- unique(c(lung_tfbs$TF, colon_tfbs$TF))
hyper_lung <- lung_tfbs[lung_tfbs$adj_p_val<0.05 & lung_tfbs$direction=="Enriched" & lung_tfbs$variable=='SEX2' & lung_tfbs$methylation=='Hyper',]
hyper_colon <- colon_tfbs[colon_tfbs$adj_p_val<0.05 & colon_tfbs$direction=="Enriched" & colon_tfbs$variable=='SEX2' & colon_tfbs$methylation=='Hyper',]

hyper <- unique(c(hyper_colon$TF, hyper_lung$TF))

write.table(hyper, '~/marenostrum/Projects/GTEx_v8/Methylation/Data/TFBS_enriched_shared_hyper_sex.txt', sep = '\n', 
            col.names = F, row.names = F, quote = F)
write.table(bg, '~/marenostrum/Projects/GTEx_v8/Methylation/Data/TFBS_bg_hyper_sex.txt', sep = '\n', 
            col.names = F, row.names = F, quote = F)

test <- enrichGO(gene= hyper, universe = bg,
                 keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                 ont = "CC")
nrow(test@result[test@result$p.adjust<0.05,])
View(test@result)

test <- enrichKEGG(gene= hyper, universe = bg,
                 keyType = "SYMBOL")

### shared all
lung_tfbs <- read.csv(paste0("TFBS/tfbs_results_", 'shared', "_tested_shared.csv"))
bg <- unique(c(lung_tfbs$TF))
lung_tfbs <- lung_tfbs[lung_tfbs$adj_p_val<0.05 & lung_tfbs$direction=='Enriched',]
write.table(lung_tfbs, '~/marenostrum/Projects/GTEx_v8/Methylation/Data/enrichment_tfbs_shared_sex.txt',sep = '\t',
            quote = F, row.names = F, col.names = T)
hyper <- lung_tfbs[lung_tfbs$adj_p_val<0.05 & lung_tfbs$direction=="Enriched" & lung_tfbs$variable=='SEX2' & lung_tfbs$methylation=='Hyper',]
hypo <- lung_tfbs[lung_tfbs$adj_p_val<0.05 & lung_tfbs$direction=="Enriched" & lung_tfbs$variable=='SEX2' & lung_tfbs$methylation=='Hypo',]

test <- enrichGO(gene= hyper$TF, universe = bg,
                 keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                 ont = "BP")
test_hypo <- enrichGO(gene= hypo$TF, universe = bg,
                 keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                 ont = "BP")
nrow(test@result[test@result$p.adjust<0.05,])
View(test_hypo@result)

write.table(unique(hyper$TF), '~/marenostrum/Projects/GTEx_v8/Methylation/Data/TFBS_enriched_shared_hyper_sex.txt', sep = '\n', 
            col.names = F, row.names = F, quote = F)
write.table(unique(hypo$TF), '~/marenostrum/Projects/GTEx_v8/Methylation/Data/TFBS_enriched_shared_hypo_sex.txt', sep = '\n', 
            col.names = F, row.names = F, quote = F)
write.table(bg, '~/marenostrum/Projects/GTEx_v8/Methylation/Data/TFBS_bg_hyper_sex.txt', sep = '\n', 
            col.names = F, row.names = F, quote = F)

#Polycomb repressive complex genes:
PcG <- c("EZH2", "PCG1", "PCG2", "RYBP", "SUZ12", "RNF2", "CBX8", "KDM2B", "RING1B", "JARID2", "RYBP", "BMI1", "PHC1", "PHC2", "CBX", "RING1", "EED", "REST", "YY1")
hyper$TF[hyper$TF %in% PcG & hyper$variable=="SEX2"] #7 in lung
hypo$TF[hypo$TF %in% PcG & hypo$variable=="SEX2"] #7 in lung


parsing_dummy_variable <- function(value) {
  return(gsub(".*_", "", value)) #Replace everything after the first _
}


data <- read.csv(paste0("TFBS/tfbs_results_shared_tested_shared.csv"))
#data <- do.call('rbind.data.frame',data)
data <- data[data$variable == 'SEX2',]
signif <- data[data$adj_p_val<0.05 & data$direction=="Enriched",]
hyper_plot <- signif[signif$methylation=="Hyper",]
hypo_plot <- signif[signif$methylation=="Hypo",]
#plot <- signif[order(signif$adj_p_val),]
plot <- rbind(hyper_plot, hypo_plot) #plot together the top 10 hyper and top 10 hypo?
plot <- signif[order(signif$odds_ratio, decreasing = T),]
plot <- plot[plot$TF!='Epitope',]
if(nrow(plot)==0){next}
plot_top_10 <- plot %>%
  group_by(methylation) %>%
  slice_head(n = 15) %>%
  ungroup()
#plot_top_10 <- plot[plot$TF %in% shared$TF[shared$Trait==trait],]

#I create a dummy variable because the same TF can be in different variables and I want a different order in each facet
plot_top_10$TF_dummy <- paste0(plot_top_10$TF, "_", plot_top_10$tissue, "_", plot_top_10$methylation)
plot_top_10$TF_dummy <- paste0(plot_top_10$tissue, "_", plot_top_10$methylation, "_", plot_top_10$TF)
plot_top_10$TF_dummy <- factor(plot_top_10$TF_dummy, levels = plot_top_10$TF_dummy[order(plot_top_10$TF_dummy)])

#plot_top_10$tissue <- factor(plot_top_10$variable, levels = c('EURv1','SEX2','AGE'))
colors_traits <- list('AGE'=c('#3D7CD0','#B4D6F6'),
                      'SEX2'=c('#3B734E','#89AA94'),
                      'EURv1'=c('#F0AE21','#F9DE8B'))


g_plot <- ggplot(plot_top_10, aes(x=log2(odds_ratio), y=(TF_dummy), col=methylation)) +
geom_errorbar(aes(xmin=log2(CI_low), xmax=log2(CI_high)), width=.4) +
geom_vline(xintercept = 0) +
geom_point(size=2.5) + ylab('') + theme_bw() +
scale_colour_manual(values=colors_traits[[trait]]) +
xlab("Log(Odds ratio)") +
theme(#legend.position = "none",
  axis.text.x = element_text(colour="black", size=11),
  axis.text.y = element_text(colour="black", size=13),
  legend.text = element_text(colour="black", size=12),
  axis.title.x = element_text(size=13),
  legend.spacing.y = unit(-0.05, "cm"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(colour = "black", linewidth=1)) +
facet_wrap2(~tissue, scales = "free_y", strip = strip, nrow = 1) + #Function facet_wrap2 is different from facet_wrap only to color
scale_y_discrete(labels = parsing_dummy_variable)
ggsave(paste0("Plots/TFBS/Top_", trait, "_shared_tested.pdf"), g_plot, width = 4 + sum(trues), height = 6)


xist <- read.delim('~/marenostrum/Projects/GTEx_v8/Methylation/Data/xist_binding_proteins_suppl3_chu_2015.txt', header = F)
xist$V1 <- toupper(xist$V1)

hyper$TF[hyper$TF %in% xist$V1 & hyper$variable=="SEX2"] #7 in lung

### are xist genes and/or polycomb DE?
deg <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/DEA_GTEx.rds')
tissues

for (tissue in tissues) {
  print(tissue)
  dea <- deg[[tissue]][['Sex']]
  print(dea$gene_name[dea$gene_name %in% c(xist$V1, PcG) & dea$adj.P.Val<0.05])
}

dea <- lapply(tissues, function(tissue) deg[[tissue]][['Sex']][deg[[tissue]][['Sex']]$gene_name %in% c(xist$V1, PcG) & deg[[tissue]][['Sex']]$adj.P.Val<0.05,])
names(dea) <- tissues
dea <- do.call('rbind.data.frame', dea)

dnmt <- c('DNMT1' , 'DNMT3A', 'DNMT3B')
demethylases <- c('TET1','TET2','TET3','TDG')
PcG <- c("EZH2", "PCG1", "PCG2", "RYBP", "SUZ12", "RNF2", "CBX8", "KDM2B", "RING1B", "JARID2", "RYBP", "BMI1", "PHC1", "PHC2", "CBX", "RING1", "EED", "REST", "YY1")

tis <- names(deg)
tis <- tis[!tis %in% c('Ovary','Uterus','Testis','Prostate','Vagina')]
dea <- lapply(tis, function(tissue) deg[[tissue]][['Sex']][deg[[tissue]][['Sex']]$gene_name %in% c(xist$V1, PcG, dnmt, demethylases) & deg[[tissue]][['Sex']]$adj.P.Val<0.05,])
names(dea) <- tis
dea <- do.call('rbind.data.frame', dea)
View(as.data.frame(table(dea$gene_name)))

### plot grant ####
parsing_dummy_variable <- function(value) {
  return(gsub(".*_", "", value)) #Replace everything after the first _
}


#data <- read.csv(paste0("TFBS/tfbs_results_WholeBlood_tested_grant_dmpsbg.csv"))
data <- do.call('rbind.data.frame',data)
data <- data[data$variable == 'SEX2',]
signif <- data[data$adj_p_val<0.05 & data$direction=="Enriched",]
hyper_plot <- signif[signif$methylation=="Hyper",]
hypo_plot <- signif[signif$methylation=="Hypo",]
#plot <- signif[order(signif$adj_p_val),]
plot <- rbind(hyper_plot, hypo_plot) #plot together the top 10 hyper and top 10 hypo?
plot <- signif[order(signif$odds_ratio, decreasing = T),]
plot <- plot[plot$TF!='Epitope',]
if(nrow(plot)==0){next}
plot_top_10 <- plot %>%
  group_by(methylation) %>%
  slice_head(n = 15) %>%
  ungroup()
#plot_top_10 <- plot[plot$TF %in% shared$TF[shared$Trait==trait],]

plot_top_10 <- plot

#I create a dummy variable because the same TF can be in different variables and I want a different order in each facet
plot_top_10$TF_dummy <- paste0(plot_top_10$TF, "_", plot_top_10$tissue, "_", plot_top_10$methylation)
plot_top_10$TF_dummy <- paste0(plot_top_10$tissue, "_", plot_top_10$methylation, "_", plot_top_10$TF)
plot_top_10$TF_dummy <- factor(plot_top_10$TF_dummy, levels = plot_top_10$TF_dummy[order(plot_top_10$TF_dummy)])

#plot_top_10$tissue <- factor(plot_top_10$variable, levels = c('EURv1','SEX2','AGE'))
colors_traits <- list('AGE'=c('#3D7CD0','#B4D6F6'),
                      'SEX2'=c('#3B734E','#89AA94'),
                      'EURv1'=c('#F0AE21','#F9DE8B'))


g_plot <- ggplot(plot_top_10, aes(x=log2(odds_ratio), y=(TF_dummy), col=methylation)) +
  geom_errorbar(aes(xmin=log2(CI_low), xmax=log2(CI_high)), width=.4) +
  geom_vline(xintercept = 0) +
  geom_point(size=2.5) + ylab('') + theme_bw() +
  scale_colour_manual(values=colors_traits[[trait]]) +
  xlab("Log(Odds ratio)") +
  theme(#legend.position = "none",
    axis.text.x = element_text(colour="black", size=11),
    axis.text.y = element_text(colour="black", size=13),
    legend.text = element_text(colour="black", size=12),
    axis.title.x = element_text(size=13),
    legend.spacing.y = unit(-0.05, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", linewidth=1)) +
  facet_wrap2(~tissue, scales = "free_y", strip = strip, nrow = 1) + #Function facet_wrap2 is different from facet_wrap only to color
  scale_y_discrete(labels = parsing_dummy_variable)
ggsave(paste0("Plots/TFBS/All_", trait, "_grant_tested.pdf"), g_plot, width = 4 , height = 15)

