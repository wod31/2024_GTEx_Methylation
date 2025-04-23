#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Script to plot a heatmap with the DVGs (genes) and their direction
# @software version: R=4.2.2

rm(list=ls())
suppressPackageStartupMessages(library(ComplexHeatmap))
library(RColorBrewer)
suppressPackageStartupMessages(library(circlize))
library(plyr)
library(reshape)
library(tidyr)
library(dplyr)
library(ggh4x)

first_dir <- "~/Documents/mn5/"
first_dir <- "~/marenostrum/"

project_path <- paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Methylation/")
project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "KidneyCortex", "Testis", "WholeBlood", "MuscleSkeletal")
names <- c("Age", "Ancestry", "Sex")
data <- matrix(nrow = length(tissues), ncol=3, dimnames = list(tissues, names))

# tissue_info <- readRDS(paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Jose/00_Data/Tissue_info_whole.rds"))
tissue_info <- readRDS(paste0(project_path, "Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]
tissue_info$name <- gsub("- ", "", tissue_info$tissue_name)

n_samples <- c()
for(tissue in tissues){ 
  print(tissue)
  model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DVG_Age.rds"))
  model <- model[!is.na(model$Adj.P.Value),]
  data[tissue, "Age"] <- sum(model$Adj.P.Value<0.05)
  model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DVG_Ancestry.rds"))
  model <- model[!is.na(model$Adj.P.Value),]
  data[tissue, "Ancestry"] <- sum(model$Adj.P.Value<0.05)
  if(file.exists(paste0(project_path, "Tissues/",tissue, "/DVG_Sex.rds"))){
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DVG_Sex.rds"))
    model <- model[!is.na(model$Adj.P.Value),]
    data[tissue, "Sex"] <- sum(model$Adj.P.Value<0.05)
  } 
  
  metadata <- readRDS(paste0(project_path, "Tissues/",tissue, "/metadata_expression.rds"))
  admixture_ancestry <- read.table('~/marenostrum/Projects/GTEx_v8/Methylation/Data/admixture_inferred_ancestry.txt')
  colnames(admixture_ancestry) <- c('Donor','AFRv1','EURv1','inferred_ancestry','AFRv2','EURv2')
  metadata <- merge(metadata, admixture_ancestry[,c("Donor","EURv1")], by='Donor')
  n_samples <- c(n_samples, nrow(metadata))
}



data_up <- matrix(nrow = length(tissues), ncol=3, dimnames = list(tissues, names))
for(tissue in tissues){ 
  print(tissue)
  
  model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DVG_Age.rds"))
  model <- model[!is.na(model$Adj.P.Value),]
  data_up[tissue, "Age"] <- sum(model$Adj.P.Value<0.05 & model$DiffLevene>0)
  model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DVG_Ancestry.rds"))
  model <- model[!is.na(model$Adj.P.Value),]
  data_up[tissue, "Ancestry"] <- sum(model$Adj.P.Value<0.05 & model$DiffLevene>0)
  if(file.exists(paste0(project_path, "Tissues/",tissue, "/DVG_Sex.rds"))){
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DVG_Sex.rds"))
    model <- model[!is.na(model$Adj.P.Value),]
    data_up[tissue, "Sex"] <- sum(model$Adj.P.Value<0.05 & model$DiffLevene>0)
  }
  
}

data_down <- matrix(nrow = length(tissues), ncol=3, dimnames = list(tissues, names))
for(tissue in tissues){ 
  print(tissue)
  
  model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DVG_Age.rds"))
  model <- model[!is.na(model$Adj.P.Value),]
  data_down[tissue, "Age"] <- sum(model$Adj.P.Value<0.05 & model$DiffLevene<0)
  model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DVG_Ancestry.rds"))
  model <- model[!is.na(model$Adj.P.Value),]
  data_down[tissue, "Ancestry"] <- sum(model$Adj.P.Value<0.05 & model$DiffLevene<0)
  if(file.exists(paste0(project_path, "Tissues/",tissue, "/DVG_Sex.rds"))){
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DVG_Sex.rds"))
    model <- model[!is.na(model$Adj.P.Value),]
    data_down[tissue, "Sex"] <- sum(model$Adj.P.Value<0.05 & model$DiffLevene<0)
  }
  
}

my_pretty_num_function <- function(n){
  if(n==""){
    return(n)
  } else{
    prettyNum(n, big.mark = ",")
  }
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
data <- data[,c("Ancestry", "Sex", "Age")]
rownames(data) <- tissue_info$name[match(rownames(data), tissue_info$tissue_ID)]
pdf(paste0(project_path, "Plots/Heatmap_DVGs.pdf"), height = 4.5, width = 6.5)
create_heatmap(data, tissue_info, size=12)
dev.off()
# create_heatmap(data_up, tissue_info, size=12)
# create_heatmap(data_down, tissue_info, size=12)

#### number hyper hypo barplot ####
sex_tissues <- c('Ovary','Prostate','Testis')
names <- c("Ancestry", "Sex", "Age")

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

binom_res <- lapply(names, function(trait) lapply(tissues, function(tissue) binom_test(trait,tissue)))
names(binom_res) <- names

for (trait in names) {
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
    tissues_to_test <- tissues
  }
  if (trait == 'Ancestry') {
    tissues_to_test <- tissues[tissues != 'KidneyCortex']
  }
  binom_all[[trait]] <- read_data(tissues_to_test, binom_res, trait)
}

# bar plot
traits_cols <- c('#C49122','#4B8C61','#70A0DF')
names(traits_cols) <- c('Ancestry','Sex','Age')
strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols))

tissue_expression_variation_explained_df <- tissue_expression_variation_explained_df[tissue_expression_variation_explained_df$Trait!='BMI',]
tissue_expression_variation_explained_df$Trait <- factor(tissue_expression_variation_explained_df$Trait, levels = c('Ancestry','Sex','Age'))
tissue_expression_variation_explained_df$Tissue <- factor(tissue_expression_variation_explained_df$Tissue, 
                                                          levels = rev(tissues))
#tissue_expression_variation_explained_df$Trait <- droplevels(tissue_expression_variation_explained_df$Trait)
ggplot(tissue_expression_variation_explained_df, aes(fill=Direction, y=Tissue, x=Prop)) + 
  geom_bar(position="stack", stat="identity", alpha=0.8) + ylab('')+
  scale_fill_manual(values = c('#b23c29','#1f6db6')) + 
  facet_wrap2(~ Trait, strip = strip) + theme_bw() + scale_x_continuous(breaks=seq(0, 1, 0.5))

tissue_expression_variation_explained_df$label <- paste0(tissue_expression_variation_explained_df$Direction,
                                                         ':',tissue_expression_variation_explained_df$Trait)

## new colors 
colors_traits <- list('Hyper:Age'=c('#3D7CD0'),'Hypo:Age'=c('#B4D6F6'),
                      'Hyper:Sex'=c('#3B734E'),'Hypo:Sex'=c('#89AA94'),
                      'Hyper:Ancestry'=c('#F0AE21'),'Hypo:Ancestry'=c('#F9DE8B'))

tissue_expression_variation_explained_df <- tissue_expression_variation_explained_df[!is.na(tissue_expression_variation_explained_df$Trait),]

pdf('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/DVG.barplot.hyper.hypo.new.colors.pdf', height = 4, width = 6)
ggplot(tissue_expression_variation_explained_df, aes(fill=label, y=Tissue, x=Prop)) + 
  geom_bar(position="stack", stat="identity", alpha=0.8) + ylab('')+
  scale_fill_manual(values = colors_traits, labels = c("Old", "EA","Female",'Young','AA','Male')) + xlab('Proportion DMPs')+
  facet_wrap2(~ Trait, strip = strip) + theme_bw() + scale_x_continuous(breaks=seq(0, 1, 0.5))
dev.off()


#Tissue sharing

final_table <- data.frame(gene="gene1",DiffLevene=1, Adj.P.Value=1, trait=1, tissue=1)
for(trait in c("Age", "Ancestry", "Sex")){
  print(trait)
  probes <- data.frame(gene="gene1",DiffLevene=1, Adj.P.Value=1,trait=1, tissue=1)
  for(tissue in tissues){ 
    print(tissue)
    if(!file.exists(paste0(project_path, "Tissues/",tissue, "/DVG_", trait, ".rds"))){next}
    res <- readRDS(paste0(project_path, "Tissues/",tissue, "/DVG_", trait, ".rds"))
    if (nrow(res)<1) { next}
    res$tissue <- tissue
    res$trait <- trait
    res$sign <- sign(res$DiffLevene)
    res$gene <- rownames(res)
    probes <- rbind(probes, res[,c('gene','DiffLevene','Adj.P.Value','trait','tissue')])
  }
  probes <- probes[-1,]
  final_table <- rbind(final_table, probes)
}
final_table <- final_table[-1,]

### add gene name ###
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

genes <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","chromosome_name"),
                               values= gsub('\\..*','', final_table$gene),mart= mart)
final_table$gene <- gsub('\\..*','', final_table$gene)
final_table <- merge(final_table, genes, by.x='gene', by.y='ensembl_gene_id', all.x=TRUE)
final_table <- final_table %>% distinct()
final_table <- final_table[!is.na(final_table$DiffLevene),]
  
write.table(final_table,'~/marenostrum/Projects/GTEx_v8/Methylation/Data/diff_variable_genes.txt', sep = '\t', quote = F, 
            row.names = F, col.names = T)

library(plyr) #Counting
#final_table$name <- paste0(final_table$seqnames,':',final_table$start,';',final_table$end,';',final_table$overlapping.genes)
counts <- ddply(final_table, .(cg, trait), nrow)
names(counts) <- c("CG", "trait", "number")
final_table$sign <- as.numeric(sign(final_table$DiffLevene))
logfc <- ddply(final_table, .(cg, trait), summarise, sign=length(unique(sign))*sign(unique(sign)[1]))
logfc <- logfc[!is.na(logfc$sign),]
logfc[abs(logfc$sign)>1,'sign'] <- abs(logfc$sign)[abs(logfc$sign)>1]
names(logfc) <- c("CG", "trait", "dir")

to_plot <- merge(counts, logfc, by=c("CG", "trait"))
library(ggplot2)
to_plot$dir <- as.factor(to_plot$dir)
to_plot$trait <- factor(to_plot$trait, c("Ancestry", "Age", "Sex"))
#to_plot <- to_plot[to_plot$number>1,]
ggplot(to_plot) + geom_jitter(aes(trait, number, col=dir), alpha=0.5) +
  theme_bw() + geom_violin(aes(x = trait, y = number, fill=trait),alpha=0.8) +
  scale_fill_manual(values = c('#C49122', '#4B8C61','#70A0DF'))




table(to_plot$trait, to_plot$dir)
table(to_plot$trait, to_plot$number)
table(to_plot$trait, to_plot$dir, to_plot$number)

View(to_plot[to_plot$number>=4,])

##### Final plot ####
head(to_plot)
to_plot$number <- as.factor(to_plot$number)

traits_cols <- c('#C49122','#4B8C61','#70A0DF')
names(traits_cols) <- c('Ancestry','Sex','Age')
strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols))
to_plot$trait <- factor(to_plot$trait, levels = c('Ancestry','Sex','Age'))

ggplot(to_plot, aes(y = number, fill=dir)) +
  geom_bar(position = 'fill', alpha=0.8) + 
  geom_text(aes(label=after_stat(count), x = after_stat(count+1.5)), stat='count', position='fill', size=3,hjust=.8, angle=45) +
  theme_bw() + ylab('Nº of Tissues') + xlab('Proportion shared CpG') +
  facet_wrap2(~ trait, strip = strip, nrow = 1) + theme_bw() + scale_x_continuous(breaks=seq(0, 1, 0.5))



#New tissue sharing plots:
head(to_plot)
to_plot$number <- as.factor(to_plot$number)

traits_cols <- c('#C49122','#4B8C61','#70A0DF')
names(traits_cols) <- c('EURv1','SEX2','AGE','BMI')
strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols))
to_plot$trait <- factor(to_plot$trait, levels = c('EURv1','SEX2','AGE','BMI'))

g <- ggplot(to_plot, aes(y = number, fill=dir)) +
  geom_bar(position = 'fill', alpha=0.8) + 
  geom_text(aes(label=after_stat(count), x = after_stat(count+1.5)), stat='count', position='fill', size=3,hjust=.8, angle=45) +
  theme_bw() + ylab('Nº of Tissues') + xlab('Proportion shared CpG') +
  facet_wrap2(~ trait, strip = strip, nrow = 1) + theme_bw() + scale_x_continuous(breaks=seq(0, 1, 0.5))
ggsave(paste0(project_path, 'Plots/tissue_sharing_DVGs_ancestry.pdf'), g, dpi = 300, width = 5, height = 2, units = 'in')

saveRDS(to_plot, '~/marenostrum/Projects/GTEx_v8/Methylation/Data/sharing_DVG.rds')

### proportion of shared per variable ####
#to_plot <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Sharing_DMP.rds')
total_DEGs <- sapply(c("Ancestry","Sex","Age"), function(trait) nrow(to_plot[to_plot$trait==trait,]))

pdf('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/DVG.barplot.total.pdf', height = 4, width = 4)
barplot(total_DEGs, col = traits_cols, names.arg=c("Ancestry","Sex","Age"), ylab = 'total unique nº of DMP')
dev.off()

plot.new()
legend("center",
       c("Ancestry","Sex","Age"),
       col = traits_cols,
       pch = 15,
       bty = 'n', ncol = 4)

d <- rbind.data.frame(sapply(c("Ancestry","Sex","Age"), function(trait) 
  100*(sum(to_plot$trait == trait & to_plot$number == '1')/total_DEGs[trait])
),
sapply(c("Ancestry","Sex","Age"), function(trait) 
  100*(sum(to_plot$trait == trait & to_plot$number %in% c('2','3'))/total_DEGs[trait])
),
sapply(c("Ancestry","Sex","Age"), function(trait) 
  100*(sum(to_plot$trait == trait & to_plot$number %in% c('4','5','6','7','8'))/total_DEGs[trait])
))
colnames(d) <- c("Ancestry","Sex","Age")
rownames(d) <- c("Tissue-specific", "Moderate sharing", "High sharing")

library(reshape)
library(ggh4x)
library(RColorBrewer)
df2 <- melt(d)
df2$type <- rep(c("Tissue-specific", "Moderate sharing", "High sharing"), 3)
df2$type <- factor(df2$type, levels = rev(c("Tissue-specific", "Low sharing", "Moderate sharing", "High sharing")), order = T) 
cols <- brewer.pal(4, "Greys")[c(2:4)]
names(cols) <- c("Tissue-specific","Moderate sharing", "High sharing")

traits_cols <- c('#C49122','#4B8C61','#70A0DF')
names(traits_cols) <- c("Ancestry","Sex","Age")

strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols))
to_plot$trait <- factor(to_plot$trait, levels = c("Ancestry","Sex","Age"))

p2 <- ggplot(df2, aes(x = 1,
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

counts <- counts[counts$trait=="Ancestry",]

colors_traits <- list('AGE'=c('#3D7CD0','#B4D6F6'),
                      'SEX2'=c('#3B734E','#89AA94'),
                      'EURv1'=c('#F0AE21','#F9DE8B'))
colors <- c('#F0AE21','#F9DE8B')
counts$sign[counts$sign==1] <- "Hyper"
counts$sign[counts$sign==-1] <- "Hypo"
g <- ggplot(counts) + geom_jitter(aes(sign, V1, color=sign, alpha=0.6)) +
  theme_bw() +
  geom_violin(aes(sign, V1, fill=sign),col = "black") +
  scale_color_manual(values = colors) + 
  scale_fill_manual(values = colors) + xlab("") + ylab("Number of tissues") +
  theme(legend.position = "none",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10))

ggsave(paste0(project_path, 'Plots/tissue_sharing_DVPs_ancestry.pdf'), g, dpi = 300, width = 2, height = 2, units = 'in')

##### overlap DVG with DVPs #####
data_meth <- list()
for(tissue in tissues){ 
  print(tissue)
  
  model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DVP_Age.rds"))
  model <- model[!is.na(model$Adj.P.Value),]
  data_meth[[tissue]][['Age']] <- gsub('\\..','',rownames(model[model$Adj.P.Value<0.05,]))
  model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DVP_Ancestry.rds"))
  model <- model[!is.na(model$Adj.P.Value),]
  data_meth[[tissue]][['Ancestry']]  <- gsub('\\..','',rownames(model[model$Adj.P.Value<0.05,]))
  if(file.exists(paste0(project_path, "Tissues/",tissue, "/DVP_Sex.rds"))){
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DVP_Sex.rds"))
    model <- model[!is.na(model$Adj.P.Value),]
    data_meth[[tissue]][['Sex']] <- gsub('\\..','',rownames(model[model$Adj.P.Value<0.05,]))
  }
  
}

data_expr <- list()
for(tissue in tissues){ 
  print(tissue)
  
  model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DVG_Age.rds"))
  model <- model[!is.na(model$Adj.P.Value),]
  data_expr[[tissue]][['Age']] <- rownames(model[model$Adj.P.Value<0.05,])
  model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DVG_Ancestry.rds"))
  model <- model[!is.na(model$Adj.P.Value),]
  data_expr[[tissue]][['Ancestry']]  <-rownames(model[model$Adj.P.Value<0.05,])
  if(file.exists(paste0(project_path, "Tissues/",tissue, "/DVG_Sex.rds"))){
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DVG_Sex.rds"))
    model <- model[!is.na(model$Adj.P.Value),]
    data_expr[[tissue]][['Sex']] <- rownames(model[model$Adj.P.Value<0.05,])
  }
  
}

### merge cpgs with gene annotation ###
annotation_genes <- read.delim('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Methylation_Epic_gene_promoter_enhancer_processed.txt')

data_meth_genes <- list()
for(tissue in tissues){ 
  print(tissue)
  
  data_meth_genes[[tissue]][['Age']] <- unique(annotation_genes$UCSC_RefGene_Name[annotation_genes$IlmnID %in% data_meth[[tissue]][['Age']]])
  data_meth_genes[[tissue]][['Ancestry']] <- unique(annotation_genes$UCSC_RefGene_Name[annotation_genes$IlmnID %in% data_meth[[tissue]][['Ancestry']]])
  if(file.exists(paste0(project_path, "Tissues/",tissue, "/DVP_Sex.rds"))){
    data_meth_genes[[tissue]][['Sex']] <- unique(annotation_genes$UCSC_RefGene_Name[annotation_genes$IlmnID %in% data_meth[[tissue]][['Sex']]])
  }
  
}

# library('biomaRt')
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# data_expr_gene_symbol <- list()

gene_annotation <- read.delim("~/marenostrum/MN4/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/data/public/gencode.v26.GRCh38.genes.bed", header=F)[,c(6,7)]
colnames(gene_annotation) <- c("gene", "symbol")

#These were genes duplicated, I changed their names to their correct one
gene_annotation$symbol[gene_annotation$gene=="ENSG00000253972.5"] <- "MAL2-AS1" #I found this a posteriori
gene_annotation$symbol[gene_annotation$gene=="ENSG00000283992.1"] <- "SLURP2" #Insted of LYNX1
gene_annotation$symbol[gene_annotation$gene=="ENSG00000235271.5"] <- "GC22P027299"
gene_annotation$symbol[gene_annotation$gene=="ENSG00000229694.6"] <- "C9orf73"
gene_annotation$symbol[gene_annotation$gene=="ENSG00000228741.2"] <- "GC13P024553" 

for(tissue in tissues){ 
  print(tissue)
  
  data_expr_gene_symbol[[tissue]][['Age']] <- gene_annotation$symbol[gene_annotation$gene %in% data_expr[[tissue]][['Age']]]
  if(length(data_expr[[tissue]][['Ancestry']])>0){
  data_expr_gene_symbol[[tissue]][['Ancestry']] <- gene_annotation$symbol[gene_annotation$gene %in% data_expr[[tissue]][['Ancestry']]]
  }
  if(file.exists(paste0(project_path, "Tissues/",tissue, "/DVP_Sex.rds"))){
    data_expr_gene_symbol[[tissue]][['Sex']] <- gene_annotation$symbol[gene_annotation$gene %in% data_expr[[tissue]][['Sex']]]
  }
  
}

overlap <- matrix(nrow = length(tissues), ncol=3, dimnames = list(tissues, names))
for(tissue in tissues){ 
  print(tissue)
  
  overlap[tissue, "Age"] <- sum(data_meth_genes[[tissue]][['Age']] %in% data_expr_gene_symbol[[tissue]][['Age']])
  overlap[tissue, "Ancestry"] <-  sum(data_meth_genes[[tissue]][['Ancestry']] %in% data_expr_gene_symbol[[tissue]][['Ancestry']])
  if(file.exists(paste0(project_path, "Tissues/",tissue, "/DVG_Sex.rds"))){
    overlap[tissue, "Sex"] <-  sum(data_meth_genes[[tissue]][['Sex']] %in% data_expr_gene_symbol[[tissue]][['Sex']])
  }
  
}

create_heatmap <- function(data, tissue_info, size=12){ #It takes as input the whole data frame, the whole information on tissues, the subset of tissues and  diseases to be plotted, and the font size
  #data, tissue_info, tissues_plot = rownames(data), diseases_plot = colnames(data), size=12
  without_NA <- replace(data, is.na(data), "")
  
  tissues_cols <- tissue_info[, 3]
  names(tissues_cols) <- tissue_info[, 1]
  tissues_cols <- tissues_cols[tissues]
  
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
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(my_pretty_num_function(without_NA[i, j]), x, y, gp = gpar(fontsize = size))}
          
  )
}
overlap <- overlap[,c("Ancestry", "Sex", "Age")]
rownames(overlap) <- tissue_info$name[match(rownames(overlap), tissue_info$tissue_ID)]
pdf(paste0(project_path, "Plots/Heatmap_DVGs_DMG.pdf"), height = 4.5, width = 6.5)
create_heatmap(overlap, tissue_info, size=12)
dev.off()


 #### enrichment overlap? 
