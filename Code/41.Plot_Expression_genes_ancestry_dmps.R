### compare expression levels of genes with DMPs in TSS ####
############# Plot expression of genes ########
#!/usr/bin/env Rscript
set.seed(1)
library(ggpubr)
library(vroom)
library(ggplot2)
library(ggtext)
library(ggpubr)

first_dir <- "~/marenostrum/"

project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Age", "Ancestry", "BMI", "Sex")

#tissue_info <- readRDS(paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Jose/00_Data/Tissue_info_whole.rds"))
tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Methylation/Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

# Gene annotation ----
# PCG and lincRNA genes with matched biotype annotation in gencode v38
# PAR genes excluded
#gene_annotation <- read.delim("/gpfs/projects/bsc83/Data/GTEx_v8/GeneAnnotation/gencode.v26.GRCh38.genes.biotype_matched_v38.bed")
gene_annotation <- read.delim("~/marenostrum/Projects/GTEx_v8/Methylation/Data/gencode.v26.gene_annotation.bed", header = F)



### plots like in cell genomics
get_box_stats <- function(y, upper_limit = max(y) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "n =", length(y)
    )
  ))
}

admixture_ancestry <- read.table('~/marenostrum_scratch/MN4/bsc83/bsc83535/GTEx/v8/genotype_data/admixture_inferred_ancestry.txt')
colnames(admixture_ancestry) <- c('Donor','AFRv1','EURv1','inferred_ancestry','AFRv2','EURv2')
admixture_ancestry$Ancestry <- 'ADX'
admixture_ancestry$Ancestry[admixture_ancestry$EURv1>=0.5] <- 'EUR'
admixture_ancestry$Ancestry[admixture_ancestry$EURv1<0.5] <- 'AFR'

metadata <- lapply(tissues, function(tissue) readRDS(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/",tissue,"/","metadata_expression.rds")))
names(metadata) <- tissues

metadata <- lapply(tissues, function(tissue) merge(metadata[[tissue]], admixture_ancestry[,c("Donor","Ancestry")], by='Donor')
)
names(metadata) <- tissues

get_gene_data <- function(tissue, df, gene_name){
  #df <- cbind.data.frame(t(tpm[gene_name,]))
  colnames(df) <- "TPM"
  df$Ancestry <- sapply(rownames(df), function(i) metadata[[tissue]][metadata[[tissue]]$Sample==i,"Ancestry.y"])
  df$Ancestry <- gsub("EUR", "EA", df$Ancestry)
  df$Ancestry <- gsub("AFR", "AA", df$Ancestry)
  df$Ancestry <- factor(df$Ancestry, levels = c("EA", "AA"), order = T)
  df$Sex <- sapply(rownames(df), function(i) metadata[[tissue]][metadata[[tissue]]$Sample==i,"Sex"])
  df$Sex <- gsub("1", "Male", df$Sex)
  df$Sex <- gsub("2", "Female", df$Sex)
  df$Sex <- factor(df$Sex, levels = c("Male", "Female"), order = T)
  df$Age_int <- sapply(rownames(df), function(i) metadata[[tissue]][metadata[[tissue]]$Sample==i,"Age"])
  df$Age <- sapply(df$Age_int, function(i) ifelse(i<45, "[20-45)", "[45-70]"))
  df$Age <- factor(df$Age, 
                   levels = c("[20-45)", "[45-70]"),
                   order = T)
  
  df$Tissue <- rep(tissue_info[tissue_info$tissue_ID==tissue, "tissue_abbrv"], nrow(df))
  df$Gene <- rep(gene_name, nrow(df))
  return(df)
}

traits_cols <- c('#C49122','#4B8C61','#70A0DF','#A76595')
names(traits_cols) <- c("Ancestry","Sex","Age","BMI")

get_gene_plot <- function(trait, tissue, gene_name){
  print(paste0(trait, " - ", tissue, " - ", gene_name))
  # Gene TPM --
  gene <- gene_annotation[gene_annotation$V6 == gene_name, "V4"]
  tpm <- cbind.data.frame(t(readRDS(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/", tissue, "/", "tpm.rds"))[gene,]))
  # tpm: matrix with gene TPM expression data for spleen samples used in this study. T
  
  # get gene data --
  data <- get_gene_data(tissue, tpm, gene_name)
  cols <- rep(traits_cols[trait], length(levels(data[,trait])))
  names(cols) <- as.character(levels(data[,trait]))
  
  # plot 1 --
  p1 <- ggplot(data = data,
               aes(x = eval(parse(text=trait)),
                   y = log10(TPM),
                   fill =  eval(parse(text=trait))),
  ) +
    geom_violin(col = "black") +
    geom_boxplot(col = "black",
                 fill = "white",
                 outlier.shape = NA,
                 notch = T,
                 width = 0.25) +
    geom_jitter(col = "black", 
                alpha = 0.1,
                size = 0.8) +
    xlab("") +
    ylab("log10(TPM)") +
    scale_fill_manual(values = cols) +
    stat_summary(fun.data = get_box_stats, geom = "text",
                 hjust = 0.5, vjust = 0.9, size = 3) +
    labs(title=paste0(gene_name, " (", tissue_info[tissue_info$tissue_ID==tissue, "tissue_abbrv"], ")")) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.title.y = ggtext::element_markdown(size=14),
          plot.title = element_text(hjust = 0.5,
                                    size = 15),
          strip.background = element_rect(fill="#B3B3B3"),
          legend.position = "none") 
  
  return(p1)
}

### plot FC 
## read DMPs

first_dir <- "~/marenostrum/"
annotation <- read.csv("~/marenostrum_scratch/MN4/bsc83/bsc83535/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")

project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("MuscleSkeletal","Testis",'KidneyCortex',"BreastMammaryTissue","Prostate","Lung", "ColonTransverse", "Ovary", "WholeBlood")
names <- c("Ancestry")
traits_to_use <- c('EURv1')

results_DML <- lapply(tissues, function(tis) 
  readRDS(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/",tis,"/DML_results_5_PEERs_continous.rds")))
names(results_DML) <- tissues

files <- list.files('~/marenostrum/Projects/GTEx_v8/Methylation/Data/EpiMap/', pattern='.bed.gz',full.names=T)

names_chrom <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "Breast", "MuscleSkeletal", "KidneyCortex", "Testis","PBMC")

chromhmm <- lapply(names_chrom, function(tis)
  read.delim(files[grep(tis, files)], sep='\t', header=F))
names(chromhmm) <- c("Lung","ColonTransverse",'Ovary',"Prostate","BreastMammaryTissue","MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")

### overlap chromhmm and cpgs tested #### 
ann_bed <- annotation[!is.na(annotation$MAPINFO) & !is.na(annotation$CHR),] %>%
  dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>% distinct()
ann_bed$chrom <- paste0('chr',ann_bed$chrom)
head(ann_bed)
ann_bed$start <- ann_bed$start-1

library(valr)
chromhmm_cpgs <- lapply(tissues, function(tis) {
  chrom_df <- chromhmm[[tis]][,c(1:4)]
  colnames(chrom_df) <- c('chrom','start','end','region')
  bed_intersect(ann_bed, chrom_df, suffix = c("_ann", "_chromhmm"))})
names(chromhmm_cpgs) <- tissues

meth_genes <- read.delim('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Methylation_Epic_gene_promoter_enhancer_processed.txt')

DEA_GTEx <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/DEA_GTEx.rds')

library(ggpubr)
for (tissue in tissues[tissues != 'KidneyCortex']) {
  print(tissue)
  res <- results_DML[[tissue]][[traits_to_use]]
  chrom_tissue <- as.data.frame(chromhmm_cpgs[[tissue]])  
  chrom_tissue$region_chromhmm_new <- chrom_tissue$region_chromhmm
  chrom_tissue <- chrom_tissue[chrom_tissue$.overlap ==1,]
  # 
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TssFlnkD", "TssFlnk", "TssFlnkU","TssA")] <- "TSS"
  tss_dmps <- rownames(res[rownames(res) %in% chrom_tissue$name_ann[chrom_tissue$region_chromhmm_new == 'TSS'] & res$adj.P.Val<0.05 & res$logFC<0,])
  other_tss <- rownames(res[rownames(res) %in% chrom_tissue$name_ann[chrom_tissue$region_chromhmm_new == 'TSS'] & res$adj.P.Val>0.05,])
  
  genes_tss_dmps <- unique(meth_genes$UCSC_RefGene_Name[meth_genes$IlmnID %in% tss_dmps])
  genes_tss_other <- unique(meth_genes$UCSC_RefGene_Name[meth_genes$IlmnID %in% other_tss])
  
  dea_tissue <- DEA_GTEx[[tissue]][['Ancestry']]
  dea_tissue <- dea_tissue[dea_tissue$gene_name %in% c(genes_tss_other, genes_tss_dmps),]
  
  dea_tissue$type <- 'Not Hyper AFR'
  dea_tissue$type[dea_tissue$gene_name %in% genes_tss_dmps] <- 'AFR hyper'
  pdf(paste0('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/Fc_expression_TSS_hyper_',tissue,'.pdf'), width = 4, height = 5)
  print(ggboxplot(dea_tissue, x = "type", y = "logFC",
                 color = "type", palette = "jco",
                 add = "jitter") + stat_compare_means() + geom_hline(yintercept = 0))
  dev.off()
  
  
}

### calculate IQR 

get_gene_data <- function(tissue, df){
  #df <- cbind.data.frame(t(tpm[gene_name,]))
  df$Ancestry <- sapply(rownames(df), function(i) metadata[[tissue]][metadata[[tissue]]$Sample==i,"Ancestry.y"])
  df$Ancestry <- gsub("EUR", "EA", df$Ancestry)
  df$Ancestry <- gsub("AFR", "AA", df$Ancestry)
  df$Ancestry <- factor(df$Ancestry, levels = c("EA", "AA"), order = T)
  df$Sex <- sapply(rownames(df), function(i) metadata[[tissue]][metadata[[tissue]]$Sample==i,"Sex"])
  df$Sex <- gsub("1", "Male", df$Sex)
  df$Sex <- gsub("2", "Female", df$Sex)
  df$Sex <- factor(df$Sex, levels = c("Male", "Female"), order = T)
  df$Age_int <- sapply(rownames(df), function(i) metadata[[tissue]][metadata[[tissue]]$Sample==i,"Age"])
  df$Age <- sapply(df$Age_int, function(i) ifelse(i<45, "[20-45)", "[45-70]"))
  df$Age <- factor(df$Age, 
                   levels = c("[20-45)", "[45-70]"),
                   order = T)
  
  df$Tissue <- rep(tissue_info[tissue_info$tissue_ID==tissue, "tissue_abbrv"], nrow(df))
  return(df)
}

rownames(gene_annotation) <- gene_annotation$V4

### correlations 

get_corr <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "SEX2"){
    NA
  }else{
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/",trait,'_Correlations_probes_genes_DEG_DMP.pnominal.rds'))
    #rownames(model[[trait]][model[[trait]]$adj.P.Val<0.05,])
    model
  }
}
DMPs_Res <- lapply(c('EURv1','SEX2','AGE','BMI'), function(trait) lapply(tissues, function(tissue) get_corr(tissue, trait)))
names(DMPs_Res) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(DMPs_Res[[trait]]) <- tissues}


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

all_cor <- unlist(signif, recursive = FALSE)
all_cor_df <- do.call("rbind", all_cor)
all_cor_df[,c('Trait','Tissue','Gene')] <- stringr::str_split_fixed(rownames(all_cor_df), '\\.', 3)
all_cor_df$label <- paste0(all_cor_df$gene,':',all_cor_df$probe, ':', all_cor_df$Trait, ':',all_cor_df$Tissue)

for (tissue in tissues[tissues != 'KidneyCortex']) {
  print(tissue)
  res <- results_DML[[tissue]][[traits_to_use]]
  chrom_tissue <- as.data.frame(chromhmm_cpgs[[tissue]])  
  chrom_tissue$region_chromhmm_new <- chrom_tissue$region_chromhmm
  chrom_tissue <- chrom_tissue[chrom_tissue$.overlap ==1,]
  # 
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TssFlnkD", "TssFlnk", "TssFlnkU","TssA")] <- "TSS"
  tss_dmps <- rownames(res[rownames(res) %in% chrom_tissue$name_ann[chrom_tissue$region_chromhmm_new == 'TSS'] & res$adj.P.Val<0.05 & res$logFC<0,])
  other_tss <- rownames(res[rownames(res) %in% chrom_tissue$name_ann[chrom_tissue$region_chromhmm_new == 'TSS'] & res$adj.P.Val>0.05,])
  
  genes_cor <- all_cor_df$gene[all_cor_df$probe %in% tss_dmps & all_cor_df$Tissue == tissue & all_cor_df$Trait=='Ancestry']
  genes_tss_dmps <- unique(meth_genes$UCSC_RefGene_Name[meth_genes$IlmnID %in% tss_dmps])
  
  if (length(genes_cor)>0) {
     gene <- gene_annotation[gene_annotation$V6 %in% genes_cor, "V4"]
  tpm <- cbind.data.frame(t(readRDS(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/", tissue, "/", "tpm.rds"))[gene,]))
  colnames(tpm) <- gene_annotation[colnames(tpm), "V6"]
  
  # get gene data --
  data <- get_gene_data(tissue, tpm)
  
  iqr <- lapply(genes_cor, function(gene) {
    iqr_afr <- IQR(data[data$Ancestry=='AA',gene], na.rm = T)
    iqr_eur <- IQR(data[data$Ancestry=='EA',gene], na.rm = T)
    return(list(AA=iqr_afr, EA=iqr_eur))
  })
  names(iqr) <- genes_cor
  iqr_df <- do.call('rbind.data.frame',iqr)
  library(data.table)
  iqr_df$gene <- rownames(iqr_df)
  long <- reshape2::melt(iqr_df, variable.name = "Ancestry", id.vars="gene")
  long$logIQR <- log10(long$value)
  pdf(paste0('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/IQR_expression_TSS_hyper_',tissue,'.pdf'), width = 4, height = 5)
  # print(ggplot(long,aes(Ancestry, logIQR)) + geom_violin(scale="width", fill="lightgray") + xlab("") + ylab("logIQR") +
  #         geom_boxplot(col = "black",
  #                      fill = "white",
  #                      outlier.shape = NA,
  #                      notch = T,
  #                      width = 0.25) +
  #         geom_jitter(col = "black", 
  #                     alpha = 0.1,
  #                     size = 0.8) +
  #         stat_compare_means(label = "p.format", paired = TRUE) +
  #         theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
  #                            panel.grid.major = element_blank(),
  #                            panel.grid.minor = element_blank(),))
  print(ggpaired(long, x = "Ancestry", y = "logIQR",
                 color = "Ancestry", line.color = "gray", line.size = 0.4,
                 palette = "jco")+
          stat_compare_means(paired = TRUE))
  dev.off()
  
  } else {
    next
  }
  
 

  
}
## plot levels of expression
for (tissue in tissues[tissues != 'KidneyCortex']) {
  print(tissue)
  res <- results_DML[[tissue]][[traits_to_use]]
  chrom_tissue <- as.data.frame(chromhmm_cpgs[[tissue]])  
  chrom_tissue$region_chromhmm_new <- chrom_tissue$region_chromhmm
  chrom_tissue <- chrom_tissue[chrom_tissue$.overlap ==1,]
  # 
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TssFlnkD", "TssFlnk", "TssFlnkU","TssA")] <- "TSS"
  tss_dmps <- rownames(res[rownames(res) %in% chrom_tissue$name_ann[chrom_tissue$region_chromhmm_new == 'TSS'] & res$adj.P.Val<0.05 & res$logFC<0,])
  other_tss <- rownames(res[rownames(res) %in% chrom_tissue$name_ann[chrom_tissue$region_chromhmm_new == 'TSS'] & res$adj.P.Val>0.05,])
  
  genes_cor <- all_cor_df$gene[all_cor_df$probe %in% tss_dmps & all_cor_df$Tissue == tissue & all_cor_df$Trait=='Ancestry']
  genes_tss_dmps <- unique(meth_genes$UCSC_RefGene_Name[meth_genes$IlmnID %in% tss_dmps])
  
  if (length(genes_cor)>0) {
    gene <- gene_annotation[gene_annotation$V6 %in% genes_cor, "V4"]
    tpm <- cbind.data.frame(t(readRDS(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/", tissue, "/", "tpm.rds"))[gene,]))
    colnames(tpm) <- gene_annotation[colnames(tpm), "V6"]
    
    # get gene data --
    data <- get_gene_data(tissue, tpm)
    
    iqr <- lapply(genes_cor, function(gene) {
      iqr_afr <- mean(data[data$Ancestry=='AA',gene], na.rm = T)
      iqr_eur <- mean(data[data$Ancestry=='EA',gene], na.rm = T)
      return(list(AA=iqr_afr, EA=iqr_eur))
    })
    names(iqr) <- genes_cor
    iqr_df <- do.call('rbind.data.frame',iqr)
    library(data.table)
    iqr_df$gene <- rownames(iqr_df)
    long <- reshape2::melt(iqr_df, variable.name = "Ancestry", id.vars="gene")
    long$logIQR <- log10(long$value)
    pdf(paste0('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/mean_expr_expression_TSS_hyper_',tissue,'.pdf'), width = 4, height = 5)
    # print(ggplot(long,aes(Ancestry, logIQR)) + geom_violin(scale="width", fill="lightgray") + xlab("") + ylab("logIQR") +
    #         geom_boxplot(col = "black",
    #                      fill = "white",
    #                      outlier.shape = NA,
    #                      notch = T,
    #                      width = 0.25) +
    #         geom_jitter(col = "black", 
    #                     alpha = 0.1,
    #                     size = 0.8) +
    #         stat_compare_means(label = "p.format", paired = TRUE) +
    #         theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    #                            panel.grid.major = element_blank(),
    #                            panel.grid.minor = element_blank(),))
    print(ggpaired(long, x = "Ancestry", y = "logIQR",
                   color = "Ancestry", line.color = "gray", line.size = 0.4,
                   palette = "jco")+
            stat_compare_means(paired = TRUE))
    dev.off()
    
  } else {
    next
  }
  
  
  
  
}
### Plot specific examples

gene_examples <- list(c("Ancestry", "Lung", "PLA2G4C"),
                      c("Ancestry", "ColonTransverse", "PLA2G4C"),
                      c("Ancestry", "Ovary", "STX1B"))
gene_examples <- list(c("Ancestry", "BreastMammaryTissue", "DHX58"),
                      c("Ancestry", "Lung", "CAMK1D"),
                      c("Ancestry", "Ovary", "CAMK1D"),
                      c("Ancestry", "Lung", "HEATR3"))
gene_examples <- list(c("Sex", "ColonTransverse", "SPESP1"),
                      c("Sex", "Lung", "SPESP1"),
                      c("Sex", "Lung", "MIR4458HG"))

library(ggpubr)
plot_path <- '~/marenostrum/Projects/GTEx_v8/Methylation/Plots/'
for(i in 1:length(gene_examples)){
  p <- get_gene_plot(gene_examples[[i]][1], gene_examples[[i]][2], gene_examples[[i]][3])  
  ggexport(p,filename = paste0(plot_path, "Figure_2E.", gene_examples[[i]][1], "_", gene_examples[[i]][2], "_", gene_examples[[i]][3], ".pdf"),
           width = 3, height = 4)
}
