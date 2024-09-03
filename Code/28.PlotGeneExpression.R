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
tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Jose/00_Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

# Gene annotation ----
# PCG and lincRNA genes with matched biotype annotation in gencode v38
# PAR genes excluded
#gene_annotation <- read.delim("/gpfs/projects/bsc83/Data/GTEx_v8/GeneAnnotation/gencode.v26.GRCh38.genes.biotype_matched_v38.bed")
gene_annotation <- read.delim("~/marenostrum/Data/gene_annotation/gencode/release_26/gencode.v26.gene_annotation.bed", header = F)



### plots like in cell genomics
get_box_stats <- function(y, upper_limit = max(y) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "n =", length(y)
    )
  ))
}

admixture_ancestry <- read.table('~/marenostrum_scratch/GTEx/v8/genotype_data/admixture_inferred_ancestry.txt')
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
