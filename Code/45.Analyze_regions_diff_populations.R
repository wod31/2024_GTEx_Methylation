#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Analyze chromatin regions that are specific of each population
# @software version: R=4.2.2

###### Analysis pop specific regions #####

#### read catalogue info####
catalogue <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/catalogue_population_specific_wholeblood.rds')

### perform enrichment ####
library(missMethyl)
res <- gometh(catalogue$name[(catalogue$European == 'Enh' | catalogue$`African-American` =='Enh') & catalogue$region!='Same' & catalogue$DMP=='DMP'],
              all.cpg=catalogue$name[(catalogue$European == 'Enh' | catalogue$`African-American` =='Enh') & catalogue$region!='Same'], collection="GO", array.type="EPIC", sig.genes = TRUE)

res <- gometh(catalogue$name[(catalogue$European == 'Enh' | catalogue$`African-American` =='Enh') & catalogue$region!='Same' & catalogue$DMP=='DMP'],
              all.cpg=catalogue$name, collection="GO", array.type="EPIC", sig.genes = TRUE)

res_all <- gometh(catalogue$name[catalogue$region!='Same' & catalogue$DMP=='DMP'],
              all.cpg=catalogue$name, collection="GO", array.type="EPIC", sig.genes = TRUE)


### TF enrichment ###

suppressMessages(library(dplyr))
library(parallel)
library(ggplot2)
suppressMessages(library(data.table))
library(tidyr)

#Set path 
setwd(system("pwd", intern = T)) #If in linux
setwd('~/marenostrum/Projects/GTEx_v8/Methylation/')

n_cores <- 2

first_path <- "~/marenostrum/"
data_path <- paste0(first_path, "Projects/GTEx_v8/Methylation/TFBS/")

fisher_function <- function(target, dmps){
  target_positions <- unique(annotation$name_ann[annotation$TF==target])
  non_target_positions <- rownames(results)[!rownames(results) %in% target_positions]
  non_dmp <- rownames(results)[!rownames(results) %in% dmps]# & results$adj.P.Val<0.05]
  m <- matrix(c(sum(dmps %in% target_positions),
                sum(dmps %in% non_target_positions),
                sum(non_dmp %in% target_positions),
                sum(non_dmp %in% non_target_positions)), nrow=2)
  t <- fisher.test(m)
  return(list(t$p.value, t$estimate, t$conf.int[1], t$conf.int[2]))
}

Sys.time()

tissues <- c('WholeBlood')
files <- list.files('~/marenostrum/Projects/GTEx_v8/Methylation/Data/EpiMap/', pattern='.bed.gz',full.names=T)

names_chrom <- c("PBMC")

chromhmm <- lapply(names_chrom, function(tis)
  read.delim(files[grep(tis, files)], sep='\t', header=F))
names(chromhmm) <- tissues

annotation <- read.csv("~/marenostrum_scratch/MN4/bsc83/bsc83535/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")
ann_bed <- annotation[!is.na(annotation$MAPINFO) & !is.na(annotation$CHR),] %>%
  dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>% distinct()
ann_bed$chrom <- paste0('chr',ann_bed$chrom)
head(ann_bed)
ann_bed$start <- ann_bed$start-1

names_table <- cbind(tissues, "short"=c("Bld"))

catalogue <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/catalogue_population_specific_wholeblood.rds')

### run enrichment
  #  print(chrom)
  short <- names_table[names_table[,1]==tissues[1],2]
  
  
  print("Reading data")
  results_DML <- readRDS(paste0("Tissues/", "WholeBlood" , "/DML_results_5_PEERs_continous.rds"))
  
  annotation <- read.csv(paste0(data_path, "chip_seq_", short,".csv"))[,-1]
  colnames(annotation)[5] <- c("TF")
  print(paste("We will be testing", length(unique(annotation$TF)), "TFs"))
  
  to_save <- data.frame(TF="1", tissue="1",state="1", variable="1", methylation="1", direction="1", p_value="1", odds_ratio="1", CI_low="1", CI_high="1")
  to_save$p_value <- as.numeric(to_save$p_value) #These 4 lines are needed to run in previous versions of R
  to_save$odds_ratio <- as.numeric(to_save$odds_ratio)
  to_save$CI_low <- as.numeric(to_save$CI_low)
  to_save$CI_high <- as.numeric(to_save$CI_high)

  variable <- 'EURv1'
    print(variable)
    results <- results_DML[[variable]]
    
    diff_dmps <- catalogue$name[catalogue$region!='Same' & catalogue$DMP=='DMP']

    results_hyper <- rownames(results[rownames(results) %in% diff_dmps & results$logFC>0,])
    results_hypo <- rownames(results[rownames(results) %in% diff_dmps & results$logFC<0,])
    
    if(length(results_hyper)>0){
      print("Computing fisher's for hypermethylated positions")
      Sys.time()
      hyper <- mclapply(unique(annotation$TF), function(target) fisher_function(target,results_hyper),
                        mc.cores=n_cores) 
      names(hyper) <- unique(annotation$TF)
      hyper <- do.call(rbind, hyper)
      hyper <- as.data.frame(hyper)
      hyper[,1] <- unlist(hyper[,1])
      hyper[,2] <- unlist(hyper[,2])
      hyper[,3] <- unlist(hyper[,3])
      hyper[,4] <- unlist(hyper[,4])
      hyper[,5] <- NULL
      colnames(hyper) <- c("p_value", "odds_ratio", "CI_low", "CI_high")
      hyper$TF <- rownames(hyper)
      #saveRDS(hyper, paste0("Tissues/", tissue , "/TFBS_hyper_", variable,".rds"))
      Sys.time()
      
      # hyper_enriched[tissue, variable] <- sum(hyper$adj_p_val<0.05 & hyper$odds_ratio>1)
      direction <- ifelse(hyper$odds_ratio>1, "Enriched", "Depleted")
      to_save <- rbind(to_save, cbind(tissue=tissue, state='shared', variable=variable, methylation="Hyper", direction=direction, hyper))
    }
    
    if(length(results_hypo)>0){
      print("Computing fisher's for hypomethylated positions")
      Sys.time()
      hypo <- mclapply(unique(annotation$TF), function(target) fisher_function(target, results_hypo),
                       mc.cores=n_cores)
      names(hypo) <- unique(annotation$TF)
      hypo <- do.call(rbind, hypo)
      hypo <- as.data.frame(hypo)
      hypo[,1] <- unlist(hypo[,1])
      hypo[,2] <- unlist(hypo[,2])
      hypo[,3] <- unlist(hypo[,3])
      hypo[,4] <- unlist(hypo[,4])
      hypo[,5] <- NULL
      colnames(hypo) <- c("p_value", "odds_ratio", "CI_low", "CI_high")
      hypo$TF <- rownames(hypo)
      
      # hypo$adj_p_val <- p.adjust(hypo$p_value, "BH")
      # saveRDS(hypo, paste0("Tissues/", tissue , "/TFBS_hypo_", variable,".rds"))
      Sys.time()
      
      # hypo_enriched[tissue, variable] <- sum(hypo$adj_p_val<0.05 & hypo$odds_ratio>1)
      direction <- ifelse(hypo$odds_ratio>1, "Enriched", "Depleted")
      to_save <- rbind(to_save, cbind(tissue=tissue, state='shared', variable=variable, methylation="Hypo", direction=direction, hypo))
    }
    
    all <- rownames(results[rownames(results) %in% diff_dmps,])
    if(length(all)>0){
      print("Computing fisher's for hypomethylated positions")
      Sys.time()
      hypo <- mclapply(unique(annotation$TF), function(target) fisher_function(target, all),
                       mc.cores=n_cores)
      names(hypo) <- unique(annotation$TF)
      hypo <- do.call(rbind, hypo)
      hypo <- as.data.frame(hypo)
      hypo[,1] <- unlist(hypo[,1])
      hypo[,2] <- unlist(hypo[,2])
      hypo[,3] <- unlist(hypo[,3])
      hypo[,4] <- unlist(hypo[,4])
      hypo[,5] <- NULL
      colnames(hypo) <- c("p_value", "odds_ratio", "CI_low", "CI_high")
      hypo$TF <- rownames(hypo)
      
      # hypo$adj_p_val <- p.adjust(hypo$p_value, "BH")
      # saveRDS(hypo, paste0("Tissues/", tissue , "/TFBS_hypo_", variable,".rds"))
      Sys.time()
      
      # hypo_enriched[tissue, variable] <- sum(hypo$adj_p_val<0.05 & hypo$odds_ratio>1)
      direction <- ifelse(hypo$odds_ratio>1, "Enriched", "Depleted")
      to_save <- rbind(to_save, cbind(tissue=tissue, state='shared', variable=variable, methylation="Hypo", direction=direction, hypo))
    }
  
  to_save <- to_save[-1,]
  to_save$adj_p_val <- p.adjust(to_save$p_value, "BH")
  to_save$direction <- factor(to_save$direction, levels=c("Enriched", "Depleted"))
  to_save$adj_p_val <- as.numeric(to_save$adj_p_val)
  to_save <- to_save[order(to_save$variable, to_save$methylation, to_save$direction, to_save$adj_p_val), ]
  write.csv(to_save, paste0(data_path, "tfbs_results_", tissue, "_tested_",'pop_specific',"_dmpsbg.csv"))
  print("done")

  
  ### look at how enhancer-specific DMPs affect expression #####
  gene_annotation <- read.delim("~/marenostrum/Projects/GTEx_v8/Methylation/Data/gencode.v26.gene_annotation.bed", header = F)
  
  tissue_info <- readRDS('Data/Tissue_info_whole.rds')
  
  tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
  tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

  
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
  
  meth_genes <- read.delim('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Methylation_Epic_gene_promoter_enhancer_processed.txt')
  
  diff_dmps <- catalogue$name[(catalogue$European == 'Enh' | catalogue$`African-American` =='Enh') & catalogue$region!='Same' & catalogue$DMP=='DMP']
  
  results_hyper <- rownames(results[rownames(results) %in% diff_dmps & results$logFC>0,])
  results_hypo <- rownames(results[rownames(results) %in% diff_dmps & results$logFC<0,])
  
  genes_hyper <- unique(meth_genes$UCSC_RefGene_Name[meth_genes$IlmnID %in% results_hyper])
  genes_hypo <- unique(meth_genes$UCSC_RefGene_Name[meth_genes$IlmnID %in% results_hypo])
  
  #cpgs <- catalogue$name[(catalogue$European == 'Enh' | catalogue$`African-American` =='Enh') & catalogue$region!='Same' & catalogue$DMP=='DMP']

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
  library(ggpubr)
  for (tissue in tissues) {
    print(tissue)
    
    if (length(genes_hyper)>0) {
      gene <- gene_annotation[gene_annotation$V6 %in% genes_hyper, "V4"]
      tpm <- cbind.data.frame(t(readRDS(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/", tissue, "/", "tpm.rds"))[gene,]))
      colnames(tpm) <- gene_annotation[colnames(tpm), "V6"]
      
      # get gene data --
      data <- get_gene_data(tissue, tpm)
      
      iqr <- lapply(genes_hyper, function(gene) {
        iqr_afr <- mean(data[data$Ancestry=='AA',gene], na.rm = T)
        iqr_eur <- mean(data[data$Ancestry=='EA',gene], na.rm = T)
        return(list(AA=iqr_afr, EA=iqr_eur))
      })
      names(iqr) <- genes_hyper
      iqr_df <- do.call('rbind.data.frame',iqr)
      library(data.table)
      iqr_df$gene <- rownames(iqr_df)
      long <- reshape2::melt(iqr_df, variable.name = "Ancestry", id.vars="gene")
      long$logIQR <- log10(long$value)
      pdf(paste0('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/mean_expr_expression_enhancer_dmp_ea_',tissue,'.pdf'), width = 4, height = 5)
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
    
    if (length(genes_hypo)>0) {
      gene <- gene_annotation[gene_annotation$V6 %in% genes_hypo, "V4"]
      tpm <- cbind.data.frame(t(readRDS(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/", tissue, "/", "tpm.rds"))[gene,]))
      colnames(tpm) <- gene_annotation[colnames(tpm), "V6"]
      
      # get gene data --
      data <- get_gene_data(tissue, tpm)
      
      iqr <- lapply(genes_hypo, function(gene) {
        iqr_afr <- mean(data[data$Ancestry=='AA',gene], na.rm = T)
        iqr_eur <- mean(data[data$Ancestry=='EA',gene], na.rm = T)
        return(list(AA=iqr_afr, EA=iqr_eur))
      })
      names(iqr) <- genes_hypo
      iqr_df <- do.call('rbind.data.frame',iqr)
      library(data.table)
      iqr_df$gene <- rownames(iqr_df)
      long <- reshape2::melt(iqr_df, variable.name = "Ancestry", id.vars="gene")
      long$logIQR <- log10(long$value)
      pdf(paste0('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/mean_expr_expression_enhancer_dmp_aa_',tissue,'.pdf'), width = 4, height = 5)
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
