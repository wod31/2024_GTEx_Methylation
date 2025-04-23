#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez; Modified by Winona Oliveros 
# @E-mail: jose.ramirez1@bsc.es
# @Description: Methylation replication in previous studies
# @software version: R=4.2.2

#Set path 
# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

library(ggplot2)

print("Reading results")
results <- readRDS("Tissues/WholeBlood/DML_results_5_PEERs_continous.rds")

#Papers to compare:
annotation <- read.csv("~/marenostrum_scratch/MN4/bsc83/bsc83535/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")
annotation$probe <- paste0(annotation$CHR,'_',annotation$MAPINFO)
background <- annotation$probe

eurv1 <- results$EURv1
signif <- eurv1[eurv1$adj.P.Val<0.05,]
hypo <- annotation$probe[annotation$Name %in% rownames(signif[signif$logFC<0,])]
hyper <- annotation$probe[annotation$Name %in% rownames(signif[signif$logFC>0,])]

#Blood https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5267325/#SD3 
name <- "Blood: Aracena et al. 2024"
data <- "WGBS"
sample_size <- 7463164 # CpG #S + NS

library(readxl)
all_tested_cpgs <- read.delim('~/Downloads/methylation_popDE_results.txt', sep = ' ')[,1]
blood <- read_excel("Data/Barreriro_Suppl_41588_2024_1668_MOESM8_ESM.xlsx", sheet = 14)
blood <- read_excel("Data/Barreriro_Suppl_41588_2024_1668_MOESM8_ESM.xlsx", sheet = 10)

library(tidyr)
blood <- blood %>% separate(feature, c('chrom', 'start','end'))
blood$start <- as.integer(blood$start)
blood$end <- as.integer(blood$end)

ann_bed <- annotation[!is.na(annotation$MAPINFO) & !is.na(annotation$CHR),] %>%
  dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>% distinct()
ann_bed$chrom <- paste0('chr',ann_bed$chrom)
head(ann_bed)
 
library(valr)
chromhmm_cpgs <- bed_intersect(blood, ann_bed, suffix = c("_ann", "_chromhmm"))

sum(blood$qvalue_NI<0.1) # 39684 enhancers
sum(blood$pvalue_NI<0.05)
##for methylatiom
blood$probe <- blood$feature
blood <- blood[blood$probe%in% background,] #Only 4160 out of 76956 (%5.4) of the loci are in the EPIC

blood_up <- blood[blood$diff_NI<0,] # EUR
blood_down <- blood[blood$diff_NI>0,] # AFR
common <- c(blood_up$probe[blood_up$probe %in% hyper], blood_down$probe[blood_down$probe %in% hypo])
only_blood <- table(c(blood_up$probe, blood_down$probe) %in% common)[1]
only_ours <- table(c(hypo, hyper)[c(hypo, hyper) %in% all_tested_cpgs] %in% common)[1]
bg <- background[!background %in% common]
bg <- bg[!bg %in% only_blood]
bg <- bg[!bg %in% only_ours]
bg <- bg[bg %in% all_tested_cpgs]
m <- matrix(c(length(common), only_blood, only_ours, length(bg)), nrow=2)
# t <- fisher.test(m, alternative = "greater") #common, signif only in blood, signif only in lung, not signig
t <- fisher.test(m) #common, signif only in blood, signif only in lung, not signig
p_vals <- t$p.value
ors <- t$estimate
CI_down <- t$conf.int[1]
CI_up <- t$conf.int[2]
overlaps <- length(common)
totals <- length(common) + only_blood
names <- name
samples <- length(background[background %in% all_tested_cpgs])
dmps <- nrow(blood) #positions they reported that are available in the EPIC array (mention in the legend)

### for chip-seq - enhancer
# 113584 peaks, 23529 sig
blood_sig <- chromhmm_cpgs[chromhmm_cpgs$qvalue_NI_ann<0.05 & chromhmm_cpgs$lfsr_NI_ann<0.1,] # 69718
#blood_up <- blood_sig[blood_sig$logFC_NI_ann<0,] # EUR
#blood_down <- blood_sig[blood_sig$logFC_NI_ann>0,] # AFR

signif <- eurv1[eurv1$adj.P.Val<0.05,]
#hypo <- rownames(signif[signif$logFC<0,])
#hyper <- rownames(signif[signif$logFC>0,])

common <- blood_sig[blood_sig$name_chromhmm %in% rownames(signif),]

result <- phyper(q = nrow(common) - 1,
                 m = length(unique(blood_sig$name_chromhmm)),
                 n = length(unique(chromhmm_cpgs$name_chromhmm)) - length(unique(blood_sig$name_chromhmm)),
                 k = sum(rownames(signif) %in% chromhmm_cpgs$name_chromhmm),
                 lower.tail = FALSE)


### sex whole blood ###
eurv1 <- results$SEX2
signif <- eurv1[eurv1$adj.P.Val<0.05,]
hypo <- annotation$Name[annotation$Name %in% rownames(signif[signif$logFC<0,])]
hyper <- annotation$Name[annotation$Name %in% rownames(signif[signif$logFC>0,])]

#Blood https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5267325/#SD3 
name <- "Blood: GRANT"
data <- "WGBS"
sample_size <- 7463164 # CpG #S + NS

library(readxl)
#all_tested_cpgs <- read.delim('~/Downloads/methylation_popDE_results.txt', sep = ' ')[,1]
blood <- read.delim("Data/papers/Grant_DMPs_autosomal_13148_2022_1279_MOESM1_ESM.csv", sep=',')

#blood <- blood[blood$probe%in% background,] #Only 4160 out of 76956 (%5.4) of the loci are in the EPIC

blood_up <- blood[blood$logFC<0,] # EUR
blood_down <- blood[blood$logFC>0,] # AFR
common <- c(blood_up$Row.names[blood_up$Row.names %in% hyper], blood_down$Row.names[blood_down$Row.names %in% hypo])
only_blood <- table(c(blood_up$Row.names, blood_down$Row.names) %in% common)[1]
only_ours <- table(c(hypo, hyper) %in% common)[1]
background <- rownames(eurv1)
bg <- background[!background %in% common]
bg <- bg[!bg %in% only_blood]
bg <- bg[!bg %in% only_ours]
#bg <- bg[bg %in% all_tested_cpgs]
m <- matrix(c(length(common), only_blood, only_ours, length(bg)), nrow=2)
# t <- fisher.test(m, alternative = "greater") #common, signif only in blood, signif only in lung, not signig
t <- fisher.test(m) #common, signif only in blood, signif only in lung, not signig

### fisher blood results ####
annotation <- read.csv("~/marenostrum_scratch/MN4/bsc83/bsc83535/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")
tissues <- c("WholeBlood")
names <- c("Sex")
traits_to_use <- c('SEX2')

files <- list.files('Data/EpiMap/', pattern='.bed.gz',full.names=T)

names_chrom <- c("PBMC")
chromhmm <- lapply(names_chrom, function(tis)
  read.delim(files[grep(tis, files)], sep='\t', header=F))
names(chromhmm) <- tissues

### overlap chromhmm and cpgs tested ####
# ann_granges <- annotation[!is.na(annotation$MAPINFO),] %>%
#   dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>% distinct() %>%
#   makeGRangesFromDataFrame(keep.extra.columns=T)
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

my_fisher <- function(type, tissue){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  res <- eurv1
  chrom_tissue <- as.data.frame(chromhmm_cpgs[[tissue]])
  chrom_tissue$region_chromhmm_new <- chrom_tissue$region_chromhmm
  chrom_tissue <- chrom_tissue[chrom_tissue$.overlap ==1,]
  #
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TssFlnkD", "TssFlnk", "TssFlnkU","TssA")] <- "TSS"
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("EnhA2", "EnhA1","EnhWk","EnhG1", "EnhG2")] <- "Enh"
  
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("ReprPCWk","ReprPC")] <- "ReprPC"
  
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TxWk","Tx")] <- "Tx"
  
  #all_cpgs <- unique(c(unique(rownames(genes_sex)), unique(rownames(genes_age))))
  all_cpgs <- rownames(res)
  all_cpgs <- blood$Row.names
  
  type_df <- chrom_tissue[chrom_tissue$region_chromhmm_new == type & chrom_tissue$name_ann %in% all_cpgs,]
  other_type <- chrom_tissue[chrom_tissue$region_chromhmm_new != type & chrom_tissue$name_ann %in% all_cpgs,]
  type_diff <- nrow(type_df[type_df$name_ann %in% (blood$Row.names[blood$logFC<0]),])
  type_notdiff <- nrow(type_df) - type_diff
  other_type_diff <- nrow(other_type[other_type$name_ann %in% (blood$Row.names[blood$logFC<0]),])
  other_type_notdiff <- nrow(other_type) - other_type_diff
  
  ### test significance
  m <- matrix(c(type_diff, type_notdiff, other_type_diff, other_type_notdiff), 2,2, byrow = T)
  print(m)
  
  m[is.na(m)] <- 0
  #m <- m[c(type,paste0('No ',type)),]
  rownames(m) <- c(type, "Other")
  colnames(m) <- c("Hyper","Not Hyper")
  print(m)
  f <- fisher.test(m)
  print(f)
  return(list("f" = f, "m" = type_diff))
}
# Two-tailed Fisher test
#families <- as.vector(unique(shared_cpgs$region_chromhmm))
families <- c('Enh','EnhBiv','Het','Quies','ReprPC','TSS','TssBiv','Tx','ZNF/Rpts')
fisher_results <- lapply(families, function(region) my_fisher(region,'WholeBlood'))
names(fisher_results) <- families


saveRDS(fisher_results, 'Tissues/enrichment_chromhmm_hyperFemale_Blood_Grant_vs_DMP.rds')

read_data <- function(variables, data){ #Function to prepare data to plot and compute adjusted p value
  
  odds_ratio <- lapply(variables, function(type) data[[type]][['f']]$estimate)
  adj.P.Val <- p.adjust(sapply(variables, function(type) data[[type]][['f']]$p.value), method = "BH")
  CI_down <- lapply(variables, function(type) data[[type]][['f']]$conf.int[1])
  CI_up <- lapply(variables, function(type) data[[type]][['f']]$conf.int[2])
  sample_size <- lapply(variables, function(type) data[[type]][['m']])
  
  
  names(odds_ratio) <- variables
  names(adj.P.Val) <- variables
  names(CI_down) <- variables
  names(CI_up) <- variables
  names(sample_size) <- variables
  
  odds_ratio_df <- as.data.frame(unlist(odds_ratio))
  odds_ratio_df$label <- variables
  odds_ratio_df$type <- deparse(substitute(data)) #Either hypo or hyper
  colnames(odds_ratio_df) <- c('oddsRatio', 'region','type')
  
  adj.P.Val_df <- as.data.frame(unlist(adj.P.Val))
  adj.P.Val_df$label <- variables
  adj.P.Val_df$type <- deparse(substitute(data))
  colnames(adj.P.Val_df) <- c('adjPvalue','region','type')
  
  CI_down_df <- as.data.frame(unlist(CI_down))
  CI_down_df$label <- variables
  CI_down_df$type <- deparse(substitute(data))
  colnames(CI_down_df) <- c('CI_down','region','type')
  
  CI_up_df <- as.data.frame(unlist(CI_up))
  CI_up_df$label <- variables
  CI_up_df$type <- deparse(substitute(data))
  colnames(CI_up_df) <- c('CI_up','region','type')
  
  sample_size_df <- as.data.frame(unlist(sample_size))
  sample_size_df$label <- variables
  sample_size_df$type <- deparse(substitute(data))
  colnames(sample_size_df) <- c('sample_size','region','type')
  
  all <- Reduce(function(x, y) merge(x, y, all=TRUE), list(odds_ratio_df, adj.P.Val_df, CI_down_df, CI_up_df, sample_size_df))
  head(all)
  all$sig <- 'not Sig'
  all$sig[all$adjPvalue<0.05] <- 'Sig'
  all <- all[,c("region","oddsRatio","adjPvalue","CI_down","CI_up","sig","type", "sample_size")]
  return(all)
}

##### plot shared positions
female <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hyperFemale_Blood_Grant_vs_DMP.rds')
male <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hyperMale_Blood_Grant_vs_DMP.rds')

library(ggpubr)

#'Sex'=c('#3B734E','#89AA94')

hypo_d <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), male)
hyper_d <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), female)
hyper_hypo <- rbind(hypo_d, hyper_d)
hyper_hypo$sig[hyper_hypo$sig =="Sig"] <- "FDR < 0.05"
hyper_hypo$sig[hyper_hypo$sig =="not Sig"] <- "FDR >= 0.05"
hyper_hypo$sig <- factor(hyper_hypo$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
hyper_hypo$region <- factor(hyper_hypo$region, levels=rev(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts")))
hyper_hypo$type[hyper_hypo$type =="female"] <- "F"
hyper_hypo$type[hyper_hypo$type =="male"] <- "M"
g <- ggplot(hyper_hypo, aes(x=log2(oddsRatio), y=region, colour=type, alpha=sig)) +
  geom_errorbar(aes(xmin=log2(CI_down), xmax=log2(CI_up)), width=.3) +
  geom_vline(xintercept = 0) +
  #xlim(0,20) + #Only for Lung to show the 0
  geom_point(size=3) + ylab('') + theme_bw() +
  scale_colour_manual(values=(c('#3B734E','#89AA94'))) +
  xlab("log2(Odds ratio)") +
  scale_alpha_discrete(range = c(0.4, 1), drop = FALSE) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(colour="black", size=13),
        axis.text.y = element_text(colour="black", size=14),
        legend.text = element_text(colour="black", size=13),
        axis.title.x = element_text(size=16),
        legend.spacing.y = unit(-0.05, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth=1)) +
  scale_y_discrete(breaks=c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"),
                   labels=c("Enhancer","Enhancer Bivalent","Heterochromatin","Quiescent","Repressed Polycomb","TSS","TSS Bivalent","Transcription","ZNF & Repeats"))# + xlim(0, 3)
# pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/genomic_location_shared",'_',trait,".pdf"), w = 6, h = 3.5)
# print(g)
# dev.off()

#Plot sample sizes:

g2 <- ggplot(hyper_hypo) + geom_col(aes(sample_size, region, fill=type), width = 0.6) +
  theme_classic() + xlab("Number of DMPs") + ylab("") +
  scale_fill_manual(values=(c('#3B734E','#89AA94'))) +
  theme(legend.position = "none",
        axis.text.x = element_text(colour="black", size=13),
        axis.text.y=element_blank(),
        axis.title.x = element_text(size=16)) +
  scale_x_continuous(n.breaks=3)+
  scale_y_discrete(breaks=c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"),
                   labels=c("Enhancer","Enhancer Bivalent","Heterochromatin","Quiescent","Repressed Polycomb","TSS","TSS Bivalent","Transcription","ZNF & Repeats")) #+
#scale_x_continuous(breaks=c(0, 20000, 40000)) #Only for lung
# pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/genomic_location_shared",'_',trait,"_sample_size.pdf"), w = 4, h = 3.5)
# print(g2)
# dev.off()

p <- ggarrange(g, g2, labels = c("A", "B"),
               common.legend = TRUE, legend = "right", widths = c(0.8,0.3))
pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/enrichment_Sex_grant_vs_DMP.pdf"), w = 8, h = 4)
print(p)
dev.off()


### enrichment miss methyl #####
library(missMethyl)
res_f <- gometh(blood$Row.names[blood$logFC<0], all.cpg=rownames(eurv1),
              collection="GO", array.type="EPIC")
