#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Plot Results cis-driven analysis
# @software version: R=4.2.2

set.seed(1)

# Libraries ----
library(RColorBrewer)
library(ComplexHeatmap)
# GO enrichment --
#library(clusterProfiler)
#library(WebGestaltR)
library(org.Hs.eg.db)
library(ggplot2)
library(cowplot)
# Box plots
library(rcompanion)
library(gtools)

# ---- Data ---- ####

# Individual Traits ----
traits_cols <- c("Age" = "#56B4E9",
                 "Ancestry" = "#E69F00",
                 "Sex" =  "#009E73",
                 "BMI" = "#CC79A7")
traits <- names(traits_cols)

# Tissues ---
first_dir <- "/gpfs/projects/bsc83/"
first_dir <- "~/marenostrum/"

project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Age", "Ancestry", "BMI", "Sex")
#data <- matrix(nrow = length(tissues), ncol=4, dimnames = list(tissues, names))

tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Methylation/Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

# Tissue metadata ----
metadata <- lapply(tissues, function(tissue) readRDS(paste0(project_path,"Tissues/", tissue, "/metadata.rds")))
names(metadata) <- tissues

# ---- Analysis ---- ####
# 1.1 mGene data ----
# inpath_mqtls <- "~/marenostrum_scratch/GTEx/v9/mQTLs/"
# mCpGs <- lapply(tissues, function(tis) {
#   mgene_data <- as.data.frame(data.table::fread(paste0(inpath_mqtls,tis,".mQTLs.conditional.txt.gz")))[mgene_data$V7<0.05 & abs(mgene_data$V3)<250000,]
#   #mgene_data <- mgene_data[mgene_data$V7<0.05 & abs(mgene_data$V3)<250000,]
#   unique(gsub(':.*','',mgene_data$V1[mgene_data$V7<0.05 & abs(mgene_data$V3)<250000]))})
# names(mCpGs) <- tissues


# 1.2 Differential methylation results ----
results_DML <- lapply(tissues, function(tis) 
  readRDS(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/",tis,"/DML_results_5_PEERs_continous.rds")))
names(results_DML) <- tissues

# DMP ---
results_DML_all <- lapply(tissues, function(tis) 
  results_DML[[tis]][['EURv1']])
names(results_DML_all) <- tissues

## sig results 
DMP <- lapply(tissues, function(tis) 
  rownames(results_DML_all[[tis]][results_DML_all[[tis]]$adj.P.Val<0.05,]))
names(DMP) <- tissues

# 1.3 cis-driven classification ----
inpath <- "~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/"
for(tissue in tissues){
  if(!file.exists(paste0(inpath, tissue, "/", tissue, ".Ancestry_DMP.Classification_summary.qval005.rds"))){print(tissue)}
}
d <- lapply(tissues, function(tissue)
  readRDS(paste0(inpath, tissue, "/", tissue, ".Ancestry_DMP.Classification_summary.qval005.rds")))
names(d) <- tissues

# Create dataframe --
data <- do.call(rbind.data.frame,d)
rownames(data) <- tissues
colnames(data) <- as.character(names(d[[tissues[1]]]))

# DMP classified as either cis-driven or not cis-driven --
results <-  lapply(tissues, function(tissue) readRDS(paste0(inpath, tissue, "/", "Ancestry_DMP.Classified.qval005.rds")))
names(results) <- tissues

# Expand dataframe --
data <- cbind.data.frame(data,
                         sapply(tissues, function(tissue) sum(results[[tissue]]$Class=="Cis-driven", na.rm = T)),
                         sapply(tissues, function(tissue) sum(results[[tissue]]$Class=="Not_cis-driven", na.rm = T)),
                         sapply(tissues, function(tissue) sum(is.na(results[[tissue]]$Class), na.rm = T))
)
colnames(data)[10] <- "Cis-driven"
colnames(data)[11] <- "Not_cis-driven"
colnames(data)[12] <- "Too_Many_mQTL"
exprs.data <- data

# Ancestry DEG data.frame  ----
dd <- data[,c("Ancestry:DMP",
              "Ancestry:DMP:eGpG",
              #"Ancestry:DEG:eGenes-ieQTL",
              #"Ancestry:DEG:eGenes-ieQTL:MAF01",
              #"Ancestry:DEG:Modelled",
              "Cis-driven",
              "Not_cis-driven", 
              "Too_Many_mQTL",
              "Ancestry:DMP:NotModelled")]

dd$`Ancestry:DMP:not_eCpG` <- dd$`Ancestry:DMP` - dd$`Ancestry:DMP:eGpG`
dd <- dd[,c(1,7,2,3,4,5,6)]
dd$perc_cis_driven <- 100*(dd$`Cis-driven`/dd$`Ancestry:DMP:eGpG`)
dd$perc_not_cis_driven <- 100*(dd$`Not_cis-driven`/dd$`Ancestry:DMP:eGpG`)
round(mean(dd$perc_cis_driven), 2)

plot(sapply(tissues, function(tissue) nrow(metadata[[tissue]])),
     dd$`Cis-driven`/dd$`Ancestry:DMP:eGpG`,
     col = tissue_info$colcodes,
     pch = 16)
cor.test(sapply(tissues, function(tissue) nrow(metadata[[tissue]])),
         dd$`Cis-driven`/dd$`Ancestry:DMP:eGpG`)
cor.test(sapply(tissues, function(tissue) nrow(metadata[[tissue]])),
         dd$`Cis-driven`/dd$`Ancestry:DMP:eGpG`, method = "spearman")

# Data for bar plot --
ddp <- cbind.data.frame("not mGene" = dd$`Ancestry:DMP` - dd$`Ancestry:DMP:eGpG`, # not sGenes
                        "not cis-driven" = dd$`Not_cis-driven`,
                        "cis-driven" = dd$`Cis-driven`,
                        "not classified" = dd$`Too_Many_mQTL` # ot modelled: No isQTL with MAF 001, no isQTL with var, no isQTL with no dependance
) 
rownames(ddp) <- tissue_info[tissues,'tissue_abbrv']
splic.ddp <- ddp
d2 <- as.data.frame(sapply(1:4, function(i) apply(splic.ddp, 1, function(x) x[i]/sum(x))))
colnames(d2) <- colnames(splic.ddp)
rownames(splic.ddp) <- tissue_info$tissue_abbrv
splic.ddp$tissue <- rownames(ddp)
library(reshape2)
d2$tissue <- tissue_info[tissues,'tissue_abbrv']
d1 <- melt(splic.ddp)
colnames(d1)[3] <- "count"
d11 <- melt(d2)
colnames(d11)[3] <- "proportion"
cis_data <- cbind.data.frame(d1, d11[,3])
colnames(cis_data)[4] <- "proportion"
cis_data$proportion <- 100*cis_data$proportion

cis_data$tissue <- factor(cis_data$tissue, levels = rev(tissue_info[tissues,'tissue_abbrv']), order = T)
cis_data$variable <- factor(cis_data$variable, levels = rev(levels(cis_data$variable)), order = T)
cols <- c(brewer.pal(11,"PRGn")[c(2)], "Red",brewer.pal(11,"PRGn")[c(9)],"light grey")

pdf("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/ancestry_DMPs.cis_driven.bar_plot.pdf",
    width = 8, height = 5)
names(cols) <- c("not mGene", "not cis-driven", "cis-driven", "not classified")
ggplot(cis_data, mapping = aes(x = 100*proportion, 
                               fill = variable, #actor(var_2, levels = rev(x_levs)), 
                               y = tissue, #factor(var_1, levels = rev(y_levs)),
                               label = count)) +
  geom_bar( stat = "identity") + 
  geom_text(size = 3, position = position_stack(vjust = 0.5)) + # Add number labelsfacet_grid(~trait) +
  ylab("") + xlab("DMPs (%)") +
  labs(fill = "") +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10))
dev.off()  

### dotplot ####
head(dd)
rownames(dd) <- tissue_info[tissues,'tissue_abbrv']
dd$perc_cis_driven_2 <- 100*(dd$`Cis-driven`/(dd$`Cis-driven`+dd$`Not_cis-driven`))
library(reshape2)
dd$tissue <- tissue_info[tissues,'tissue_abbrv']
d1 <- melt(dd)

tissues_cols <- tissue_info[, 3]
names(tissues_cols) <- tissue_info[, 1]
tissues_cols <- tissues_cols[tissues]
names(tissues_cols) <- tissue_info$tissue_abbrv

pdf('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/cis_driven_DMPs_dotplot.pdf', width = 1.2, height = 2)
ggplot(data = d1[d1$variable %in% c('perc_cis_driven_2'),],
             aes(x = variable,
                 y = value,
                 color =  tissue),
) +geom_boxplot(col = "black",
               fill = "white",
               outlier.shape = NA,
               width = 0.8) +
  geom_jitter(alpha = 1,
              size = 2) +
  xlab("") +
  ylab("% cis-driven DMPs") +
  ylim(c(0,100.5))+
  scale_x_discrete(breaks=c("perc_cis_driven_2"),
                   labels=c(""))+
  scale_fill_manual(values = tissues_cols) +
  scale_color_manual(values = tissues_cols) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.title.y = ggtext::element_markdown(size=12),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none") 
dev.off()
#ggsave('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/cis_driven_DMPs_dotplot.pdf', plot = p1, width = 4, height = 5)
# % of DMPs with cis-mQTLs
to_plot <- dd[,c(1,3)]
to_plot$percentage <- to_plot[,2]/to_plot[,1]*100
to_plot$Tissue <- rownames(to_plot)
to_plot$dummy <- "Dummy"
pdf('Plots/cis_driven_DMPs_with_mQTLs.pdf', width = 1.2, height = 2)
ggplot(data = to_plot,
       aes(x = dummy,
           y = percentage,
           color =  Tissue),
) +geom_boxplot(col = "black",
                fill = "white",
                outlier.shape = NA,
                width = 0.8) +
  geom_jitter(alpha = 1,
              size = 2) +
  xlab("") +
  ylab("% DVPs with mQTLs") +
  ylim(c(0, 100.5)) +
  scale_x_discrete(breaks=c("Dummy"),
                   labels=c(""))+
  scale_fill_manual(values = tissues_cols) +
  scale_color_manual(values = tissues_cols) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.title.y = ggtext::element_markdown(size=12),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none") 
dev.off()

# 1.5 Tissue sharing ----
deg_sharing <- readRDS("~/marenostrum/Projects/GTEx_v8/Methylation/Data/Sharing_DMP.rds")
deg_sharing <- deg_sharing[deg_sharing$trait=='EURv1',]

res_all <- do.call('rbind.data.frame', results)
cis_driven <- unique(res_all$deg[res_all$Class=='Cis-driven'])
deg_sharing$Type <- 'Not Cis-driven'
deg_sharing$Type[deg_sharing$CG %in% unique(res_all$deg)] <- 'Cis-mCpG'
deg_sharing$Type[deg_sharing$CG %in% cis_driven] <- 'Cis-driven'

### Plot jitter for trait specific results #####
# to_plot$type <- 'Hyper'
# to_plot$type[to_plot$dir==-1] <- 'Hypo'

colors <- c('#A39A92','#9C731C')
to_plot$label <- 'No'
to_plot$label[to_plot$number >= 7] <- 'Yes'

library(ggrepel)
#pdf('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/tissue_sharing_ancestry_violin.pdf', width = 3, height = 3)
#png('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/tissue_sharing_ancestry_violin.png', width = 900, height = 700, res = 300)
p <- (ggplot(deg_sharing[deg_sharing$Type %in% c("Cis-driven",'Not Cis-driven'),]) + #geom_jitter(aes(Type, number, col=Type), alpha=0.5) +
        theme_bw() + #geom_violin(aes(x = type, y = number, fill=type),alpha=0.8) +
        geom_boxplot(aes(Type, number),col = "black",
                     fill = "white",
                     outlier.shape = NA,
                     notch = T,
                     width = 0.25) +
        geom_violin(aes(Type, number,fill=Type),col = "black")+
        scale_fill_manual(values = rev(colors)) +
        scale_color_manual(values = rev(colors)) + xlab('') + ylab('Nº Tissues') +
        scale_x_discrete(breaks=c("Cis-driven",'Not Cis-driven'),
                         labels=c("Cis-\nDriven",'Not cis\nmCpG'))+
        #ggrepel::geom_text_repel(aes(type, number,label=ifelse(label=='Yes',as.character(UCSC_RefGene_Name),'')), max.overlaps = Inf, )+
        theme(legend.position = "none",
              axis.text = element_text(size = 9),
              axis.title = element_text(size = 9)))
#dev.off()
ggsave('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/tissue_sharing_ancestry_violin_cis_driven.pdf', p, dpi = 300, width = 1.5, height = 2, units = 'in')


### functional enrichment cis-driven cpgs #####
head(res_all)
res_all$tissue <- gsub('\\..*','',rownames(res_all))
head(results_DML_all$Lung)
GO_res <- list()
library(missMethyl)
for (tissue in tissues) {
  res <- res_all[res_all$Class=='Cis-driven' & res_all$tissue==tissue,]
  res <- res[!is.na(res$deg),]
  if (nrow(res) < 1 ) {
    next
  }
  tryCatch(
    {res <- gometh(as.character(res$deg), all.cpg=rownames(results_DML_all[[tissue]]),
                   collection="GO", array.type="EPIC")
    res <- res[res$ONTOLOGY=="BP",]
    print(table(res$FDR<0.05))
    GO_res[[tissue]] <- res
    },  error=function(cond) {
      message("Error")
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    })
}

#### enrichment chromhmm #####
files <- list.files('~/marenostrum/Projects/GTEx_v8/Methylation/Data/EpiMap/', pattern='.bed.gz',full.names=T)

names_chrom <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "Breast", "MuscleSkeletal", "KidneyCortex", "Testis", "PBMC")
chromhmm <- lapply(names_chrom, function(tis)
  read.delim(files[grep(tis, files)], sep='\t', header=F))
names(chromhmm) <- tissues

library(dplyr)
annotation <- read.csv("~/marenostrum_scratch/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")
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
  print(tissue)
  res <- res_all[res_all$tissue==tissue,]
  chrom_tissue <- as.data.frame(chromhmm_cpgs[[tissue]])
  chrom_tissue$region_chromhmm_new <- chrom_tissue$region_chromhmm
  chrom_tissue <- chrom_tissue[chrom_tissue$.overlap ==1,]
  # 
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TssFlnkD", "TssFlnk", "TssFlnkU","TssA")] <- "TSS"
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("EnhA2", "EnhA1","EnhWk","EnhG1", "EnhG2")] <- "Enh"
  
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("ReprPCWk","ReprPC")] <- "ReprPC"
  
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TxWk","Tx")] <- "Tx"
  
  type_df <- chrom_tissue[chrom_tissue$region_chromhmm_new == type,]# & chrom_tissue$name_ann %in% res$deg,]
  other_type <- chrom_tissue[chrom_tissue$region_chromhmm_new != type,]# & chrom_tissue$name_ann %in% res$deg,]
  type_diff <- nrow(type_df[type_df$name_ann %in% res$deg[res$Class == 'Cis-driven'],])
  type_notdiff <- nrow(type_df) - type_diff
  other_type_diff <- nrow(other_type[other_type$name_ann %in% res$deg[res$Class == 'Cis-driven'],])
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
#families <- as.vector(unique(chromhmm_cpgs$Lung$region_chromhmm))
families <- c('Enh','EnhBiv','Het','Quies','ReprPC','TSS','TssBiv','Tx','ZNF/Rpts')
fisher_results <- lapply(tissues, function(tis) lapply(families, function(region) my_fisher(region, tis)))
names(fisher_results) <- tissues

for (name in tissues) {
  names(fisher_results[[name]]) <- families
}

saveRDS(fisher_results, '~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_cis_driven_vs_all.rds')

## plot together in a heatmap
read_data <- function(variables, data, tissue){ #Function to prepare data to plot and compute adjusted p value
  
  odds_ratio <- lapply(variables, function(type) data[[tissue]][[type]][['f']]$estimate)
  adj.P.Val <- p.adjust(sapply(variables, function(type) data[[tissue]][[type]][['f']]$p.value), method = "BH")
  CI_down <- lapply(variables, function(type) data[[tissue]][[type]][['f']]$conf.int[1])
  CI_up <- lapply(variables, function(type) data[[tissue]][[type]][['f']]$conf.int[2])
  sample_size <- lapply(variables, function(type) data[[tissue]][[type]][['m']])
  
  
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

fisher_results <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_cis_driven_vs_all.rds')
hyper_all <- lapply( tissues[!tissues %in% c('Testis','BreastMammaryTissue')], function(x) read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), fisher_results, x))
names(hyper_all) <- tissues[!tissues %in% c('Testis','BreastMammaryTissue')]
hyper_all <- do.call('rbind.data.frame', hyper_all)

hyper_all$Tissue <- gsub('\\..*','',rownames(hyper_all))

pdf('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/heatmap_cisdriven_enrichment_ancestry.chromhmm.pdf', width = 5, height = 4)
ggplot(hyper_all, aes(x=Tissue, y=region, size=-log10(adjPvalue), fill=log2(oddsRatio))) +
  geom_point(alpha=0.9, shape=21, color="black") +
  scale_size(range = c(0.5, 12), name="-log10(FDR)", limits = c(0,max(-log10(hyper_all$adjPvalue)))) +
  scale_fill_gradient2(low="#1982C4", high="#690500",  mid = '#F4EBBE', midpoint = 0.7,name='log2(oddsRatio)')+
  #theme_ipsum() +
  theme(legend.position="right") +
  ylab("") +
  xlab("") +
  scale_x_discrete(breaks= tissues[!tissues %in% c('Testis','BreastMammaryTissue')],
                   labels=tissue_info[tissues[!tissues %in% c('Testis','BreastMammaryTissue')],'tissue_abbrv'])+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        panel.grid = element_line(colour = 'light grey'),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
dev.off()

### enrichment cis vs not cis ####
read_data <- function(variables, data, type){ #Function to prepare data to plot and compute adjusted p value
  
  # if (tissue %in% sex_tissues & trait == "SEX2") {
  #   print(NA)
  # } else {
  
  odds_ratio <- lapply(variables, function(tissue) data[[tissue]][[type]][['f']]$estimate)
  adj.P.Val <- p.adjust(sapply(variables, function(tissue) data[[tissue]][[type]][['f']]$p.value), method = "BH")
  CI_down <- lapply(variables, function(tissue) data[[tissue]][[type]][['f']]$conf.int[1])
  CI_up <- lapply(variables, function(tissue) data[[tissue]][[type]][['f']]$conf.int[2])
  sample_size <- lapply(variables, function(tissue) data[[tissue]][[type]][['m']])
  
  
  names(odds_ratio) <- variables
  names(adj.P.Val) <- variables
  names(CI_down) <- variables
  names(CI_up) <- variables
  names(sample_size) <- variables
  
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
  
  sample_size_df <- as.data.frame(unlist(sample_size))
  sample_size_df$label <- variables
  sample_size_df$type <- deparse(substitute(data))
  colnames(sample_size_df) <- c('sample_size','tissue','type')
  
  all <- Reduce(function(x, y) merge(x, y, all=TRUE), list(odds_ratio_df, adj.P.Val_df, CI_down_df, CI_up_df, sample_size_df))
  head(all)
  all$sig <- 'not Sig'
  all$sig[all$adjPvalue<0.05] <- 'Sig'
  all <- all[,c("tissue","oddsRatio","adjPvalue","CI_down","CI_up","sig","type", "sample_size")]
  return(all)
}

##### plot shared positions
female <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_cis_driven_vs_all.rds')
male <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_cis_not_driven_vs_all.rds')

library(ggpubr)

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "KidneyCortex", "Testis", "WholeBlood","MuscleSkeletal")

library(ggpubr)
sex_tissues <- c('Ovary','Prostate','Testis')
for (type in c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts")) {
    hypo_d <- read_data(tissues, male, type)
    hyper_d <- read_data(tissues, female, type)
    hyper_hypo <- rbind(hypo_d, hyper_d)
    hyper_hypo$sig[hyper_hypo$sig =="Sig"] <- "FDR < 0.05"
    hyper_hypo$sig[hyper_hypo$sig =="not Sig"] <- "FDR >= 0.05"
    hyper_hypo$sig <- factor(hyper_hypo$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
    hyper_hypo$tissue <- factor(hyper_hypo$tissue, levels=tissues)
    hyper_hypo$type[hyper_hypo$type =="male"] <- "not-cis-driven"
    hyper_hypo$type[hyper_hypo$type =="female"] <- "cis-driven"
    g <- ggplot(hyper_hypo, aes(x=log2(oddsRatio), y=tissue, colour=type, alpha=sig)) +
      geom_errorbar(aes(xmin=log2(CI_down), xmax=log2(CI_up)), width=.3) +
      geom_vline(xintercept = 0) +
      #xlim(0,20) + #Only for Lung to show the 0
      geom_point(size=3) + ylab('') + theme_bw() +
      scale_colour_manual(values= c('#A39A92','#F0AE21')) +
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
            panel.border = element_rect(colour = "black", linewidth=1)) #+
    # scale_y_discrete(breaks=c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"),
    #                  labels=c("Enhancer","Enhancer Bivalent","Heterochromatin","Quiescent","Repressed Polycomb","TSS","TSS Bivalent","Transcription","ZNF & Repeats"))# + xlim(0, 3)
    # pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/genomic_location_", gsub('\\/','_',type),'_',trait,".pdf"), w = 6, h = 3.5)
    # print(g)
    # dev.off()
    
    #Plot sample sizes:
    
    g2 <- ggplot(hyper_hypo) + geom_col(aes(sample_size, tissue, fill=type), width = 0.6) +
      theme_classic() + xlab("Number of DMPs") + ylab("") +
      scale_fill_manual(values= c('#A39A92','#F0AE21')) +
      theme(legend.position = "none",
            axis.text.x = element_text(colour="black", size=13),
            axis.text.y=element_blank(),  #remove y axis labels,
            axis.title.x = element_text(size=16)) +
      scale_x_continuous(n.breaks=3)
    #   scale_y_discrete(breaks=c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"),
    #                    labels=c("Enhancer","Enhancer Bivalent","Heterochromatin","Quiescent","Repressed Polycomb","TSS","TSS Bivalent","Transcription","ZNF & Repeats")) #+
    # #scale_x_continuous(breaks=c(0, 20000, 40000)) #Only for lung
    # pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/genomic_location_", gsub('\\/','_',type),'_',trait,"_sample_size.pdf"), w = 4, h = 3.5)
    # print(g2)
    # dev.off()
    
    p <- ggarrange(g, g2, labels = c("A", "B"),
                   common.legend = TRUE, legend = "right", widths = c(0.8,0.3))
    pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/enrichment", gsub('\\/','_',type),".all.cis_driven.pdf"), w = 8, h = 4)
    print(p)
    dev.off()
  }



### Enrichment cis-driven 
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
  colnames(odds_ratio_df) <- c('oddsRatio', 'tissue')
  
  adj.P.Val_df <- as.data.frame(unlist(adj.P.Val))
  adj.P.Val_df$label <- variables
  colnames(adj.P.Val_df) <- c('adjPvalue','tissue')
  
  CI_down_df <- as.data.frame(unlist(CI_down))
  CI_down_df$label <- variables
  colnames(CI_down_df) <- c('CI_down','tissue')
  
  CI_up_df <- as.data.frame(unlist(CI_up))
  CI_up_df$label <- variables
  colnames(CI_up_df) <- c('CI_up','tissue')
  
  sample_size_df <- as.data.frame(unlist(sample_size))
  sample_size_df$label <- variables
  colnames(sample_size_df) <- c('sample_size','tissue')
  
  all <- Reduce(function(x, y) merge(x, y, all=TRUE), list(odds_ratio_df, adj.P.Val_df, CI_down_df, CI_up_df, sample_size_df))
  head(all)
  all$sig <- 'not Sig'
  all$sig[all$adjPvalue<0.05] <- 'Sig'
  all <- all[,c("tissue","oddsRatio","adjPvalue","CI_down","CI_up","sig", "sample_size")]
  return(all)
}
cis_driven_DMPs <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_mcpgs_dmps.fisher.rds')
hyper_all <- read_data(tissues[!tissues %in% c('Testis','BreastMammaryTissue')],cis_driven_DMPs)

hyper_all$sig[hyper_all$sig =="Sig"] <- "FDR < 0.05"
hyper_all$sig[hyper_all$sig =="not Sig"] <- "FDR >= 0.05"
hyper_all$sig <- factor(hyper_all$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
hyper_all$tissue <- factor(hyper_all$tissue, levels=(c('MuscleSkeletal','WholeBlood','KidneyCortex','Prostate','Ovary','ColonTransverse','Lung')))

colors <- c('#A39A92','#F0AE21')
g <- ggplot(hyper_all, aes(x=(oddsRatio), y=tissue, colour=sig)) +
  geom_errorbar(aes(xmin=(CI_down), xmax=(CI_up)), width=.3) +
  geom_vline(xintercept = 1) +
  #xlim(0,8) + #Only for Lung to show the 0
  geom_point(size=3) + ylab('') + theme_bw() +
  scale_colour_manual(values=(colors)) +
  xlab("Odds ratio") +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(colour="black", size=13),
        axis.text.y = element_text(colour="black", size=14),
        legend.text = element_text(colour="black", size=13),
        axis.title.x = element_text(size=16),
        legend.spacing.y = unit(-0.05, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth=1))
# pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Ovary_menopause/genomic_location_", tissue,'_',trait,".pdf"), w = 6, h = 3.5)
# print(g)
# dev.off()

#Plot sample sizes:

g2 <- ggplot(hyper_all) + geom_col(aes(sample_size, tissue, fill=sig), width = 0.6) +
  theme_classic() + xlab("Number of DMPs") + ylab("") +
  scale_fill_manual(values=(colors)) +
  theme(legend.position = "none",
        axis.text.x = element_text(colour="black", size=13),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=16)) +
  scale_x_continuous(n.breaks=3)
#scale_x_continuous(breaks=c(0, 20000, 40000)) #Only for lung
library(ggpubr)
p <- ggarrange(g, g2, labels = c("A", "B"),
               common.legend = TRUE, legend = "right", widths = c(0.8,0.3))
pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/enrichment_cis_driven.pdf"), w = 8, h = 4)
print(p)
dev.off()

#### cis-driven vs not driven enrichment #####
