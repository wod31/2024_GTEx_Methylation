#### Enrichment of enhancers demethylated in colon, Lung and Ovary

library(minfi)
# library(lumi)
library(dplyr)
library(tidyr)
library(tidyverse)
library(missMethyl)
## Input data #####
#library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

#first_dir <- "/gpfs/projects/bsc83/"
first_dir <- "~/marenostrum/"
#annotation <- read.csv("/gpfs/scratch/bsc83/bsc83535/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")
annotation <- read.csv("~/marenostrum_scratch/MN4/bsc83/bsc83535/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")

project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary")
#tissues <- c("MuscleSkeletal","Testis",'KidneyCortex',"BreastMammaryTissue","Prostate","Lung", "ColonTransverse", "Ovary", "WholeBlood")
tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis","WholeBlood")

names <- c("Age")
names <- c("Ancestry")
traits_to_use <- c('AGE')
traits_to_use <- c('EURv1')

results_DML <- lapply(tissues, function(tis) 
  readRDS(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/",tis,"/DML_results_5_PEERs_continous.rds")))
names(results_DML) <- tissues

files <- list.files('~/marenostrum/Projects/GTEx_v8/Methylation/Data/EpiMap/', pattern='.bed.gz',full.names=T)

names_chrom <- c("Lung", "ColonTransverse", "Ovary")
names_chrom <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "Breast", "MuscleSkeletal", "KidneyCortex", "Testis","PBMC")

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

enrichment <- function(type, tissue, trait){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  print(tissue)
  res <- results_DML[[tissue]][[trait]]
  chrom_tissue <- as.data.frame(chromhmm_cpgs[[tissue]])
  #chrom_tissue$region_chromhmm <- as.factor(chrom_tissue$region_chromhmm)
  chrom_tissue$region_chromhmm_new <- chrom_tissue$region_chromhmm
  chrom_tissue <- chrom_tissue[chrom_tissue$.overlap ==1,]
  # 
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TssFlnkD", "TssFlnk", "TssFlnkU","TssA")] <- "TSS"
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("EnhA2", "EnhA1","EnhWk","EnhG1", "EnhG2")] <- "Enh"
  
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("ReprPCWk","ReprPC")] <- "ReprPC"
  
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TxWk","Tx")] <- "Tx"
  
  type_df <- chrom_tissue[chrom_tissue$region_chromhmm_new == type & chrom_tissue$name_ann %in% rownames(res[res$adj.P.Val<0.05,]),]
  ## associated to gene
  #print(table(ann$Type[ann$IlmnID %in% (type_df$name_ann)]))
  #type_df <- type_df[type_df$name_ann %in% ann$IlmnID,]
  
    if (nrow(type_df) < 1 ) {
      return(NA)
    }
    tryCatch(
      {res <- gometh(type_df$name_ann, all.cpg=chrom_tissue$name_ann,#chrom_tissue$name_ann[chrom_tissue$region_chromhmm_new == type],
                     collection="GO", array.type="EPIC", sig.genes = TRUE)
      res <- res#[res$ONTOLOGY=="BP",]
      print(table(res$FDR<0.05))
      return(res)
      },  error=function(cond) {
        message("Error")
        message("Here's the original error message:")
        message(cond)
        # Choose a return value in case of error
        return(NA)
      })
}

### run
families <- c('Enh','ReprPC','TssBiv','EnhBiv','Quies')
families <- c('Enh')
families <- c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts")
#ann <- read.delim('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Methylation_Epic_gene_promoter_enhancer_processed.txt')

fisher_results <- lapply(tissues, function(tis) lapply(traits_to_use, function(trait) lapply(families, function(region) enrichment(region, tis,trait))))
names(fisher_results) <- tissues

for (name in tissues) {
  names(fisher_results[[name]]) <- traits_to_use
}

for (name in tissues) {
  for (trait in names(fisher_results[[name]])) {
    names(fisher_results[[name]][[trait]]) <- families
  }
}

##all final table
final_table_all <- data.frame(name="X",term="X",ONTOLOGY='X', GeneRatio=1, FDR=1, P.DE=1, trait=1, tissue=1, loc=1, genes="X")
for (tissue in tissues) {
  print(tissue)
  probes <- data.frame(name="X",term="X",ONTOLOGY='X',GeneRatio=1, FDR=1, P.DE=1,trait=1, tissue=1, loc=1, genes="X")
  for (trait in traits_to_use) {
    for (fam in families) {
      res <- as.data.frame(fisher_results[[tissue]][[trait]][[fam]])
      if (nrow(res)<2) {
        next}
      res$name <- rownames(res)
      ###for GO 
      res$term <- res$TERM
      res$tissue <- tissue
      res$trait <- trait
      res$loc <- fam
      res$GeneRatio <- res$DE/res$N
      res$genes <- res$SigGenesInSet
      probes <- rbind(probes, res[,c('name','term','ONTOLOGY','GeneRatio','FDR','P.DE','trait','tissue','loc','genes')])
    }
  }
  probes <- probes[-1,]
  final_table_all <- rbind(final_table_all, probes)
}

final_table_all <- final_table_all[-1,]

head(final_table_all)
final_table_all <- final_table_all[final_table_all$P.DE<0.05,]
hyper <- final_table_all[final_table_all$FDR<0.05,]
hyper$direction <- 'Hypemerthylation'

hypo <- final_table_all[final_table_all$FDR<0.05,]
hypo$direction <- 'Hypomethylation'

all <- rbind(hyper, hypo)

write.table(all, '~/marenostrum/Projects/GTEx_v8/Methylation/Data/SupplementaryTablex_enrichment_chrom_state_aging.txt',
            sep = '\t', row.names = F, col.names = T, quote = F)

write.table(hyper, '~/marenostrum/Projects/GTEx_v8/Methylation/Data/SupplementaryTable4_enrichment_enhancers_ancestry.txt',
            sep = '\t', row.names = F, col.names = T, quote = F)


### where are quiescent located? gene body? enhancer-associated?
ann <- read.delim('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Methylation_Epic_gene_promoter_enhancer_processed.txt')
my_fisher <- function(type){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  
  type_df <- res[res$Type == type,]
  other_type <-  res[res$Type != type,]
  type_diff <- nrow(res_df[rownames(res_df) %in% rownames(type_df),])
  type_notdiff <- nrow(type_df) - type_diff
  other_type_diff <- nrow(res_df[rownames(res_df) %in% rownames(other_type),])
  other_type_notdiff <- nrow(other_type) - other_type_diff
  
  ### test significance
  m <- matrix(c(type_diff, other_type_diff, type_notdiff, other_type_notdiff), 2,2, byrow = T)
  print(m)  
  
  m[is.na(m)] <- 0
  #m <- m[c(type,paste0('No ',type)),]
  rownames(m) <- c('DMP', "Not DMP")
  colnames(m) <- c(type,"Other")
  print(m)
  f <- fisher.test(m)
  print(f)
  return(list("f" = f, "m" = type_diff))
}

res_fisher <- list()
for (tissue in tissues) {
  print(tissue)
  res <- results_DML[[tissue]][['EURv1']]
  chrom_tissue <- as.data.frame(chromhmm_cpgs[[tissue]])
  chrom_tissue <- chrom_tissue[chrom_tissue$.overlap ==1,]
  res_df <- res[rownames(res) %in% chrom_tissue$name_ann[chrom_tissue$region_chromhmm == 'Quies'] & res$adj.P.Val<0.05,]
  res <- res[rownames(res) %in% chrom_tissue$name_ann[chrom_tissue$region_chromhmm == 'Quies'],]
  print(paste0('Nº of DMPs in Quies',nrow(res_df)))
  print(table(ann$Type[ann$IlmnID %in% rownames(res_df)]))
  res$Type <- 'Other'
  res$Type[rownames(res) %in% ann$IlmnID[ann$Type=="Promoter_Associated"]] <- "Promoter_Associated"
  res$Type[rownames(res) %in% ann$IlmnID[ann$Type=="Enhancer_Associated"]] <- "Enhancer_Associated"
  res$Type[rownames(res) %in% ann$IlmnID[ann$Type=="Gene_Associated"]] <-"Gene_Associated"
  table(res$Type)
  #enrichment
  families <- unique(res$Type)
  fisher_results <- lapply(families, function(region) my_fisher(region))
  names(fisher_results) <- families
  res_fisher[[tissue]] <- fisher_results
  
}
read_data <- function(tissues, data,type){ #Function to prepare data to plot and compute adjusted p value
  
  odds_ratio <- lapply(tissues, function(tissue) data[[tissue]][[type]][['f']]$estimate)
  adj.P.Val <- p.adjust(sapply(tissues, function(tissue) data[[tissue]][[type]][['f']]$p.value), method = "BH")
  CI_down <- lapply(tissues, function(tissue) data[[tissue]][[type]][['f']]$conf.int[1])
  CI_up <- lapply(tissues, function(tissue) data[[tissue]][[type]][['f']]$conf.int[2])
  sample_size <- lapply(tissues, function(tissue) data[[tissue]][[type]][['m']])
  
  names(odds_ratio) <- tissues
  names(adj.P.Val) <- tissues
  names(CI_down) <- tissues
  names(CI_up) <- tissues
  names(sample_size) <- tissues
  
  odds_ratio_df <- as.data.frame(unlist(odds_ratio))
  odds_ratio_df$label <- tissues
  odds_ratio_df$type <- deparse(substitute(data)) #Either hypo or hyper
  colnames(odds_ratio_df) <- c('oddsRatio', 'tissue','type')
  
  adj.P.Val_df <- as.data.frame(unlist(adj.P.Val))
  adj.P.Val_df$label <- tissues
  adj.P.Val_df$type <- deparse(substitute(data))
  colnames(adj.P.Val_df) <- c('adjPvalue','tissue','type')
  
  CI_down_df <- as.data.frame(unlist(CI_down))
  CI_down_df$label <- tissues
  CI_down_df$type <- deparse(substitute(data))
  colnames(CI_down_df) <- c('CI_down','tissue','type')
  
  CI_up_df <- as.data.frame(unlist(CI_up))
  CI_up_df$label <- tissues
  CI_up_df$type <- deparse(substitute(data))
  colnames(CI_up_df) <- c('CI_up','tissue','type')
  
  sample_size_df <- as.data.frame(unlist(sample_size))
  sample_size_df$label <- tissues
  sample_size_df$type <- deparse(substitute(data))
  colnames(sample_size_df) <- c('sample_size','tissue','type')
  
  all <- Reduce(function(x, y) merge(x, y, all=TRUE), list(odds_ratio_df, adj.P.Val_df, CI_down_df, CI_up_df, sample_size_df))
  head(all)
  all$sig <- 'not Sig'
  all$sig[all$adjPvalue<0.05] <- 'Sig'
  all <- all[,c("tissue","oddsRatio","adjPvalue","CI_down","CI_up","sig","type", "sample_size")]
  return(all)
}

results_fisher <- lapply(unique(res$Type), function(type) read_data(tissues, res_fisher, type))
names(results_fisher) <- unique(res$Type)

## plot numbers
numbers_quiesc <- list()
for (tissue in tissues) {
  print(tissue)
  res <- results_DML[[tissue]][['EURv1']]
  chrom_tissue <- as.data.frame(chromhmm_cpgs[[tissue]])
  chrom_tissue <- chrom_tissue[chrom_tissue$.overlap ==1,]
  res_df <- res[rownames(res) %in% chrom_tissue$name_ann[chrom_tissue$region_chromhmm == 'Quies'] & res$adj.P.Val<0.05,]
  #res <- res[rownames(res) %in% chrom_tissue$name_ann[chrom_tissue$region_chromhmm == 'Quies'],]
  print(paste0('Nº of DMPs in Quies',nrow(res_df)))
  print(table(ann$Type[ann$IlmnID %in% rownames(res_df)]))
  res_df$Type <- 'Other'
  res_df$Type[rownames(res_df) %in% ann$IlmnID[ann$Type=="Gene_Associated"]] <-"Gene_Associated"
  res_df$Type[rownames(res_df) %in% ann$IlmnID[ann$Type=="Enhancer_Associated"]] <- "Enhancer_Associated"
  res_df$Type[rownames(res_df) %in% ann$IlmnID[ann$Type=="Promoter_Associated"]] <- "Promoter_Associated"
  numbers_quiesc[[tissue]]<- as.data.frame(table(res_df$Type))
  
}
numbers_quiesc <- do.call('rbind.data.frame',numbers_quiesc)
numbers_quiesc$Var1 <- factor(numbers_quiesc$Var1, levels = c('Gene_Associated','Enhancer_Associated','Promoter_Associated','Other'))
numbers_quiesc$tissue <- gsub('\\..*','',rownames(numbers_quiesc))
colors_types <- c('#FB6107','#D1AC00','#F1E07E','#69747C')
pdf('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/quiescent_location_ancestry.pdf', width = 5, height = 3)
ggplot(numbers_quiesc, aes(y = tissue,x=Freq, fill=Var1)) +
  geom_bar(position = 'fill', alpha=0.8, stat = 'identity') + 
  geom_text(aes(label=Freq, x = Freq+30), stat='identity', position='fill', size=3,hjust=1.6, angle=20) +
  theme_bw() + ylab('') + xlab('Proportion correlated DMPs') + scale_fill_manual(values = colors_types)#+
  #facet_wrap2(~ Trait, strip = strip, nrow = 1) + theme_classic() + scale_x_continuous(breaks=seq(0, 1, 0.5))
dev.off()

#### get genes quiescent #####

genes_quiesc_tissue <- list()
for (tissue in tissues) {
  print(tissue)
  res <- results_DML[[tissue]][['EURv1']]
  chrom_tissue <- as.data.frame(chromhmm_cpgs[[tissue]])
  chrom_tissue <- chrom_tissue[chrom_tissue$.overlap ==1,]
  res_df <- res[rownames(res) %in% chrom_tissue$name_ann[chrom_tissue$region_chromhmm == 'Quies'] & res$adj.P.Val<0.05 & res$logFC>0,]
  res <- res[rownames(res) %in% chrom_tissue$name_ann[chrom_tissue$region_chromhmm == 'Quies'],]
  print(paste0('Nº of DMPs in Quies',nrow(res_df)))
  print(table(ann$Type[ann$IlmnID %in% rownames(res_df)]))

  genes_quiesc_tissue[[tissue]]<- unique(ann$UCSC_RefGene_Name[ann$IlmnID %in% rownames(res_df)])
  write.table(genes_quiesc_tissue[[tissue]], paste0('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Genes_Quiescent_Ancestry_',tissue,'.txt'),
              sep = '\n', row.names = F, col.names = F, quote = F)
  genes_quiesc <- unique(ann$UCSC_RefGene_Name[ann$IlmnID %in% rownames(res)])
  write.table(genes_quiesc, paste0('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Genes_Quiescent_Ancestry_',tissue,'_bg.txt'),
              sep = '\n', row.names = F, col.names = F, quote = F)
}

genes_quiesc_tissue_all <- unique(unlist(genes_quiesc_tissue))
write.table(genes_quiesc_tissue_all, paste0('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Genes_Quiescent_Ancestry_all','.txt'),
            sep = '\n', row.names = F, col.names = F, quote = F)
write.table(unique(ann$UCSC_RefGene_Name), paste0('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Genes_cpgs','.txt'),
            sep = '\n', row.names = F, col.names = F, quote = F)

