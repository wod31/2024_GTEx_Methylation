#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Get additive effects of Age with other demographic traits and perform plots
# @software version: R=4.2.2

##########################################
#### Additive effects ####################
##########################################

#!/usr/bin/env Rscript

# Libraries ----
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)
library(circlize)

# ---- Data ----- ####
# Individual Traits ----
traits_cols <- c("Age" = "#56B4E9",
                 "Ancestry" = "#E69F00",
                 "Sex" =  "#009E73",
                 "BMI" = "#CC79A7")
traits <- names(traits_cols)

# Tissues ----
### read methylation results ####
first_dir <- "~/marenostrum/"
project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Age", "Ancestry", "BMI", "Sex")

tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Methylation/Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

sex_tissues <- c('Ovary','Prostate','Testis')

tissues <- tissue_info$tissue_ID

get_DMPs <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "SEX2"){
    NA
  }else{
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.rds"))[[trait]]
    #rownames(model[[trait]][model[[trait]]$adj.P.Val<0.05,])
    model
  }
}
DMPs_Res <- lapply(c('EURv1','SEX2','AGE','BMI'), function(trait) lapply(tissues, function(tissue) get_DMPs(tissue, trait)))
names(DMPs_Res) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(DMPs_Res[[trait]]) <- tissues}


# Lists of DEGs ----
get_DMPs <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "Sex"){
    NA
  }else{
    rownames(DMPs_Res[[trait]][[tissue]][DMPs_Res[[trait]][[tissue]]$adj.P.Val<0.05,])
  }
}
DMPs <- lapply(tissues, function(tissue) lapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait) get_DMPs(tissue, trait)))
names(DMPs) <- tissues
for(trait in tissues){names(DMPs[[trait]]) <- c("Ancestry", "Sex", "Age", "BMI")}

# Read tissue sharing? ----


# ---- Analysis ----- ####

# 1. Overlap between pairwise combinations of traits ----
# Genes expressed per tissue --
tested_cpgs <- lapply(tissues, function(tissue) rownames(DMPs_Res[['Age']][[tissue]]))
names(tested_cpgs) <- tissues

# Pairwise combination of traits --
pw.traits <- c("Age:Ancestry",
               "Age:Sex"
)

# List of DEGs per tissue
deg.pw.traits <- lapply(pw.traits, function(pw)
  lapply(tissues, function(tissue) 
    intersect(DMPs[[tissue]][[unlist(strsplit(pw, split = ":"))[[1]]]], DMPs[[tissue]][[unlist(strsplit(pw, split = ":"))[[2]]]])
  )
)
names(deg.pw.traits) <- pw.traits
for(pw in pw.traits){
  names(deg.pw.traits[[pw]]) <- tissues
}

length(unique(unlist(lapply(names(deg.pw.traits), function(pw) unique(unlist(deg.pw.traits[[pw]]))))))

# Number of DEGs
number_of_DEGs <- sapply(pw.traits, function(pw) 
  sapply(tissues, function(tissue)
    length(deg.pw.traits[[pw]][[tissue]])
  )
)

total_unique_pw_DEGs <- sapply(pw.traits, function(pw) 
  length(unique(unlist(lapply(tissues, function(tissue) deg.pw.traits[[pw]][[tissue]]))))
)
total_pw_DEGs <- sapply(pw.traits, function(pw) 
  length(unlist(lapply(tissues, function(tissue) deg.pw.traits[[pw]][[tissue]])))
)
total_unique_pw_DEGs
total_pw_DEGs

# Number of tissues with genes with additive effects ----
apply(number_of_DEGs, 2, function(x) sum(x>0, na.rm=T))


rownames(number_of_DEGs) <- tissue_info$tissue_abbrv
tissues_cols <- tissue_info$colcodes 
names(tissues_cols) <- tissue_info$tissue_abbrv
row_ha_left <- HeatmapAnnotation("Tissues" = rownames(number_of_DEGs),
                                 col = list("Tissues" = tissues_cols),
                                 show_legend = F, 
                                 show_annotation_name = F,
                                 simple_anno_size = unit(0.3,"cm"),
                                 which = "row")

colnames(number_of_DEGs) <- c("Ancestry-age DEGs",
                              "Sex-age DEGs")
column_ha_top <- HeatmapAnnotation("Number of DEGs" = anno_barplot(total_unique_pw_DEGs,
                                                                   gp = gpar(fill = "light gray",
                                                                             col = "light gray"),
                                                                   border = F),
                                   annotation_name_side = "left",
                                   annotation_name_rot = 90,
                                   annotation_name_gp = gpar(fontsize = 8),
                                   height = unit(2.5, "cm"),
                                   gap = unit(0.25,"cm"))
column_ha_bottom <- HeatmapAnnotation("variable1" = sapply(pw.traits, function(i) unlist(strsplit(i, split = ":"))[[1]]),
                                      "variable2" = sapply(pw.traits, function(i) unlist(strsplit(i, split = ":"))[[2]]),
                                      col = list("variable1" = traits_cols,
                                                 "variable2" = traits_cols),
                                      show_legend = F, 
                                      show_annotation_name = F,
                                      simple_anno_size = unit(0.3,"cm"))


ht <- Heatmap(apply(number_of_DEGs, 2, function(x) x/sum(x)),
              col = brewer.pal(8, "BuPu"),
              na_col = "white",
              cluster_rows = F,
              cluster_columns = F,
              name = "DE signal",
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(prettyNum(number_of_DEGs[i, j], big.mark = ","), x, y, gp = gpar(fontsize = 10))},
              row_names_side = "left",
              column_names_side = "bottom",
              left_annotation = row_ha_left,
              row_names_gp = gpar(fontsize = 10),
              bottom_annotation = column_ha_bottom,
              column_names_gp = gpar(fontsize = 10),
              top_annotation = column_ha_top,
              heatmap_legend_param = list(direction = "horizontal"))





# Fisher's exact test to test if there are more DEGs with 2 traits than expected ----
overlap.enrichment.fun <- function(tissue, trait.1, trait.2){
  if(tissue %in% sex_tissues & (trait.1 == "Sex" | trait.2 == "Sex")){
    return(list("overlap" = NA,
                "odds.ratio" = NA,
                "p.value" = NA))
  }else{
    #                   de.trait.2  not.de.trait.2
    # de.trait.1
    # not.de.trait.1
    x11 <- length(intersect(DMPs[[tissue]][[trait.1]], DMPs[[tissue]][[trait.2]]))
    #x11 <- length(deg.pw.traits[[paste0(trait.1,":", trait.2)]][[tissue]])
    x12 <- length(DMPs[[tissue]][[trait.1]][! DMPs[[tissue]][[trait.1]] %in% DMPs[[tissue]][[trait.2]]])
    x21 <- length(DMPs[[tissue]][[trait.2]][! DMPs[[tissue]][[trait.2]] %in% DMPs[[tissue]][[trait.1]]])
    x22 <- length(tested_cpgs[[tissue]][! tested_cpgs[[tissue]] %in% unique(c(DMPs[[tissue]][[trait.1]], DMPs[[tissue]][[trait.2]]))])
    m <- matrix(c(x11,x12,x21,x22),2,2,byrow = T)
    rownames(m) <- c(paste0(trait.1,".deg"), paste0(trait.1, ".not_deg"))
    colnames(m) <- c(paste0(trait.2,".deg"), paste0(trait.2, ".not_deg"))
    #sum(m) == length(exprs.genes[[tissue]])
    f <- fisher.test(m)
    f$estimate
    return(list("overlap" = x11,
                "odds.ratio" = f$estimate,
                "p.value" = f$p.value,
                "counts.matrix" = m,
                "lower_CI" = f$conf.int[1],
                "upper_CI" = f$conf.int[2]))
  }
}


# Enrichment analysis  --
fisher.results <- lapply(pw.traits, function(pw)
  lapply(tissues, function(tissue) 
    overlap.enrichment.fun(tissue, unlist(strsplit(pw, split = ":"))[[1]], unlist(strsplit(pw, split = ":"))[[2]])
  ))
names(fisher.results) <- pw.traits
for(pw in pw.traits){names(fisher.results[[pw]]) <- tissues}
#fisher.results$`Ancestry:BMI`$WholeBlood$counts.matrix

# Enrichment statistics --
p.value <- sapply(pw.traits, function(pw)
  sapply(tissues, function(tissue) 
    fisher.results[[pw]][[tissue]][["p.value"]]
  )
)

# Multiple testing correction by tissue --
fdr <- t(apply(p.value, 1, function(x) p.adjust(x, method = "BH")))
fdr[fdr>=0.05] <- NA

# multiple testing correction across tissues and number of pairwise combinations of traits
adj_p_values <- p.adjust(unlist(p.value), method = "BH")
adjPVal_matrix <- matrix(adj_p_values, 
                         nrow = length(tissues), ncol = length(pw.traits),
                         byrow = F)

# Number of DEGs in overlap --
overlap <- sapply(pw.traits, function(pw)
  sapply(tissues, function(tissue) 
    fisher.results[[pw]][[tissue]][["overlap"]]
  )
)
rownames(overlap) <- tissues

# Odds ratio --
odds.ratio <- sapply(pw.traits, function(pw)
  sapply(tissues, function(tissue) 
    fisher.results[[pw]][[tissue]][["odds.ratio"]]
  )
)
rownames(odds.ratio) <- tissues
or <- odds.ratio
or[is.na(fdr)] <- NA

log2_OR <- log2(or)
range(log2_OR, na.rm = T)
col_fun <- colorRamp2(seq(0,3,length.out=9),
                      (brewer.pal(9, "Reds")))

#log2_OR[log2_OR==0] <- NA
#number_of_DEGs <- number_of_DEGs[,pw.traits[c(1,4,3,2,5,6)]]
#number_of_DEGs[number_of_DEGs==0] <- ""
colnames(log2_OR) <- c("Ancestry-age DEGs",
                       "Sex-age DEGs")
rownames(log2_OR) <- tissue_info$tissue_abbrv
tissues_cols <- tissue_info$colcodes 
names(tissues_cols) <- tissue_info$tissue_abbrv
row_ha_left <- HeatmapAnnotation("Tissues" = rownames(log2_OR),
                                 col = list("Tissues" = tissues_cols),
                                 show_legend = F, 
                                 show_annotation_name = F,
                                 simple_anno_size = unit(0.3,"cm"),
                                 which = "row")

m <- as.matrix(rbind(apply(or, 2, function(x) sum(x<1, na.rm=T)),
                     apply(or, 2, function(x) sum(x>1, na.rm=T))))

column_ha_top_ht2 <- HeatmapAnnotation("Number of tissues" = anno_barplot(t(m),
                                                                          gp = gpar(fill = brewer.pal(11, "RdBu")[c(10,2)],
                                                                                    col = brewer.pal(11, "RdBu")[c(10,2)]),
                                                                          border = F),
                                       annotation_name_side = "right",
                                       annotation_name_rot = 90,
                                       annotation_name_gp = gpar(fontsize = 8),
                                       height = unit(2.5, "cm"),
                                       gap = unit(0.25,"cm"))

pw.traits[1] <- "Ancestry:Age"
pw.traits[2] <- "Sex:Age"
column_ha_bottom <- HeatmapAnnotation("variable1" = sapply(pw.traits, function(i) unlist(strsplit(i, split = ":"))[[1]]),
                                      "variable2" = sapply(pw.traits, function(i) unlist(strsplit(i, split = ":"))[[2]]),
                                      col = list("variable1" = traits_cols,
                                                 "variable2" = traits_cols),
                                      show_legend = F, 
                                      show_annotation_name = F,
                                      simple_anno_size = unit(0.3,"cm"))

ht2 <- Heatmap(log2_OR,
               col = col_fun,
               na_col = "white",
               cluster_rows = F,
               cluster_columns = F,
               name = "Odds ratio (log2)",
               #cell_fun = function(j, i, x, y, width, height, fill) {
               #  grid.text(round(-log10(fdr[i, j]), 2), x, y, gp = gpar(fontsize = 10))},
               row_names_side = "left",
               column_names_side = "bottom",
               #left_annotation = row_ha_left,
               row_names_gp = gpar(fontsize = 10),
               bottom_annotation = column_ha_bottom,
               column_names_gp = gpar(fontsize = 10),
               top_annotation = column_ha_top_ht2,
               heatmap_legend_param = list(direction = "horizontal"))


log10_fdr_sign <- -log10(fdr)
for(i in 1:nrow(log10_fdr_sign)){
  for(j in 1:ncol(log10_fdr_sign)){
    if(!is.na(log2_OR[i,j]) & log2_OR[i,j]<0){
      log10_fdr_sign[i, j] <- -log10_fdr_sign[i, j]
    }
  }
}
#log10_fdr_sign[log2_OR<0] <-log10_fdr_sign*-1 
col_fun <- colorRamp2(seq(0,70,length.out=9),
                      (brewer.pal(9, "Reds")))
rownames(log10_fdr_sign) <- tissue_info$tissue_abbrv
ht <- Heatmap(log10_fdr_sign,
              col = col_fun,
              na_col = "white",
              cluster_rows = F,
              cluster_columns = F,
              name = "-log10(FDR)",
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(round(log2_OR[i, j], 2), x, y, gp = gpar(fontsize = 10))},
              row_names_side = "left",
              column_names_side = "bottom",
              left_annotation = row_ha_left,
              row_names_gp = gpar(fontsize = 10),
              bottom_annotation = column_ha_bottom,
              column_names_gp = gpar(fontsize = 10),
              top_annotation = column_ha_top,
              heatmap_legend_param = list(direction = "horizontal"))


#######################################

# Build summary data.table --
# Tissue, pairwise:traits, overlap, odds.ratio, FDR 
data <- do.call(rbind.data.frame,
                lapply(tissues, function(tissue) 
                  data.frame("Tissue" = rep(tissue_info[tissue_info$tissue_ID==tissue, "tissue_abbrv"], length(pw.traits)),
                             "Tissue_ID" = rep(tissue,  length(pw.traits)),
                             "Traits" = pw.traits,
                             "Overlap" = overlap[tissue,], # all tissue and pairwise trait combinations
                             "OR" = or[tissue, ], # only fdr < 0.05
                             "fdr" = fdr[tissue,] # only fdr < 0.05
                  )))
data$FDR <- -log10(data$fdr)
data$FDR <- ifelse(data$OR < 1, (-1)*data$FDR, data$FDR) # if depletion, blue (negative)
data$OddsRatio <- log2(data$OR)
data$Overlap <- ifelse(data$Overlap=="0","",data$Overlap)
data$Overlap <- as.numeric(data$Overlap)
data$Tissue <- factor(data$Tissue, levels = rev(tissue_info$tissue_abbrv), order = T)
data <- data[!is.na(data$fdr),]
data$lower_CI <- sapply(1:nrow(data), function(row) 
  fisher.results[[ data[row, "Traits"] ]][[ data[row, "Tissue_ID"] ]][['lower_CI']])
data$upper_CI <- sapply(1:nrow(data), function(row) 
  fisher.results[[ data[row, "Traits"] ]][[ data[row, "Tissue_ID"] ]][['upper_CI']])
data <- data[data$OR>1 & data$Overlap >= 10,]
# Limit color scale (FDR) --
data$FDR_bounded <- ifelse(is.na(data$FDR), NA,
                           ifelse(data$FDR < -5, -5,
                                  ifelse(data$FDR > 5, 5, data$FDR))) # if depletion, blue (negative)
data$y_dummy <- paste0(data$Tissue, "_", data$Traits)

data$Traits <- gsub("Age:Ancestry", "Ancestry:Age", data$Traits)
data$Traits <- factor(data$Traits, levels = gsub("Age:Ancestry", "Ancestry:Age", pw.traits), order = T)
data <- data[order(data$Traits),]

data$y_dummy <- paste0(data$Tissue, "_", data$Traits)
data$y_dummy <- factor(data$y_dummy, levels = rev(paste0(data$Tissue, "_", data$Traits)), order = T)
#data <- data[data$Overlap>=10 & data$OR>1,]

#### enrichment overlap positions ancestry - Age colon transverse
shared_positions <- read.delim('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Overlap_Ancestry_Age_COLON.txt')
shared_positions <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Additive_Sex_Age_colon.rds')
shared_positions <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Additive_Sex_Age_lung.rds')
shared_positions <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Additive_Ancestry_Age_lung.rds')
shared_positions <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Additive_Ancestry_Age_ovary.rds')

annotation <- read.csv("~/marenostrum_scratch/MN4/bsc83/bsc83535/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")

tissues <- c("ColonTransverse")
names <- c("Age", "Sex")
traits_to_use <- c('SEX2','AGE')

results_DML <- lapply(tissues, function(tis)
  readRDS(paste0("Tissues/",tis,"/DML_results_5_PEERs_continous.rds")))
names(results_DML) <- tissues

files <- list.files('Data/EpiMap/', pattern='.bed.gz',full.names=T)

names_chrom <- c("ColonTransverse")
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

genes_age <- results_DML$Ovary$AGE[results_DML$Ovary$AGE$adj.P.Val<0.05,]
genes_sex <- results_DML$Ovary$EURv1[results_DML$Ovary$EURv1$adj.P.Val<0.05,]

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
  res <- shared_positions[shared_positions$class==1,]
  chrom_tissue <- as.data.frame(chromhmm_cpgs[[tissue]])
  chrom_tissue$region_chromhmm_new <- chrom_tissue$region_chromhmm
  chrom_tissue <- chrom_tissue[chrom_tissue$.overlap ==1,]
  #
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TssFlnkD", "TssFlnk", "TssFlnkU","TssA")] <- "TSS"
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("EnhA2", "EnhA1","EnhWk","EnhG1", "EnhG2")] <- "Enh"

  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("ReprPCWk","ReprPC")] <- "ReprPC"

  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TxWk","Tx")] <- "Tx"

  all_cpgs <- unique(c(unique(rownames(genes_sex)), unique(rownames(genes_age))))
  all_cpgs <- rownames(shared_positions)

  type_df <- chrom_tissue[chrom_tissue$region_chromhmm_new == type & chrom_tissue$name_ann %in% all_cpgs,]
  other_type <- chrom_tissue[chrom_tissue$region_chromhmm_new != type & chrom_tissue$name_ann %in% all_cpgs,]
  type_diff <- nrow(type_df[type_df$name_ann %in% rownames(res[res$logFC_Sex<0,]),])
  type_notdiff <- nrow(type_df) - type_diff
  other_type_diff <- nrow(other_type[other_type$name_ann %in% rownames(res[res$logFC_Sex<0,]),])
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
fisher_results <- lapply(families, function(region) my_fisher(region,'ColonTransverse'))
names(fisher_results) <- families


saveRDS(fisher_results, 'Tissues/enrichment_chromhmm_hyperMale_hyperYoung_CI.continous.colon.intersection.rds')

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
EUR_YOUNG <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hyperMale_hyperYoung_CI.continous.colon.intersection.rds')
FEM_OLD <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hyperFemale_hyperOLD_CI.continous.colon.intersection.rds')

library(ggpubr)

hypo_d <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), EUR_YOUNG)
hyper_d <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), FEM_OLD)
hyper_hypo <- rbind(hypo_d, hyper_d)
hyper_hypo$sig[hyper_hypo$sig =="Sig"] <- "FDR < 0.05"
hyper_hypo$sig[hyper_hypo$sig =="not Sig"] <- "FDR >= 0.05"
hyper_hypo$sig <- factor(hyper_hypo$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
hyper_hypo$region <- factor(hyper_hypo$region, levels=rev(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts")))
hyper_hypo$type[hyper_hypo$type =="EUR_YOUNG"] <- "Male - Young"
hyper_hypo$type[hyper_hypo$type =="FEM_OLD"] <- "Female - Old"
g <- ggplot(hyper_hypo, aes(x=log2(oddsRatio), y=region, colour=type, alpha=sig)) +
  geom_errorbar(aes(xmin=log2(CI_down), xmax=log2(CI_up)), width=.3) +
  geom_vline(xintercept = 0) +
  #xlim(0,20) + #Only for Lung to show the 0
  geom_point(size=3) + ylab('') + theme_bw() +
  scale_colour_manual(values=rev(c("#00b159", "#f37735"))) +
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
  scale_fill_manual(values=rev(c("#00b159", "#f37735"))) +
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
pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/enrichment_overlap_Age_Sex_colon.intersection.pdf"), w = 8, h = 4)
print(p)
dev.off()

