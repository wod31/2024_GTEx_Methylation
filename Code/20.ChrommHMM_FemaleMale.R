#### plot enrichments sex in female/male chromhmmm ####
#### Load results and data
annotation <- read.csv("marenostrum_scratch///GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")
first_dir <- "/gpfs/projects/bsc83/"
first_dir <- "marenostrum/"

project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("ColonTransverse")
names <- c("Sex")
traits_to_use <- c('SEX2')

results_DML <- lapply(tissues, function(tis) 
  readRDS(paste0("marenostrum/Projects/GTEx_v8/Methylation/Tissues/",tis,"/DML_results_5_PEERs_continous.rds")))
names(results_DML) <- tissues

files <- list.files('~/marenostrum/Projects/GTEx_v8/Methylation/Data/EpiMap/Sex-sp/', pattern='ColonTr.*.bed.gz',full.names=T)
#files <- files[2:5]

names_chrom <- c("Female.BSS01848",'Female.BSS01849','Male.BSS01850','Male.BSS01851')
chromhmm <- lapply(names_chrom, function(tis)
  read.delim(files[grep(tis, files)], sep='\t', header=F))
names(chromhmm) <- c('Female1','Female2','Male1','Male2')

### overlap chromhmm and cpgs tested #### 
ann_bed <- annotation[!is.na(annotation$MAPINFO) & !is.na(annotation$CHR),] %>%
  dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>% distinct()
ann_bed$chrom <- paste0('chr',ann_bed$chrom)
head(ann_bed)
ann_bed$start <- ann_bed$start-1

library(valr)
chromhmm_cpgs <- lapply(c('Female1','Female2','Male1','Male2'), function(tis) {
  chrom_df <- chromhmm[[tis]][,c(1:4)]
  colnames(chrom_df) <- c('chrom','start','end','region')
  bed_intersect(ann_bed, chrom_df, suffix = c("_ann", "_chromhmm"))})
names(chromhmm_cpgs) <- c('Female1','Female2','Male1','Male2')

#### enrichment -- chi-squared first ####
#From 18 states to 14

my_fisher <- function(type, tissue, trait){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  print(tissue)
  res <- results_DML[['ColonTransverse']][[trait]]
  chrom_tissue <- as.data.frame(chromhmm_cpgs[[tissue]])
  chrom_tissue <- chrom_tissue[chrom_tissue$.overlap ==1,]
  chrom_tissue$region_chromhmm_new <- chrom_tissue$region_chromhmm
  # 
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TssFlnkD", "TssFlnk", "TssFlnkU","TssA")] <- "TSS"
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("EnhA2", "EnhA1","EnhWk","EnhG1", "EnhG2")] <- "Enh"
  
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("ReprPCWk","ReprPC")] <- "ReprPC"
  
  chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TxWk","Tx")] <- "Tx"
  
  type_df <- chrom_tissue[chrom_tissue$region_chromhmm_new == type & chrom_tissue$name_ann %in% rownames(res),]
  other_type <- chrom_tissue[chrom_tissue$region_chromhmm_new != type & chrom_tissue$name_ann %in% rownames(res),]
  type_diff <- nrow(type_df[type_df$name_ann %in% rownames(res[res$adj.P.Val<0.05 & res$logFC>0,]),])
  type_notdiff <- nrow(type_df) - type_diff
  other_type_diff <- nrow(other_type[other_type$name_ann %in% rownames(res[res$adj.P.Val<0.05 & res$logFC>0,]),])
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
#families <- as.vector(unique(chromhmm_cpgs$Ovary$region_chromhmm))
families <- c('Enh','EnhBiv','Het','Quies','ReprPC','TSS','TssBiv','Tx','ZNF/Rpts')
fisher_results <- lapply(c('Female1','Female2','Male1','Male2'), function(tis) lapply(c('SEX2'), function(trait) lapply(families, function(region) my_fisher(region, tis,trait))))
names(fisher_results) <- c('Female1','Female2','Male1','Male2')

for (name in c('Female1','Female2','Male1','Male2')) {
  names(fisher_results[[name]]) <- 'Sex'
}

for (name in c('Female1','Female2','Male1','Male2')) {
  for (trait in c('Sex')) {
    names(fisher_results[[name]][[trait]]) <- families
  }
}

saveRDS(fisher_results, 'marenostrum/Projects//GTEx_v8/Methylation/Tissues/ColonTransverse///enrichment_chromhmm_hyper_batch_CI.simple.continous.Sex.rds')

## plot ###
read_data <- function(variables, data, tissue, trait){ #Function to prepare data to plot and compute adjusted p value
  
  odds_ratio <- lapply(variables, function(type) data[[tissue]][[trait]][[type]][['f']]$estimate)
  adj.P.Val <- p.adjust(sapply(variables, function(type) data[[tissue]][[trait]][[type]][['f']]$p.value), method = "BH")
  CI_down <- lapply(variables, function(type) data[[tissue]][[trait]][[type]][['f']]$conf.int[1])
  CI_up <- lapply(variables, function(type) data[[tissue]][[trait]][[type]][['f']]$conf.int[2])
  sample_size <- lapply(variables, function(type) data[[tissue]][[trait]][[type]][['m']])
  
  
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

hypo <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/ColonTransverse//enrichment_chromhmm_hypo_batch_CI.simple.continous.Sex.rds')
hyper <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/ColonTransverse//enrichment_chromhmm_hyper_batch_CI.simple.continous.Sex.rds')

colors_traits <- list('AGE'=c('#3D7CD0','#B4D6F6'),
                      'Sex'=c('#3B734E','#89AA94'),
                      'EURv1'=c('#F0AE21','#F9DE8B'))


for (tissue in  c('Female1','Female2','Male1','Male2')) {
  for (trait in c('Sex')) {
    
    hypo_d <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), hypo, tissue, trait)
    hyper_d <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), hyper, tissue, trait)
    hyper_hypo <- rbind(hypo_d, hyper_d)
    hyper_hypo$sig[hyper_hypo$sig =="Sig"] <- "FDR < 0.05"
    hyper_hypo$sig[hyper_hypo$sig =="not Sig"] <- "FDR >= 0.05"
    hyper_hypo$sig <- factor(hyper_hypo$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
    hyper_hypo$region <- factor(hyper_hypo$region, levels=rev(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts")))
    hyper_hypo$type[hyper_hypo$type =="hypo"] <- "Hypomethylation"
    hyper_hypo$type[hyper_hypo$type =="hyper"] <- "Hypermethylation"
    g <- ggplot(hyper_hypo, aes(x=(oddsRatio), y=region, colour=type, alpha=sig)) +
      geom_errorbar(aes(xmin=(CI_down), xmax=(CI_up)), width=.3) +
      geom_vline(xintercept = 1) +
      #xlim(0,8) + #Only for Lung to show the 0
      geom_point(size=3) + ylab('') + theme_bw() +
      scale_colour_manual(values=colors_traits[[trait]]) +
      xlab("Odds ratio") +
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
    # pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Ovary_menopause/genomic_location_", tissue,'_',trait,".pdf"), w = 6, h = 3.5)
    # print(g)
    # dev.off()
    
    #Plot sample sizes:
    
    g2 <- ggplot(hyper_hypo) + geom_col(aes(sample_size, region, fill=type), width = 0.6) +
      theme_classic() + xlab("Number of DMPs") + ylab("") +
      scale_fill_manual(values=colors_traits[[trait]]) +
      theme(legend.position = "none",
            axis.text.x = element_text(colour="black", size=13),
            axis.text.y = element_blank(),
            axis.title.x = element_text(size=16)) +
      scale_x_continuous(n.breaks=3)+
      scale_y_discrete(breaks=c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"),
                       labels=c("Enhancer","Enhancer Bivalent","Heterochromatin","Quiescent","Repressed Polycomb","TSS","TSS Bivalent","Transcription","ZNF & Repeats")) #+
    #scale_x_continuous(breaks=c(0, 20000, 40000)) #Only for lung
    p <- ggarrange(g, g2, labels = c("A", "B"),
                   common.legend = TRUE, legend = "right", widths = c(0.8,0.3))
    pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/genomic_location_Colon_", tissue,'_',trait,"_5peers.pdf"), w = 8, h = 4)
    print(p)
    dev.off()
  }
}

#### get specific regions ######
# Create catalogue 
catalog <- ann_bed
for (tissue in c('Female1','Female2','Male1','Male2')) {
  cat <- chromhmm_cpgs[[tissue]]
  cat <- cat[cat$.overlap ==1,]
  cat <- as.data.frame(cat)
  rownames(cat) <- cat$name_ann  
  cat$region_chromhmm[cat$region_chromhmm %in% c("TssFlnkD", "TssFlnk", "TssFlnkU","TssA")] <- "TSS"
  cat$region_chromhmm[cat$region_chromhmm %in% c("EnhA2", "EnhA1","EnhWk","EnhG1", "EnhG2")] <- "Enh"
  
  cat$region_chromhmm[cat$region_chromhmm %in% c("ReprPCWk","ReprPC")] <- "ReprPC"
  
  cat$region_chromhmm[cat$region_chromhmm %in% c("TxWk","Tx")] <- "Tx"
  catalog[,tissue] <- cat[catalog$name,'region_chromhmm']
}

catalog$Female <- catalog$Female1
catalog$Female[catalog$Female1 != catalog$Female2] <- 'NotConsistent'
catalog$Male <- catalog$Male1
catalog$Male[catalog$Male1 != catalog$Male2] <- 'NotConsistent'

catalog$DMP <- 'No'
catalog$DMP[catalog$name %in% rownames(results_DML$ColonTransverse$SEX2[results_DML$ColonTransverse$SEX2$adj.P.Val<0.05,])] <- 'DMP'
catalog_filt <- catalog[catalog$Female != 'NotConsistent',]
catalog_filt <- catalog_filt[catalog_filt$Male != 'NotConsistent',]
catalog_filt$region <- 'Same'
catalog_filt$region[catalog_filt$Female != catalog_filt$Male] <- 'Sex-specific'

results_t <- as.data.frame(table(catalog_filt$DMP, catalog_filt$region))

pdf('marenostrum/Projects/GTEx_v8/Methylation/Plots/n_sex_specific_colon_dmps.pdf', width = 4, height = 3)
ggplot(catalog_filt, aes(fill=region,x=DMP)) + 
  geom_bar(position="fill", stat="count")
dev.off()

### Enrichment sex-specific regions #####
my_fisher <- function(type, trait){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  print(trait)
  res <- results_DML[['ColonTransverse']][[trait]]
  dmps <- catalog_filt[catalog_filt$DMP == 'DMP',]
  
  type_df <- dmps[dmps$Male == type & dmps$region == 'Sex-specific',]
  other_type <- dmps[dmps$Male != type | dmps$region == 'Same',]
  type_diff <- nrow(type_df[type_df$name %in% rownames(res[res$adj.P.Val<0.05 & res$logFC>0,]),])
  type_notdiff <- nrow(type_df) - type_diff
  other_type_diff <- nrow(other_type[other_type$name %in% rownames(res[res$adj.P.Val<0.05 & res$logFC>0,]),])
  other_type_notdiff <- nrow(other_type) - other_type_diff
  
  ### test significance
  m <- matrix(c(type_diff, other_type_diff, type_notdiff, other_type_notdiff), 2,2, byrow = T)
  print(m)  
  
  m[is.na(m)] <- 0
  #m <- m[c(type,paste0('No ',type)),]
  rownames(m) <- c(type, "Other")
  colnames(m) <- c("Hyper","Hypo")
  print(m)
  f <- fisher.test(m)
  print(f)
  return(list("f" = f, "m" = type_diff))
}
### now vs all CpGs, not hyper vs hypo
my_fisher <- function(type, trait){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  print(trait)
  res <- results_DML[['ColonTransverse']][[trait]]
  dmps <- catalog_filt
  
  type_df <- dmps[dmps$Male == type & dmps$region == 'Sex-specific',]
  other_type <- dmps[dmps$Male != type | dmps$region == 'Same',]
  type_diff <- nrow(type_df[type_df$name %in% rownames(res[res$adj.P.Val<0.05 & res$logFC>0,]),])
  type_notdiff <- nrow(type_df) - type_diff
  other_type_diff <- nrow(other_type[other_type$name %in% rownames(res[res$adj.P.Val<0.05 & res$logFC>0,]),])
  other_type_notdiff <- nrow(other_type) - other_type_diff
  
  ### test significance
  m <- matrix(c(type_diff, other_type_diff, type_notdiff, other_type_notdiff), 2,2, byrow = T)
  print(m)  
  
  m[is.na(m)] <- 0
  #m <- m[c(type,paste0('No ',type)),]
  rownames(m) <- c(type, "Other")
  colnames(m) <- c("Hyper","Hypo")
  print(m)
  f <- fisher.test(m)
  print(f)
  return(list("f" = f, "m" = type_diff))
}
# Two-tailed Fisher test
#families <- as.vector(unique(chromhmm_cpgs$Ovary$region_chromhmm))
families <- c('Enh','EnhBiv','Het','Quies','ReprPC','TSS','TssBiv','Tx','ZNF/Rpts')
fisher_results <- lapply(families, function(region) lapply(c('SEX2'), function(trait) my_fisher(region,trait)))
names(fisher_results) <- families

for (name in families) {
  names(fisher_results[[name]]) <- 'Sex'
}

read_data <- function(variables, data, trait){ #Function to prepare data to plot and compute adjusted p value
  
  odds_ratio <- lapply(variables, function(type) data[[type]][[trait]][['f']]$estimate)
  adj.P.Val <- p.adjust(sapply(variables, function(type) data[[type]][[trait]][['f']]$p.value), method = "BH")
  CI_down <- lapply(variables, function(type) data[[type]][[trait]][['f']]$conf.int[1])
  CI_up <- lapply(variables, function(type) data[[type]][[trait]][['f']]$conf.int[2])
  sample_size <- lapply(variables, function(type) data[[type]][[trait]][['m']])
  
  
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
# 
# female_vs_male <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), fisher_results, 'Sex')
# male_vs_female <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), fisher_results, 'Sex')
# male_vs_female_malesp <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), fisher_results, 'Sex')
# female_vs_male_malesp <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), fisher_results, 'Sex')
# female_sp <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), fisher_results, 'Sex')
# female_sp_male <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), fisher_results, 'Sex')
# male_sp <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), fisher_results, 'Sex')
# male_sp_female <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), fisher_results, 'Sex')

#### merge data with genes ######
ann <- read.delim('marenostrum/Projects/GTEx_v8/Methylation/Data/Methylation_Epic_gene_promoter_enhancer_processed.txt')
catalog_filt <- merge(catalog_filt, ann, by.x='name',by.y='IlmnID')
hyper <- rownames(results_DML$ColonTransverse$SEX2[results_DML$ColonTransverse$SEX2$logFC>0 & results_DML$ColonTransverse$SEX2$adj.P.Val<0.05,])
hypo <- rownames(results_DML$ColonTransverse$SEX2[results_DML$ColonTransverse$SEX2$logFC<0 & results_DML$ColonTransverse$SEX2$adj.P.Val<0.05,])

catalog_filt$dir <- 'NA'
catalog_filt$dir[catalog_filt$name %in% hyper] <- 'Female-biased'
catalog_filt$dir[catalog_filt$name %in% hypo] <- 'Male-biased'

View(catalog_filt[catalog_filt$dir == 'Female-biased',])

## plot together in a heatmap
hyper_all <- lapply( c('Female1','Female2','Male1','Male2'), function(x) read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), hyper, x, 'Sex'))
names(hyper_all) <- c('Female1','Female2','Male1','Male2')

hypo_all <- lapply( c('Female1','Female2','Male1','Male2'), function(x) read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), hypo, x, 'Sex'))
names(hypo_all) <- c('Female1','Female2','Male1','Male2')

hyper_all <- do.call('rbind.data.frame', hyper_all)
hypo_all <- do.call('rbind.data.frame', hypo_all)

#all <- rbind(hyper_all, hypo)
hyper_all$chromhmm <- gsub('\\..*','',rownames(hyper_all))

pdf('marenostrum/Projects/GTEx_v8/Methylation/Plots/hypermethylation_enrichment.heatmap.pdf', width = 4.3, height = 3)
ggplot(hyper_all, aes(x=chromhmm, y=region, size=-log10(adjPvalue), fill=log2(oddsRatio))) +
  geom_point(alpha=0.8, shape=21, color="black") +
  scale_size(range = c(0.5, 10), name="-log10(FDR)", limits = c(0,max(-log10(hyper_all$adjPvalue)))) +
  scale_fill_gradient2(low="black", mid = "white" ,high="#3B734E", name='log2(oddsRatio)', midpoint = -0.5,)+
  #theme_ipsum() +
  theme(legend.position="right") +
  ylab("") +
  xlab("") +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        panel.grid = element_line(colour = 'light grey'),
        axis.text.x = element_text(colour="black", size=10),
        axis.text.y = element_text(colour="black", size=12),
        legend.text = element_text(colour="black", size=8),
        axis.title.x = element_text(size=12))
dev.off()



### plot sex enrichments by tissue

read_data <- function(variables, data, tissue, trait){ #Function to prepare data to plot and compute adjusted p value
  
  odds_ratio <- lapply(variables, function(type) data[[tissue]][[trait]][[type]][['f']]$estimate)
  adj.P.Val <- p.adjust(sapply(variables, function(type) data[[tissue]][[trait]][[type]][['f']]$p.value), method = "BH")
  CI_down <- lapply(variables, function(type) data[[tissue]][[trait]][[type]][['f']]$conf.int[1])
  CI_up <- lapply(variables, function(type) data[[tissue]][[trait]][[type]][['f']]$conf.int[2])
  sample_size <- lapply(variables, function(type) data[[tissue]][[trait]][[type]][['m']])
  
  
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

hypo <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hypo_batch_CI.simple.continous.rds')
hyper <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hyper_batch_CI.simple.continous.rds')

colors_traits <- list('AGE'=c('#3D7CD0','#B4D6F6'),
                      'SEX2'=c('#3B734E','#89AA94'),
                      'EURv1'=c('#F0AE21','#F9DE8B'))

tissue_info <- readRDS(paste0('~/marenostrum/', "Projects/GTEx_v8/Jose/00_Data/Tissue_info_whole.rds"))

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "KidneyCortex", "Testis", "WholeBlood","MuscleSkeletal")

sex_tissues <- c('Ovary','Prostate','Testis')

library(ggplot2)
library(ggpubr)

for (tissue in names(hypo)) {
  for (trait in c('SEX2')) {
    if (tissue %in% sex_tissues & trait == "SEX2") {
      print(NA)
    } else {
    
    hypo_d <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), hypo, tissue, trait)
    hyper_d <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), hyper, tissue, trait)
    hyper_hypo <- rbind(hypo_d, hyper_d)
    hyper_hypo$sig[hyper_hypo$sig =="Sig"] <- "FDR < 0.05"
    hyper_hypo$sig[hyper_hypo$sig =="not Sig"] <- "FDR >= 0.05"
    hyper_hypo$sig <- factor(hyper_hypo$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
    hyper_hypo$region <- factor(hyper_hypo$region, levels=rev(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts")))
    hyper_hypo$type[hyper_hypo$type =="hypo"] <- "Hypomethylation"
    hyper_hypo$type[hyper_hypo$type =="hyper"] <- "Hypermethylation"
    g <- ggplot(hyper_hypo, aes(x=log2(oddsRatio), y=region, colour=type, alpha=sig)) +
      geom_errorbar(aes(xmin=log2(CI_down), xmax=log2(CI_up)), width=.3) +
      geom_vline(xintercept = 0) +
      #xlim(0,8) + #Only for Lung to show the 0
      geom_point(size=3) + ylab('') + theme_bw() +
      scale_colour_manual(values=colors_traits[[trait]]) +
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
    # pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Ovary_menopause/genomic_location_", tissue,'_',trait,".pdf"), w = 6, h = 3.5)
    # print(g)
    # dev.off()
    
    #Plot sample sizes:
    
    g2 <- ggplot(hyper_hypo) + geom_col(aes(sample_size, region, fill=type), width = 0.6) +
      theme_classic() + xlab("Number of DMPs") + ylab("") +
      scale_fill_manual(values=colors_traits[[trait]]) +
      theme(legend.position = "none",
            axis.text.x = element_text(colour="black", size=13),
            axis.text.y = element_blank(),
            axis.title.x = element_text(size=16)) +
      scale_x_continuous(n.breaks=3)
      # scale_y_discrete(breaks=c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"),
      #                  labels=c("Enhancer","Enhancer Bivalent","Heterochromatin","Quiescent","Repressed Polycomb","TSS","TSS Bivalent","Transcription","ZNF & Repeats")) #+
    #scale_x_continuous(breaks=c(0, 20000, 40000)) #Only for lung
    p <- ggarrange(g, g2, labels = c("A", "B"),
                   common.legend = TRUE, legend = "right", widths = c(0.8,0.3))
    pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/genomic_location_", tissue,'_',trait,"_5peers.new.pdf"), w = 8, h = 4)
    print(p)
    dev.off()
    }
  }
}


#### number of peaks per sample ####
for (ind in c('Female1','Female2','Male1','Male2')) {
  print(ind)
  print(table(chromhmm[[ind]]$V4))
}

for (ind in c('Female1','Female2','Male1','Male2')) {
  print(ind)
  chromhmm[[ind]]$length <- chromhmm[[ind]]$V3-chromhmm[[ind]]$V2
  print(sum(chromhmm[[ind]]$length[chromhmm[[ind]]$V4 %in% c('ReprPC','ReprPCWk')]))
}
