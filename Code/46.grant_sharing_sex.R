#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Analyze sharing sex-DMPs with Grant et al., results and perform enrichments with chromhmm
# @software version: R=4.2.2


#### sharing grant et al. #####
sharing <- readRDS('Data/Sharing_DMP.rds')
shared_cpgs <- readRDS('Data/chromHMM_shared_9tissues.rds')
blood <- read.delim("Data/papers/Grant_DMPs_autosomal_13148_2022_1279_MOESM1_ESM.csv", sep=',')

#blood <- blood[blood$probe%in% background,] #Only 4160 out of 76956 (%5.4) of the loci are in the EPIC

blood_up <- blood[blood$logFC<0,] # EUR
blood_down <- blood[blood$logFC>0,] # AFR

for (cpg in blood_up$Row.names) {
  if (cpg %in% sharing$CG[sharing$trait=='SEX2' & sharing$dir==1]) {
    sharing$number[sharing$CG == cpg & sharing$trait=='SEX2' & sharing$dir==1] <- sharing$number[sharing$CG == cpg & sharing$trait=='SEX2' & sharing$dir==1] + 1
  } else {
    print('not shared')
  }
}

for (cpg in blood_down$Row.names) {
  if (cpg %in% sharing$CG[sharing$trait=='SEX2' & sharing$dir==-1]) {
    sharing$number[sharing$CG == cpg & sharing$trait=='SEX2' & sharing$dir==-1] <- sharing$number[sharing$CG == cpg & sharing$trait=='SEX2' & sharing$dir==-1] + 1
  } else {
    print('not shared')
  }
}

my_fisher <- function(type, trait){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  res <- sharing[sharing$trait == trait,]
  chrom_tissue <- shared_cpgs
  
  type_df <- chrom_tissue[chrom_tissue$region_chromhmm == type & chrom_tissue$name_ann %in% res$CG,]
  other_type <- chrom_tissue[chrom_tissue$region_chromhmm != type & chrom_tissue$name_ann %in% res$CG,]
  type_diff <- nrow(type_df[type_df$name_ann %in% res$CG[res$number>=2 & res$dir==-1],])
  type_notdiff <- nrow(type_df) - type_diff
  other_type_diff <- nrow(other_type[other_type$name_ann %in% res$CG[res$number>=2 & res$dir==-1],])
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
fisher_results <- lapply(c("SEX2" ), function(trait) lapply(families, function(region) my_fisher(region,trait)))
names(fisher_results) <- c("Sex")

for (name in c("Sex")) {
  names(fisher_results[[name]]) <- families
}

saveRDS(fisher_results, 'Tissues/enrichment_chromhmm_hyper_shared_sex.grant.rds')
saveRDS(fisher_results, 'Tissues/enrichment_chromhmm_hypo_shared_sex.grant.rds')

read_data <- function(variables, data, trait){ #Function to prepare data to plot and compute adjusted p value
  
  odds_ratio <- lapply(variables, function(type) data[[trait]][[type]][['f']]$estimate)
  adj.P.Val <- p.adjust(sapply(variables, function(type) data[[trait]][[type]][['f']]$p.value), method = "BH")
  CI_down <- lapply(variables, function(type) data[[trait]][[type]][['f']]$conf.int[1])
  CI_up <- lapply(variables, function(type) data[[trait]][[type]][['f']]$conf.int[2])
  sample_size <- lapply(variables, function(type) data[[trait]][[type]][['m']])
  
  
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


hypo <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hypo_shared_sex.grant.rds')
hyper <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hyper_shared_sex.grant.rds')


colors_traits <- list('Age'=c('#3D7CD0','#B4D6F6'),
                      'Sex'=c('#3B734E','#89AA94'),
                      'Ancestry'=c('#F0AE21','#F9DE8B'))

library(ggpubr)
for (trait in c('Sex')) {
  
  hypo_d <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), hypo, trait)
  hyper_d <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), hyper, trait)
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
    #xlim(0,20) + #Only for Lung to show the 0
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
  # pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/genomic_location_shared",'_',trait,".pdf"), w = 6, h = 3.5)
  # print(g)
  # dev.off()
  
  #Plot sample sizes:
  
  g2 <- ggplot(hyper_hypo) + geom_col(aes(sample_size, region, fill=type), width = 0.6) +
    theme_classic() + xlab("Number of DMPs") + ylab("") +
    scale_fill_manual(values=colors_traits[[trait]]) +
    theme(legend.position = "none",
          axis.text.x = element_text(colour="black", size=13),
          axis.text.y=element_blank(),
          axis.title.x = element_text(size=16)) +
    scale_y_discrete(breaks=c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"),
                     labels=c("Enhancer","Enhancer Bivalent","Heterochromatin","Quiescent","Repressed Polycomb","TSS","TSS Bivalent","Transcription","ZNF & Repeats")) #+
  #scale_x_continuous(breaks=c(0, 20000, 40000)) #Only for lung
  # pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/genomic_location_shared",'_',trait,"_sample_size.pdf"), w = 4, h = 3.5)
  # print(g2)
  # dev.off()
  
  p <- ggarrange(g, g2, labels = c("A", "B"),
                 common.legend = TRUE, legend = "right", widths = c(0.8,0.3))
  pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/enrichment",'_',trait,"_shared.Grant.2.pdf"), w = 8, h = 4)
  print(p)
  dev.off()
}


