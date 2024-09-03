#### Analyze enrichment extreme cpgs ####
### plot enrichments sex in female/male chromhmmm ####
#### Load results and data
annotation <- read.csv("~/marenostrum_scratch/MN4/bsc83/bsc83535/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")
#first_dir <- "/gpfs/projects/bsc83/"
first_dir <- "marenostrum/"

project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Age", "Ancestry", "BMI", "Sex")

tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Methylation/Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

traits_to_use <- c('EURv1')

results_DML <- lapply(tissues, function(tis) 
  readRDS(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/",tis,"/DML_results_5_PEERs_continous.rds")))
names(results_DML) <- tissues

extremes <- readRDS(paste0(project_path,'Data/extreme_cpgs_09_01.rds'))

#### enrichment -- chi-squared first ####

my_fisher <- function(type, tissue, trait){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  print(tissue)
  res <- results_DML[[tissue]][[trait]]
  chrom_tissue <- extremes[[type]][[tissue]]
  
  type_diff <- sum(rownames(res[res$adj.P.Val<0.05 & res$logFC<0,]) %in% chrom_tissue)
  type_notdiff <- length(chrom_tissue) - type_diff
  other_type_diff <- (nrow(res[res$adj.P.Val<0.05 & res$logFC<0,]))-type_diff
  other_type_notdiff <- nrow(res) - (nrow(res[res$adj.P.Val<0.05 & res$logFC<0,])) - length(chrom_tissue)
  
  ### test significance
  m <- matrix(c(type_diff, type_notdiff, other_type_diff, other_type_notdiff), 2,2, byrow = T)
  print(m)  
  
  m[is.na(m)] <- 0
  #m <- m[c(type,paste0('No ',type)),]
  rownames(m) <- c(type, "Other")
  colnames(m) <- c("DMP","Not DMP")
  print(m)
  f <- fisher.test(m)
  print(f)
  return(list("f" = f, "m" = type_diff))
}
# Two-tailed Fisher test
#families <- as.vector(unique(chromhmm_cpgs$Ovary$region_chromhmm))
families <- c('hyper','hypo')
fisher_results <- lapply(tissues, function(tis) lapply(c('EURv1'), function(trait) lapply(families, function(region) my_fisher(region, tis,trait))))
names(fisher_results) <- tissues

for (name in tissues) {
  names(fisher_results[[name]]) <- 'EURv1'
}

for (name in tissues) {
  for (trait in c('EURv1')) {
    names(fisher_results[[name]][[trait]]) <- families
  }
}

saveRDS(fisher_results, '~/marenostrum/Projects//GTEx_v8/Methylation/Data/enrichment_extremes_AA.eurv1.rds')

## plot ###
read_data <- function(variables, data, type, trait){ #Function to prepare data to plot and compute adjusted p value
  
  odds_ratio <- lapply(variables, function(tissue) data[[tissue]][[trait]][[type]][['f']]$estimate)
  adj.P.Val <- p.adjust(sapply(variables, function(tissue) data[[tissue]][[trait]][[type]][['f']]$p.value), method = "BH")
  CI_down <- lapply(variables, function(tissue) data[[tissue]][[trait]][[type]][['f']]$conf.int[1])
  CI_up <- lapply(variables, function(tissue) data[[tissue]][[trait]][[type]][['f']]$conf.int[2])
  sample_size <- lapply(variables, function(tissue) data[[tissue]][[trait]][[type]][['m']])
  
  
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

hyper <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/enrichment_extremes_EA.eurv1.rds')
hypo <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/enrichment_extremes_AA.eurv1.rds')

colors_traits <- list('AGE'=c('#3D7CD0','#B4D6F6'),
                      'Sex'=c('#3B734E','#89AA94'),
                      'EURv1'=c('#F0AE21','#F9DE8B'))

library(ggplot2)
library(ggpubr)
for (type in c('hyper','hypo')) {
  for (trait in c('EURv1')) {
    
    hypo_d <- read_data(tissues, hypo, type, trait)
    hyper_d <- read_data(tissues, hyper, type, trait)
    hyper_hypo <- rbind(hypo_d, hyper_d)
    hyper_hypo$sig[hyper_hypo$sig =="Sig"] <- "FDR < 0.05"
    hyper_hypo$sig[hyper_hypo$sig =="not Sig"] <- "FDR >= 0.05"
    hyper_hypo$sig <- factor(hyper_hypo$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
    hyper_hypo$tissue <- factor(hyper_hypo$tissue, levels=tissues)
    hyper_hypo$type[hyper_hypo$type =="hypo"] <- "AA"
    hyper_hypo$type[hyper_hypo$type =="hyper"] <- "EA"
    g <- ggplot(hyper_hypo, aes(x=log2(oddsRatio), y=tissue, colour=type, alpha=sig)) +
      geom_errorbar(aes(xmin=log2(CI_down), xmax=log2(CI_up)), width=.3) +
      geom_vline(xintercept = 0) + ggtitle(paste0(type, 'positions in genome'))+
      #xlim(0,20) + #Only for Lung to show the 0
      geom_point(size=3) + ylab('') + theme_bw() +
      scale_colour_manual(values=rev(colors_traits[[trait]])) +
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
      scale_fill_manual(values=rev(colors_traits[[trait]])) +
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
    pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/enrichment", gsub('\\/','_',type),'_',trait,".v2.pdf"), w = 8, h = 4)
    print(p)
    dev.off()
  }
}

