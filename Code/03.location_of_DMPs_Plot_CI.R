#### plot enrichments CI 
library(tidyr)
# 
# hypo <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hypo_batch_CI.simple.continous.rds')
# 
# # Plot
# odds_ratio <- list()
# adj.P.Val <- list()
# CI_down <- list()
# CI_up <- list()
# # odds ratio
# odds_ratio <- lapply(names(hypo), function(tis) lapply(names(hypo[[tis]]), function(trait)
#   sapply(names(hypo[[tis]][[trait]]), function(chr) hypo[[tis]][[trait]][[chr]][['f']]$estimate)))
# # adj.P.Val
# adj.P.Val <- lapply(names(hypo), function(tis) lapply(names(hypo[[tis]]), function(trait)
#   p.adjust(sapply(names(hypo[[tis]][[trait]]), function(chr) hypo[[tis]][[trait]][[chr]][['f']]$p.value), method = "BH")))
# # CI down
# CI_down <- lapply(names(hypo), function(tis) lapply(names(hypo[[tis]]), function(trait)
#   sapply(names(hypo[[tis]][[trait]]), function(chr) hypo[[tis]][[trait]][[chr]][['f']]$conf.int[1])))
# # CI up
# CI_up <- lapply(names(hypo), function(tis) lapply(names(hypo[[tis]]), function(trait)
#   sapply(names(hypo[[tis]][[trait]]), function(chr) hypo[[tis]][[trait]][[chr]][['f']]$conf.int[2])))
# 
# names(odds_ratio) <- names(hypo)
# names(adj.P.Val) <- names(hypo)
# names(CI_down) <- names(hypo)
# names(CI_up) <- names(hypo)
# #
# for (name in names(hypo)) {
#   names(odds_ratio[[name]]) <- names(hypo[[name]])
# }
# for (name in names(hypo)) {
#   names(adj.P.Val[[name]]) <- names(hypo[[name]])
# }
# for (name in names(hypo)) {
#   names(CI_down[[name]]) <- names(hypo[[name]])
# }
# for (name in names(hypo)) {
#   names(CI_up[[name]]) <- names(hypo[[name]])
# }
# #
# odds_ratio_df <- as.data.frame(unlist(odds_ratio))
# odds_ratio_df$label <- rownames(odds_ratio_df)
# odds_ratio_df <- odds_ratio_df %>% separate(col=label, into=c('tissue', 'trait', 'region','type'), sep='\\.')
# colnames(odds_ratio_df) <- c('oddsRatio','tissue', 'trait', 'region','type')
# 
# adj.P.Val_df <- as.data.frame(unlist(adj.P.Val))
# adj.P.Val_df$label <- rownames(adj.P.Val_df)
# adj.P.Val_df <- adj.P.Val_df %>% separate(col=label, into=c('tissue', 'trait', 'region','type'), sep='\\.')
# adj.P.Val_df$type <- 'adj.P.value'
# colnames(adj.P.Val_df) <- c('adjPvalue','tissue', 'trait', 'region','type')
# 
# CI_down_df <- as.data.frame(unlist(CI_down))
# CI_down_df$label <- rownames(CI_down_df)
# CI_down_df <- CI_down_df %>% separate(col=label, into=c('tissue', 'trait', 'region','type'), sep='\\.')
# CI_down_df$type <- 'CI_down'
# colnames(CI_down_df) <- c('CI_down','tissue', 'trait', 'region','type')
# 
# CI_up_df <- as.data.frame(unlist(CI_up))
# CI_up_df$label <- rownames(CI_up_df)
# CI_up_df <- CI_up_df %>% separate(col=label, into=c('tissue', 'trait', 'region','type'), sep='\\.')
# CI_up_df$type <- 'CI_up'
# colnames(CI_up_df) <- c('CI_up','tissue', 'trait', 'region','type')
# 
# all <- merge(odds_ratio_df, merge(adj.P.Val_df, merge(CI_down_df, CI_up_df, by=c('tissue','trait','region')),  by=c('tissue','trait','region')),  by=c('tissue','trait','region'))
# head(all)
# all$sig <- 'not Sig'
# all$sig[all$adjPvalue<0.05] <- 'Sig'
# 
# ### hyper
# hypo <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hyper_batch_CI.simple.continous.rds')
# 
# # Plot
# odds_ratio <- list()
# adj.P.Val <- list()
# CI_down <- list()
# CI_up <- list()
# # odds ratio
# odds_ratio <- lapply(names(hypo), function(tis) lapply(names(hypo[[tis]]), function(trait)
#   sapply(names(hypo[[tis]][[trait]]), function(chr) hypo[[tis]][[trait]][[chr]][['f']]$estimate)))
# # adj.P.Val
# adj.P.Val <- lapply(names(hypo), function(tis) lapply(names(hypo[[tis]]), function(trait)
#   p.adjust(sapply(names(hypo[[tis]][[trait]]), function(chr) hypo[[tis]][[trait]][[chr]][['f']]$p.value), method = "BH")))
# # CI down
# CI_down <- lapply(names(hypo), function(tis) lapply(names(hypo[[tis]]), function(trait)
#   sapply(names(hypo[[tis]][[trait]]), function(chr) hypo[[tis]][[trait]][[chr]][['f']]$conf.int[1])))
# # CI up
# CI_up <- lapply(names(hypo), function(tis) lapply(names(hypo[[tis]]), function(trait)
#   sapply(names(hypo[[tis]][[trait]]), function(chr) hypo[[tis]][[trait]][[chr]][['f']]$conf.int[2])))
# 
# names(odds_ratio) <- names(hypo)
# names(adj.P.Val) <- names(hypo)
# names(CI_down) <- names(hypo)
# names(CI_up) <- names(hypo)
# #
# for (name in names(hypo)) {
#   names(odds_ratio[[name]]) <- names(hypo[[name]])
# }
# for (name in names(hypo)) {
#   names(adj.P.Val[[name]]) <- names(hypo[[name]])
# }
# for (name in names(hypo)) {
#   names(CI_down[[name]]) <- names(hypo[[name]])
# }
# for (name in names(hypo)) {
#   names(CI_up[[name]]) <- names(hypo[[name]])
# }
# #
# odds_ratio_df <- as.data.frame(unlist(odds_ratio))
# odds_ratio_df$label <- rownames(odds_ratio_df)
# odds_ratio_df <- odds_ratio_df %>% separate(col=label, into=c('tissue', 'trait', 'region','type'), sep='\\.')
# colnames(odds_ratio_df) <- c('oddsRatio','tissue', 'trait', 'region','type')
# 
# adj.P.Val_df <- as.data.frame(unlist(adj.P.Val))
# adj.P.Val_df$label <- rownames(adj.P.Val_df)
# adj.P.Val_df <- adj.P.Val_df %>% separate(col=label, into=c('tissue', 'trait', 'region','type'), sep='\\.')
# adj.P.Val_df$type <- 'adj.P.value'
# colnames(adj.P.Val_df) <- c('adjPvalue','tissue', 'trait', 'region','type')
# 
# CI_down_df <- as.data.frame(unlist(CI_down))
# CI_down_df$label <- rownames(CI_down_df)
# CI_down_df <- CI_down_df %>% separate(col=label, into=c('tissue', 'trait', 'region','type'), sep='\\.')
# CI_down_df$type <- 'CI_down'
# colnames(CI_down_df) <- c('CI_down','tissue', 'trait', 'region','type')
# 
# CI_up_df <- as.data.frame(unlist(CI_up))
# CI_up_df$label <- rownames(CI_up_df)
# CI_up_df <- CI_up_df %>% separate(col=label, into=c('tissue', 'trait', 'region','type'), sep='\\.')
# CI_up_df$type <- 'CI_up'
# colnames(CI_up_df) <- c('CI_up','tissue', 'trait', 'region','type')
# 
# all_hyper <- merge(odds_ratio_df, merge(adj.P.Val_df, merge(CI_down_df, CI_up_df, by=c('tissue','trait','region')),  by=c('tissue','trait','region')),  by=c('tissue','trait','region'))
# head(all)
# all_hyper$sig <- 'not Sig'
# all_hyper$sig[all_hyper$adjPvalue<0.05] <- 'Sig'
# 
# ### all hyper and hypo
# all$type <- 'hypo'
# all_hyper$type <- 'hyper'
# hyper_hypo <- rbind(all, all_hyper)
# head(hyper_hypo)
#saveRDS(all_info_hypo, '/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm.rds')

# library(ggh4x)
# tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Jose/00_Data/Tissue_info_whole.rds"))
# tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
# 
# tissues_cols <- tissue_info[unique(hyper_hypo$tissue[hyper_hypo$adjPvalue<0.05 & hyper_hypo$trait=='SEX2']), 3]
# names(tissues_cols) <- tissue_info[unique(hyper_hypo$tissue[hyper_hypo$adjPvalue<0.05 & hyper_hypo$trait=='SEX2']), 1]
# 
# strip <- strip_themed(background_x = elem_list_rect(fill = tissues_cols))
# 
# hyper_hypo <- hyper_hypo[,c("tissue","trait","region","oddsRatio","adjPvalue","CI_down","CI_up","sig","type")]
# subset <- hyper_hypo[hyper_hypo$adjPvalue<0.05 & hyper_hypo$trait=='SEX2',]
# subset$tissue <- factor(subset$tissue, levels = tissue_info[unique(hyper_hypo$tissue[hyper_hypo$adjPvalue<0.05 & hyper_hypo$trait=='SEX2']), 1])
# subset$oddsRatio[subset$oddsRatio == Inf] <- max(subset$oddsRatio[subset$oddsRatio != Inf])
# subset$CI_up[subset$CI_up == Inf] <- max(subset$CI_up[subset$CI_up != Inf])
#   
# pd <- position_dodge(0.1)
# ggplot(subset, aes(x=(oddsRatio), y=region, colour=type)) + 
#   geom_errorbar(aes(xmin=(CI_down), xmax=(CI_up)), width=.1, position=pd) + #xlim(c(0,10))+
#   geom_vline(xintercept = 1)+
#   geom_line(position=pd) +
#   geom_point(position=pd) + ylab('')+
#   facet_wrap2(~ tissue, strip = strip) + theme_bw()

#chromhmm_enrichment_simple_onlypeers_ancestry

####### Jose's code to plot both number and enrichment one plot per tissue and traitt

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

ovary_hyper <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/Ovary.filt.enrichment_chromhmm_hyper.continous.rds')
ovary_hypo <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/Ovary.filt.enrichment_chromhmm_hypo.continous.rds')

colon_hyper <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/Colon.peer.enrichment_chromhmm_hyper.continous.rds')
colon_hypo <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/Colon.peers.enrichment_chromhmm_hypo.continous.rds')

hyper <- list()
hypo <- list()

hyper[['ColonTransverse']][['AGE']] <- colon_hyper
hypo[['ColonTransverse']][['AGE']] <- colon_hypo

colors_traits <- list('AGE'=c('#3D7CD0','#B4D6F6'),
                      'SEX2'=c('#3B734E','#89AA94'),
                      'EURv1'=c('#F0AE21','#F9DE8B'))

sex_tissues <- c('Ovary','Prostate','Testis')
#for (tissue in names(hypo)) {
for (tissue in c('ColonTransverse')) {
  #for (trait in names(hypo$Lung)) {
  for (trait in c('AGE')) {
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
      # pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/genomic_location_", tissue,'_',trait,".pdf"), w = 6, h = 3.5)
      # print(g)
      # dev.off()
      
      #Plot sample sizes:
      
      g2 <- ggplot(hyper_hypo) + geom_col(aes(sample_size, region, fill=type), width = 0.6) +
        theme_classic() + xlab("Number of DMPs") + ylab("") +
        scale_fill_manual(values=colors_traits[[trait]]) +
        theme(legend.position = "none",
              axis.text.x = element_text(colour="black", size=13),
              axis.text.y=element_blank(),  #remove y axis labels,
              axis.title.x = element_text(size=16)) +
        scale_x_continuous(n.breaks=3)
      #scale_x_continuous(breaks=c(0, 20000, 40000)) #Only for lung
      # pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/genomic_location_", tissue,'_',trait,"_sample_size.pdf"), w = 4, h = 3.5)
      # print(g2)
      # dev.off()
      
      p <- ggarrange(g, g2, labels = c("A", "B"),
                     common.legend = TRUE, legend = "right", widths = c(0.8,0.3))
      pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/enrichment_", tissue,'_',trait,".v2.filt.pdf"), w = 8, h = 4)
      print(p)
      dev.off()
    }
  }
}


####### Jose's code to plot both number and enrichment one plot per tissue and traitt
####### NOW plot per state and not per Tissue #######

read_data <- function(variables, data, type, trait){ #Function to prepare data to plot and compute adjusted p value
  
  # if (tissue %in% sex_tissues & trait == "SEX2") {
  #   print(NA)
  # } else {
  
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


hypo <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hypo_batch_CI.simple.continous.rds')
hyper <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hyper_batch_CI.simple.continous.rds')

colors_traits <- list('AGE'=c('#3D7CD0','#B4D6F6'),
                      'SEX2'=c('#3B734E','#89AA94'),
                      'EURv1'=c('#F0AE21','#F9DE8B'))

tissue_info <- readRDS(paste0('~/marenostrum/', "Projects/GTEx_v8/Jose/00_Data/Tissue_info_whole.rds"))

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "KidneyCortex", "Testis", "WholeBlood","MuscleSkeletal")

library(ggpubr)
sex_tissues <- c('Ovary','Prostate','Testis')
for (type in c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts")) {
  for (trait in c('AGE','SEX2','EURv1')) {
    if (trait == 'SEX2') {
      tissues_to_test <- tissues[!tissues %in% sex_tissues]
    } 
    if (trait == 'AGE') {
      tissues_to_test <- tissues[tissues != 'MuscleSkeletal']
    }
    if (trait == 'EURv1') {
      tissues_to_test <- tissues
    }
      hypo_d <- read_data(tissues_to_test, hypo, type, trait)
      hyper_d <- read_data(tissues_to_test, hyper, type, trait)
      hyper_hypo <- rbind(hypo_d, hyper_d)
      hyper_hypo$sig[hyper_hypo$sig =="Sig"] <- "FDR < 0.05"
      hyper_hypo$sig[hyper_hypo$sig =="not Sig"] <- "FDR >= 0.05"
      hyper_hypo$sig <- factor(hyper_hypo$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
      hyper_hypo$tissue <- factor(hyper_hypo$tissue, levels=tissues_to_test)
      hyper_hypo$type[hyper_hypo$type =="hypo"] <- "Hypomethylation"
      hyper_hypo$type[hyper_hypo$type =="hyper"] <- "Hypermethylation"
      g <- ggplot(hyper_hypo, aes(x=log2(oddsRatio), y=tissue, colour=type, alpha=sig)) +
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
              panel.border = element_rect(colour = "black", linewidth=1)) #+
        # scale_y_discrete(breaks=c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"),
        #                  labels=c("Enhancer","Enhancer Bivalent","Heterochromatin","Quiescent","Repressed Polycomb","TSS","TSS Bivalent","Transcription","ZNF & Repeats"))# + xlim(0, 3)
      # pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/genomic_location_", gsub('\\/','_',type),'_',trait,".pdf"), w = 6, h = 3.5)
      # print(g)
      # dev.off()
      
      #Plot sample sizes:
      
      g2 <- ggplot(hyper_hypo) + geom_col(aes(sample_size, tissue, fill=type), width = 0.6) +
        theme_classic() + xlab("Number of DMPs") + ylab("") +
        scale_fill_manual(values=colors_traits[[trait]]) +
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
      pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/enrichment", gsub('\\/','_',type),'_',trait,".v2.pdf"), w = 8, h = 4)
      print(p)
      dev.off()
    }
}

### try simpler figure for main v1 #####
tissue_info <- readRDS(paste0('~/marenostrum/', "Projects/GTEx_v8/Jose/00_Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

tissues_cols <- tissue_info[, 3]
names(tissues_cols) <- tissue_info[, 1]
tissues_cols <- tissues_cols[tissues]

traits_cols <- c('#C49122','#4B8C61','#70A0DF')

library(ggpubr)
sex_tissues <- c('Ovary','Prostate','Testis')
for (type in c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts")) {
  for (trait in names(hypo$Lung)) {
    if (trait == 'SEX2') {
      tissues_to_test <- tissues[!tissues %in% sex_tissues]
    } else {
      tissues_to_test <- tissues
    }
    hypo_d <- read_data(tissues_to_test, hypo, type, trait)
    hyper_d <- read_data(tissues_to_test, hyper, type, trait)
    hyper_hypo <- rbind(hypo_d, hyper_d)
    hyper_hypo$sig[hyper_hypo$sig =="Sig"] <- "FDR < 0.05"
    hyper_hypo$sig[hyper_hypo$sig =="not Sig"] <- "FDR >= 0.05"
    hyper_hypo$sig <- factor(hyper_hypo$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
    hyper_hypo$tissue <- factor(hyper_hypo$tissue, levels=tissues_to_test)
    hyper_hypo$type[hyper_hypo$type =="hypo"] <- "Hypo"
    hyper_hypo$type[hyper_hypo$type =="hyper"] <- "Hyper"
    g <- ggplot(hyper_hypo[hyper_hypo$sample_size>0,], aes(y=log2(oddsRatio), x=type)) +
      geom_boxplot(outlier.shape = NA) +
      geom_hline(yintercept = 0) +
      geom_jitter(aes(color=tissue, shape=sig), size=3)+
      #xlim(0,20) + #Only for Lung to show the 0
      theme_bw() + xlab('')+
      scale_colour_manual(values = tissues_cols) +
      ylab("log2(Odds ratio)") +
      #scale_alpha_discrete(range = c(0.2, 1), drop = FALSE) +
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
    
    g2 <- ggplot(hyper_hypo[hyper_hypo$sample_size>0,]) + geom_col(aes(y=sample_size, x=type, fill=tissue), width = 0.6) +
      theme_classic() + ylab("Number of DMPs") + xlab("") +
      scale_fill_manual(values = tissues_cols) +
      theme(legend.position = "none",
            axis.text.x = element_text(colour="black", size=13),
            axis.text.y = element_text(colour="black", size=14),
            #axis.text.y=element_blank(),  #remove y axis labels,
            axis.title.x = element_text(size=16)) +
      scale_y_continuous(n.breaks=6)
    #   scale_y_discrete(breaks=c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"),
    #                    labels=c("Enhancer","Enhancer Bivalent","Heterochromatin","Quiescent","Repressed Polycomb","TSS","TSS Bivalent","Transcription","ZNF & Repeats")) #+
    # #scale_x_continuous(breaks=c(0, 20000, 40000)) #Only for lung
    # pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/genomic_location_", gsub('\\/','_',type),'_',trait,"_sample_size.pdf"), w = 4, h = 3.5)
    # print(g2)
    # dev.off()
    # 
    library(ggpubr)
    p <- ggarrange(g, g2,
                   common.legend = TRUE, legend = "right", widths = c(0.5,0.5))
    pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/enrichment", gsub('\\/','_',type),'_',trait,"mainv1.pdf"), w = 8, h = 4)
    print(p)
    dev.off()
  }
}


###### enrichment shared positions #####
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

##### plot shared positions
hypo <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hypo_shared_CI.continous.2.rds')
hyper <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hyper_shared_CI.continous.2.rds')

colors_traits <- list('Age'=c('#3D7CD0','#B4D6F6'),
                      'Sex'=c('#3B734E','#89AA94'),
                      'Ancestry'=c('#F0AE21','#F9DE8B'))

library(ggpubr)
  for (trait in c('Ancestry','Sex','Age')) {

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
    pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/enrichment",'_',trait,"_shared.2.pdf"), w = 8, h = 4)
    print(p)
    dev.off()
  }


##### enrichment overlap ageing - ancestry colon
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
EUR_YOUNG <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hyperEUR_hypoOLD_CI.continous.rds')
AFR_OLD <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hyperAFR_hyperOLD_CI.continous.rds')

library(ggpubr)
 
  hypo_d <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), EUR_YOUNG)
  hyper_d <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), AFR_OLD)
  hyper_hypo <- rbind(hypo_d, hyper_d)
  hyper_hypo$sig[hyper_hypo$sig =="Sig"] <- "FDR < 0.05"
  hyper_hypo$sig[hyper_hypo$sig =="not Sig"] <- "FDR >= 0.05"
  hyper_hypo$sig <- factor(hyper_hypo$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
  hyper_hypo$region <- factor(hyper_hypo$region, levels=rev(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts")))
  hyper_hypo$type[hyper_hypo$type =="EUR_YOUNG"] <- "EUR - Young"
  hyper_hypo$type[hyper_hypo$type =="AFR_OLD"] <- "AFR - Old"
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
  pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/enrichment_overlap_Age_ancestry_colon.pdf"), w = 8, h = 4)
  print(p)
  dev.off()
