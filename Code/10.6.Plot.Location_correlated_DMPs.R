#### plot enrichments CI 

library(ggh4x)
tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Jose/00_Data/Tissue_info_whole.rds"))
tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]


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


hypo <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/enrichment_anno_hyper_correlated.rds')
hyper <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/enrichment_anno_hypo_correlated.rds')

# sex_tissues <- c('Ovary','Prostate','Testis')
# for (tissue in names(hypo)) {
#   for (trait in names(hypo$Lung)) {
#     if (tissue %in% sex_tissues & trait == "SEX2") {
#       print(NA)
#     } else {
#       hypo_d <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), hypo, tissue, trait)
#       hyper_d <- read_data(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"), hyper, tissue, trait)
#       hyper_hypo <- rbind(hypo_d, hyper_d)
#       hyper_hypo$sig[hyper_hypo$sig =="Sig"] <- "FDR < 0.05"
#       hyper_hypo$sig[hyper_hypo$sig =="not Sig"] <- "FDR >= 0.05"
#       hyper_hypo$sig <- factor(hyper_hypo$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
#       hyper_hypo$region <- factor(hyper_hypo$region, levels=rev(c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts")))
#       hyper_hypo$type[hyper_hypo$type =="hypo"] <- "Hypomethylation"
#       hyper_hypo$type[hyper_hypo$type =="hyper"] <- "Hypermethylation"
#       g <- ggplot(hyper_hypo, aes(x=(oddsRatio), y=region, colour=type, alpha=sig)) +
#         geom_errorbar(aes(xmin=(CI_down), xmax=(CI_up)), width=.3) +
#         geom_vline(xintercept = 1) +
#         xlim(0,20) + #Only for Lung to show the 0
#         geom_point(size=3) + ylab('') + theme_bw() +
#         scale_colour_manual(values=c("#CC6677", "#88CCEE")) +
#         xlab("Odds ratio") +
#         scale_alpha_discrete(range = c(0.4, 1), drop = FALSE) +
#         theme(legend.title = element_blank(),
#               axis.text.x = element_text(colour="black", size=13),
#               axis.text.y = element_text(colour="black", size=14),
#               legend.text = element_text(colour="black", size=13),
#               axis.title.x = element_text(size=16),
#               legend.spacing.y = unit(-0.05, "cm"),
#               panel.grid.major = element_blank(),
#               panel.grid.minor = element_blank(),
#               panel.border = element_rect(colour = "black", linewidth=1)) +
#         scale_y_discrete(breaks=c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"),
#                          labels=c("Enhancer","Enhancer Bivalent","Heterochromatin","Quiescent","Repressed Polycomb","TSS","TSS Bivalent","Transcription","ZNF & Repeats"))# + xlim(0, 3)
#       pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/genomic_location_", tissue,'_',trait,".pdf"), w = 6, h = 3.5)
#       print(g)
#       dev.off()
#       
#       #Plot sample sizes:
#       
#       g2 <- ggplot(hyper_hypo) + geom_col(aes(sample_size, region, fill=type), width = 0.6) +
#         theme_classic() + xlab("Number of DMPs") + ylab("") +
#         scale_fill_manual(values=c("#CC6677", "#88CCEE")) +
#         theme(legend.position = "none",
#               axis.text.x = element_text(colour="black", size=13),
#               axis.text.y = element_text(colour="black", size=14),
#               axis.title.x = element_text(size=16)) +
#         scale_y_discrete(breaks=c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"),
#                          labels=c("Enhancer","Enhancer Bivalent","Heterochromatin","Quiescent","Repressed Polycomb","TSS","TSS Bivalent","Transcription","ZNF & Repeats")) #+
#       #scale_x_continuous(breaks=c(0, 20000, 40000)) #Only for lung
#       pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/genomic_location_", tissue,'_',trait,"_sample_size.pdf"), w = 4, h = 3.5)
#       print(g2)
#       dev.off()
#     }
#   }
# }


####### Jose's code to plot both number and enrichment one plot per tissue and traitt 
### gene body , promoter, enhancer
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


# hypo <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hypo_batch_CI.simple.continous.rds')
# hyper <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hyper_batch_CI.simple.continous.rds')

colors_traits <- list('AGE'=c('#3D7CD0','#B4D6F6'),
                      'SEX2'=c('#3B734E','#89AA94'),
                      'EURv1'=c('#F0AE21','#F9DE8B'))

library(ggpubr)
sex_tissues <- c('Ovary','Prostate','Testis')
for (type in c("promoter","enhancer","gene_body")) {
  for (trait in c('AGE','SEX2','EURv1')) {
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
    hyper_hypo$type[hyper_hypo$type =="hypo"] <- "Hypomethylation"
    hyper_hypo$type[hyper_hypo$type =="hyper"] <- "Hypermethylation"
    hyper_hypo$oddsRatio[hyper_hypo$oddsRatio == Inf] <- max(hyper_hypo$oddsRatio[hyper_hypo$oddsRatio != Inf])
    hyper_hypo$CI_up[hyper_hypo$CI_up == Inf] <- max(hyper_hypo$CI_up[hyper_hypo$oddsRatio != Inf])
    hyper_hypo$CI_down[hyper_hypo$CI_down == -Inf] <- min(hyper_hypo$CI_down[hyper_hypo$oddsRatio != -Inf])
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
            panel.border = element_rect(colour = "black")) #+
    # scale_y_discrete(breaks=c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"),
    #                  labels=c("Enhancer","Enhancer Bivalent","Heterochromatin","Quiescent","Repressed Polycomb","TSS","TSS Bivalent","Transcription","ZNF & Repeats"))# + xlim(0, 3)
    # pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/genomic_location_", gsub('\\/','_',type),'_',trait,".correlated.pdf"), w = 6, h = 3.5)
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
    # pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/genomic_location_", gsub('\\/','_',type),'_',trait,"_sample_size.correlated.pdf"), w = 4, h = 3.5)
    # print(g2)
    # dev.off()
    # 
    p <- ggarrange(g, g2, labels = c("A", "B"),
                   common.legend = TRUE, legend = "right", widths = c(0.8,0.3))
    pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/enrichment", gsub('\\/','_',type),'_',trait,".correlated.pdf"), w = 8, h = 4)
    print(p)
    dev.off()
  }

# library(ggpubr)
# sex_tissues <- c('Ovary','Prostate','Testis')
# for (type in c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts")) {
#   for (trait in names(hypo$Lung)) {
#     if (trait == 'SEX2') {
#       tissues_to_test <- tissues[!tissues %in% sex_tissues]
#     } else {
#       tissues_to_test <- tissues
#     }
#     hypo_d <- read_data(tissues_to_test, hypo, type, trait)
#     hyper_d <- read_data(tissues_to_test, hyper, type, trait)
#     hyper_hypo <- rbind(hypo_d, hyper_d)
#     hyper_hypo$sig[hyper_hypo$sig =="Sig"] <- "FDR < 0.05"
#     hyper_hypo$sig[hyper_hypo$sig =="not Sig"] <- "FDR >= 0.05"
#     hyper_hypo$sig <- factor(hyper_hypo$sig, levels=c("FDR >= 0.05", "FDR < 0.05"))
#     hyper_hypo$tissue <- factor(hyper_hypo$tissue, levels=tissues_to_test)
#     hyper_hypo$type[hyper_hypo$type =="hypo"] <- "Hypomethylation"
#     hyper_hypo$type[hyper_hypo$type =="hyper"] <- "Hypermethylation"
#     g <- ggplot(hyper_hypo, aes(x=log2(oddsRatio), y=tissue, colour=type, alpha=sig)) +
#       geom_errorbar(aes(xmin=log2(CI_down), xmax=log2(CI_up)), width=.3) +
#       geom_vline(xintercept = 0) +
#       #xlim(0,20) + #Only for Lung to show the 0
#       geom_point(size=3) + ylab('') + theme_bw() +
#       scale_colour_manual(values=c("#CC6677", "#88CCEE")) +
#       xlab("log2(Odds ratio)") +
#       scale_alpha_discrete(range = c(0.4, 1), drop = FALSE) +
#       theme(legend.title = element_blank(),
#             axis.text.x = element_text(colour="black", size=13),
#             axis.text.y = element_text(colour="black", size=14),
#             legend.text = element_text(colour="black", size=13),
#             axis.title.x = element_text(size=16),
#             legend.spacing.y = unit(-0.05, "cm"),
#             panel.grid.major = element_blank(),
#             panel.grid.minor = element_blank(),
#             panel.border = element_rect(colour = "black", linewidth=1)) #+
#     # scale_y_discrete(breaks=c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"),
#     #                  labels=c("Enhancer","Enhancer Bivalent","Heterochromatin","Quiescent","Repressed Polycomb","TSS","TSS Bivalent","Transcription","ZNF & Repeats"))# + xlim(0, 3)
#     pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/genomic_location_", gsub('\\/','_',type),'_',trait,".pdf"), w = 6, h = 3.5)
#     print(g)
#     dev.off()
#     
#     #Plot sample sizes:
#     
#     g2 <- ggplot(hyper_hypo) + geom_col(aes(sample_size, tissue, fill=type), width = 0.6) +
#       theme_classic() + xlab("Number of DMPs") + ylab("") +
#       scale_fill_manual(values=c("#CC6677", "#88CCEE")) +
#       theme(legend.position = "none",
#             axis.text.x = element_text(colour="black", size=13),
#             axis.text.y=element_blank(),  #remove y axis labels,
#             axis.title.x = element_text(size=16)) +
#       scale_x_continuous(n.breaks=3)
#     #   scale_y_discrete(breaks=c("Enh","EnhBiv","Het","Quies","ReprPC","TSS","TssBiv","Tx","ZNF/Rpts"),
#     #                    labels=c("Enhancer","Enhancer Bivalent","Heterochromatin","Quiescent","Repressed Polycomb","TSS","TSS Bivalent","Transcription","ZNF & Repeats")) #+
#     # #scale_x_continuous(breaks=c(0, 20000, 40000)) #Only for lung
#     pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/genomic_location_", gsub('\\/','_',type),'_',trait,"_sample_size.pdf"), w = 4, h = 3.5)
#     print(g2)
#     dev.off()
#     
#     p <- ggarrange(g, g2, labels = c("A", "B"),
#                    common.legend = TRUE, legend = "right", widths = c(0.8,0.3))
#     pdf(file = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/chromhmm/enrichment", gsub('\\/','_',type),'_',trait,".pdf"), w = 8, h = 4)
#     print(p)
#     dev.off()
#   }
# }
# 
