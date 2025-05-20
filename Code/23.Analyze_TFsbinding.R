#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Create summary plots for the binding of TFBS on DMPs
# @software version: R=4.2.2

library(ggplot2)
library(dplyr)
library(tidyr)

### get beds ####
annotation <- read.csv("marenostrum_scratch///GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")
#first_dir <- "/gpfs/projects/bsc83/"
first_dir <- "marenostrum/"

project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Age", "Ancestry", "BMI")
traits_to_use <- c('EURv1','AGE','BMI')

results_DML <- lapply(tissues, function(tis) 
  readRDS(paste0("marenostrum/Projects/GTEx_v8/Methylation/Tissues/",tis,"/DML_results_5_PEERs_continous.rds")))
names(results_DML) <- tissues

ann_bed <- annotation[!is.na(annotation$MAPINFO) & !is.na(annotation$CHR),] %>%
  dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>% distinct()
ann_bed$chrom <- paste0('chr',ann_bed$chrom)
head(ann_bed)
ann_bed$start <- ann_bed$start-1

for (tissue in tissues) {
  for (variable in traits_to_use) {
    hyper <- ann_bed[ann_bed$name %in% rownames(results_DML[[tissue]][[variable]][results_DML[[tissue]][[variable]]$adj.P.Val<0.05 & results_DML[[tissue]][[variable]]$logFC>0,]),]
    hypo <- ann_bed[ann_bed$name %in% rownames(results_DML[[tissue]][[variable]][results_DML[[tissue]][[variable]]$adj.P.Val<0.05 & results_DML[[tissue]][[variable]]$logFC<0,]),]
    write.table(hyper, paste0('marenostrum/Projects/GTEx_v8/Methylation/Data/up_',variable,'_cpgs_',tissue,'.bed'), sep = '\t', col.names = F, row.names = F, quote = F)
    write.table(hypo, paste0('marenostrum/Projects/GTEx_v8/Methylation/Data/down_',variable,'_cpgs_',tissue,'.bed'), sep = '\t', col.names = F, row.names = F, quote = F)
    tested <- ann_bed[ann_bed$name %in% rownames(results_DML[[tissue]][[variable]]),]
    write.table(tested, paste0('marenostrum/Projects/GTEx_v8/Methylation/Data/tested_',variable,'_cpgs_',tissue,'.bed'), sep = '\t', col.names = F, row.names = F, quote = F)
    
  }
}

######## Analize chip-atlas results #####
colors_traits <- list('AGE'=c('#3D7CD0','#B4D6F6'),
                      'SEX2'=c('#3B734E','#89AA94'),
                      'EURv1'=c('#F0AE21','#F9DE8B'))

#Plot top enriched TFBS in ggarrange
parsing_dummy_variable <- function(value) {
  return(gsub("_.*", "", value)) #Replace everything after the first _
}

library(ggh4x)
traits_cols <- c('#C49122','#4B8C61','#70A0DF')
names(traits_cols) <- c('EURv1','SEX2','AGE')
strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols))

for(tissue in tissues){
  print(tissue)
  data <- read.csv(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/TFBS/tfbs_results_", tissue, ".csv"))
  signif <- data[data$adj_p_val<0.05 & data$direction=="Enriched",]
  hyper_plot <- signif[signif$methylation=="Hyper",]
  hyper_plot$logOR <- log2(hyper_plot$odds_ratio+0.001)
  hypo_plot <- signif[signif$methylation=="Hypo",]
  hypo_plot$logOR <- log2(hypo_plot$odds_ratio+0.001)
  hypo_plot$logOR <- -hypo_plot$logOR
  # plot <- signif[order(signif$adj_p_val),]
  plot <- rbind(hyper_plot, hypo_plot) #plot together the top 10 hyper and top 10 hypo?
  if(nrow(plot)==0){next}
  plot_top_10 <- plot %>%
    group_by(variable, methylation) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  #I create a dummy variable because the same TF can be in different variables and I want a different order in each facet
  plot_top_10$TF_dummy <- paste0(plot_top_10$TF, "_", plot_top_10$variable, "_", plot_top_10$methylation)
  plot_top_10$TF_dummy <- factor(plot_top_10$TF_dummy, levels = rev(plot_top_10$TF_dummy))
  
  plot_top_10$variable <- factor(plot_top_10$variable, levels = c('EURv1','SEX2','AGE'))
  
  trues <- names(traits_cols) %in% plot_top_10$variable
  strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols[trues]))
  
  g_plot <- ggplot(plot_top_10, aes(x=log2(odds_ratio), y=TF_dummy, col=methylation)) +
    geom_errorbar(aes(xmin=log2(CI_low), xmax=log2(CI_high)), width=.4) +
    geom_vline(xintercept = 0) +
    geom_point(size=2.5) + ylab('') + theme_bw() +
    scale_colour_manual(values=c("#b23c29", "#1f6db6")) +
    xlab("Log(Odds ratio)") +
    theme(#legend.position = "none",
      axis.text.x = element_text(colour="black", size=11),
      axis.text.y = element_text(colour="black", size=13),
      legend.text = element_text(colour="black", size=12),
      axis.title.x = element_text(size=13),
      legend.spacing.y = unit(-0.05, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", linewidth=1)) +
    facet_wrap2(~variable, scales = "free_y", strip = strip, nrow = 1) + #Function facet_wrap2 is different from facet_wrap only to color
    scale_y_discrete(labels = parsing_dummy_variable)
  ggsave(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/TFBS/Top_", tissue, ".pdf"), g_plot, width = 4 + sum(trues), height = 4)
  
  for (trait in c('EURv1','SEX2','AGE')) {
    pdf(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/TFBS/Top_Enrichment_", tissue,'_',trait, ".pdf"), width = 5, height = 8)
    print(ggplot(plot[plot$variable == trait,]) + geom_bar(aes(x=logOR,y=reorder(TF, -logOR), fill=adj_p_val), width = 0.6, stat = 'identity') +
      theme_classic() + xlab("log2(meanOR)") + ylab("") +
      scale_fill_gradient2( high = colors_traits[[trait]][1] ,low=colors_traits[[trait]][2], name='meanFDR', midpoint = 0)+
      theme( axis.text.x = element_text(colour="black", size=13),
             axis.text.y = element_text(colour="black", size=13),
             axis.title.x = element_text(size=16)))
    dev.off()
  }
  
}


