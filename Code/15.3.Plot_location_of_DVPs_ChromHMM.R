#!/usr/bin/env Rscript

#Set path 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

data <- readRDS('Tissues/DVP_enrichment_chromhmm.rds')

#Put names hyper and hypo
for (name in names(data)) {for (trait in names(data[[name]])) { for(region in names(data[[name]][[trait]])){names(data[[name]][[trait]][[region]]) <- c("Hyper", "Hypo")}}}


#transform into data frame

#Parsing data. It could be optimized
nested_list_to_df <- function(nested_list) {
  if (length(nested_list)==8) { #The last list is what I want as data frame, and it is the only one with 8 elements
    return(data.frame(nested_list))
  } else {
    inner_dfs <- lapply(nested_list, nested_list_to_df)
    result_df <- do.call(rbind, inner_dfs)
    return(result_df)
  }
}

# Convert nested list to data frame
result_df <- nested_list_to_df(data)
ancestry <- result_df[grep("Ancestry", rownames(result_df)),]
age <- result_df[grep("Age", rownames(result_df)),]
sex <- result_df[grep("Sex", rownames(result_df)),]

ancestry$adj.p.value <- p.adjust(ancestry$p_value)
age$adj.p.value <- p.adjust(age$p_value)
sex$adj.p.value <- p.adjust(sex$p_value)

library(ggplot2)
pdf('Plots/DVP_heatmap_ancestry.chromhmm.pdf', width = 5, height = 4)
ggplot(ancestry, aes(x=tissue, y=region, size=-log10(adj.p.value), fill=log2(oddsRatio))) +
  geom_point(alpha=0.9, shape=21, color="black") +
  scale_size(range = c(0.5, 12), name="-log10(FDR)", limits = c(0,max(-log10(ancestry$adj.p.value)))) +
  scale_fill_gradient2(low="#1982C4", high="#690500",  mid = '#F4EBBE', midpoint = 0.7,name='log2(oddsRatio)')+
  theme(legend.position="right") +
  ylab("") +
  xlab("") +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        panel.grid = element_line(colour = 'light grey'),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf('Plots/DVP_heatmap_age.chromhmm.pdf', width = 5, height = 4)
ggplot(age, aes(x=tissue, y=region, size=-log10(adj.p.value), fill=log2(oddsRatio))) +
  geom_point(alpha=0.9, shape=21, color="black") +
  scale_size(range = c(0.5, 12), name="-log10(FDR)", limits = c(0,max(-log10(age$adj.p.value)))) +
  scale_fill_gradient2(low="#1982C4", high="#690500",  mid = '#F4EBBE', midpoint = 0.7,name='log2(oddsRatio)')+
  theme(legend.position="right") +
  ylab("") +
  xlab("") +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        panel.grid = element_line(colour = 'light grey'),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf('Plots/DVP_heatmap_sex.chromhmm.pdf', width = 5, height = 4)
ggplot(sex, aes(x=tissue, y=region, size=-log10(adj.p.value), fill=log2(oddsRatio))) +
  geom_point(alpha=0.9, shape=21, color="black") +
  scale_size(range = c(0.5, 12), name="-log10(FDR)", limits = c(0,max(-log10(sex$adj.p.value)))) +
  scale_fill_gradient2(low="#1982C4", high="#690500",  mid = '#F4EBBE', midpoint = 0.7,name='log2(oddsRatio)')+
  theme(legend.position="right") +
  ylab("") +
  xlab("") +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        panel.grid = element_line(colour = 'light grey'),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
