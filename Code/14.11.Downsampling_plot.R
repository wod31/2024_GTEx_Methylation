#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to plot expression downsampling
# @software version: R=4.2.2


#Set path 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

variable <- 30
variable <- 40
variable <- 50
variable <- 35
medians <- readRDS(paste0("Downsampling/downsampling_", variable, ".rds"))
medians <- round(medians)

library(ComplexHeatmap)
library(RColorBrewer)

medians[c("Ovary", "Prostate", "Testis"),"Sex"] <- NA
no_nas <- medians
no_nas <- apply(apply(apply(no_nas, 2, as.numeric), 2, round), 2, prettyNum, big.mark=",")
no_nas[no_nas=="NA"] <- ""

pdf(paste0("Plots/Downsampling_", variable,".pdf"), height = 2.5, width = 4)
Heatmap(medians, col = brewer.pal(9,"BuPu")[c(1,4:7)],
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(no_nas[i, j], x, y, gp = gpar(fontsize = 9))},
        na_col = "white",
        cluster_rows = F,
        column_names_side = "top",
        cluster_columns = F,
        name = "Nº DMPs",
        column_names_rot = 0,
        column_title = paste("Number of samples =", variable),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        row_names_side = "left",)
dev.off()
#Sort tissues



# #To main:
# library(reshape2)
# library(ggplot2)
# 
# df <- melt(medians)
# df$Var2 <- as.factor(df$Var2)
# cols <- tissue_info$color
# names(cols) <- tissue_info$Name
# 
# g <- ggplot(df, aes(Var2, value, col=Var1)) + geom_point() +
#   scale_color_manual(values=cols) +
#   geom_line(aes(group = Var1)) + theme_bw() +
#   theme(legend.title = element_blank(),
#         axis.text = element_text(colour="black", size=13),
#         legend.text = element_text(colour="black", size=13),
#         axis.title = element_text(size=16),
#         # legend.spacing.y = unit(-0.05, "cm"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(colour = "black", linewidth=1)) +
#   ylab("Number of DEGs") + xlab("Number of samples")
# 
# pdf("../figures/figures/figure_s1/downsampling_main.pdf", height = 2.5, width = 10.2)
# g
# dev.off()
