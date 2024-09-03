#!/usr/bin/env Rscript

#Set path 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
setwd("..")

library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
suppressPackageStartupMessages(library(circlize))

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "KidneyCortex", "Testis", "WholeBlood", "MuscleSkeletal")
sex_tissues <- c("Ovary", "Prostate", "Testis")
variable <- "Age"
variable <- "Sex"
variable <- "Ancestry"
table <- matrix(nrow = length(tissues), ncol = 5, dimnames = list(tissues, c("Overlap", "Only_DMP", "Only_DVP", "p-value", "OR")))
for(tissue in tissues){
  dmp <- readRDS(paste0("Tissues/", tissue, "/DML_results_5_PEERs_continous.rds"))
  bg <- rownames(dmp$AGE) #all tested cpgs
  if(variable=="Ancestry"){
    dmp <- rownames(dmp$EURv1[dmp$EURv1$adj.P.Val<0.05,])
  } else if(variable=="Age"){
    dmp <- rownames(dmp$AGE[dmp$AGE$adj.P.Val<0.05,])
  } else if(variable=="Sex"){
    if(tissue %in% sex_tissues){next}
    dmp <- rownames(dmp$SEX2[dmp$SEX2$adj.P.Val<0.05,])
  }
  dvp_a <- (readRDS(paste0("Tissues/", tissue, "/DVP_", variable, ".rds")))
  dvp <- rownames(dvp_a[dvp_a$Adj.P.Value<0.05,])
  table[tissue,1] <- length(dvp[dvp %in% dmp])
  table[tissue,2] <- length(dmp[!dmp %in% dvp])
  table[tissue,3] <- length(dvp[!dvp %in% dmp])
  none <- bg[!bg %in% c(table[tissue,1], table[tissue,2], table[tissue,3])]
  matrix <- matrix(c(table[tissue,1], table[tissue,2], table[tissue,3], length(none)), nrow=2) #contingency table to do fisher test
  t <- fisher.test(matrix)
  table[tissue,4] <- t$p.value
  table[tissue,5] <- t$estimate
}
#put the rectangle when it is significant, as it is relevant for age and sex

colnames(table)[1:3] <- c("Overlap", "Only DMP", "Only DVP")

#Change tissue names:

pdf("Plots/DMP_DVP_Overlap.pdf", width = 3.5, height = 3)
Heatmap(table[,c(1:3)], col = brewer.pal(9,"BuPu")[c(1,4:7)],
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (j == 1 && i %in% which(table[,4]<0.05)) { 
            grid.rect(x, y, width = width, height = height, gp = gpar(fill=fill, color = "black"))
            grid.text(table[i, j], x, y, gp = gpar(fontsize = 9))
          } else {
          grid.text(table[i, j], x, y, gp = gpar(fontsize = 9))}},
        na_col = "white",
        cluster_rows = F,
        column_names_side = "top",
        cluster_columns = F,
        name = "Nº CpGs",
        # column_names_rot = 0,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        row_names_side = "left")
dev.off()


#Plot examples (in lung)
lung_betas <- readRDS("Tissues/Lung/data.rds")
lung_metadata <- readRDS("Tissues/Lung/metadata.rds")
# admixture_ancestry <- read.table('../../../../../scratch/bsc83/bsc83535/GTEx/v8/genotype_data/admixture_inferred_ancestry.txt')
admixture_ancestry <- read.table('../../../../../scratch/bsc83/MN4/bsc83/bsc83535/GTEx/v8/genotype_data/admixture_inferred_ancestry.txt')
colnames(admixture_ancestry) <- c('SUBJID','AFRv1','EURv1','inferred_ancestry','AFRv2','EURv2')
lung_metadata <- merge(lung_metadata, admixture_ancestry[,c("SUBJID","EURv1")], by='SUBJID')

dmp <- readRDS(paste0("Tissues/Lung/DML_results_5_PEERs_continous.rds"))
dmp <- dmp$EURv1[dmp$EURv1$adj.P.Val<0.05,]
dvp <- readRDS(paste0("Tissues/Lung/DVP_Ancestry.rds"))
common <- rownames(dvp)[rownames(dvp) %in% rownames(dmp)]
only_dmp <- rownames(dmp)[!rownames(dmp) %in% rownames(dvp)]
only_dvp <- rownames(dvp)[!rownames(dvp) %in% rownames(dmp)]
#get examples with high logFCs
rownames(dmp[abs(dmp$logFC)>4,]) %in% common #cg26084667
rownames(dmp[abs(dmp$logFC)>3.5,]) %in% only_dmp #cg18232235, cg09255157, cg23386212, cg08477332, cg12262617, cg12080266

dvp[rownames(dvp) %in% only_dvp & abs(dvp$SampleVar)>3.5,] #cg16512708, cg11106864, cg17343385, cg05875700 in one direction, or cg21245975, cg12209881 in other
dvp[rownames(dvp) %in% only_dvp & dvp$DiffLevene<(-1.2),] #cg16512708, cg11106864 in one direction, or cg21245975, cg12209881, cg00874599 in other
#DiffLevene sign positive more variability in europeans
#More candidates to have the same direction as the other examples: cg11143152, cg06393909, cg11143152

#DMP and DVP
to_plot <- lung_betas[rownames(lung_betas)=="cg26084667",]

#Only DMP
# to_plot <- lung_betas[rownames(lung_betas)=="cg18232235",]
to_plot <- lung_betas[rownames(lung_betas)=="cg09255157",]

#Only DVP
# to_plot <- lung_betas[rownames(lung_betas)=="cg11106864",]
to_plot <- lung_betas[rownames(lung_betas)=="cg00874599",]


to_plot <- t(to_plot)
colnames(to_plot) <- "betas"
to_plot <- merge(to_plot, lung_metadata, by.x="row.names", by.y="SUBJID")
to_plot$betas <- as.numeric(to_plot$betas)
to_plot$Ancestry <- "EUR"
to_plot$Ancestry[to_plot$EURv1<0.5] <- "AFR" 

g <- ggplot(to_plot, aes(Ancestry, betas)) +     
  geom_violin(aes(fill = Ancestry), col = "black", scale = "width") +
  geom_jitter(col = "black",
              alpha = 0.5,
              size = 0.8) +
  theme_bw() +
  theme(        
    # strip.background = element_rect(fill=traits_cols["Ancestry"]),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", linewidth=1)
  ) +
  stat_summary(geom = "point",
    fun = "mean", col = "black",
    size = 3, shape = 24, fill = "red"
  )

ggsave("Plots/Example_DMP_DVP.pdf", g, width = 1.8, height = 2.2)
ggsave("Plots/Example_DMP_only.pdf", g, width = 1.8, height = 2.2)
ggsave("Plots/Example_DVP_only.pdf", g, width = 1.8, height = 2.2)
