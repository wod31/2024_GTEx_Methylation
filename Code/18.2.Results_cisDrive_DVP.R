#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Plot results of cis-driven DVPs
# @software version: R=4.2.2

set.seed(1)

# Libraries ----
library(RColorBrewer)
library(ComplexHeatmap)
library(ggplot2)
library(cowplot)
# Box plots
library(rcompanion)
library(gtools)

# ---- Data ---- ####
# Tissues ---
first_dir <- "~/marenostrum/"
first_dir <- "~/Documents/mn4/projects/bsc83/"

project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")

tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Jose/00_Data/Tissue_info_whole.rds"))
# 
tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

# Tissue metadata ----
metadata <- lapply(tissues, function(tissue) readRDS(paste0(project_path,"Tissues/", tissue, "/metadata.rds")))
names(metadata) <- tissues

# 1.3 cis-driven classification ----
inpath <- paste0(first_dir, "/Projects/GTEx_v8/Methylation/Tissues/")
model <- "DVP"
input <- "residuals"
correlation <- "0.8"

for(tissue in tissues){ #For me this is the same in DMP and DVP
  if(!file.exists(paste0(inpath, tissue, '/Testing_', model, '_', input, "_", correlation, '_Ancestry.Classification_summary.qval005.rds'))){print(tissue)}
}
d <- lapply(tissues, function(tissue)
  readRDS(paste0(inpath, tissue, '/Testing_', model, '_', input, "_", correlation, '_Ancestry.Classification_summary.qval005.rds')))
names(d) <- tissues

# Create dataframe --
data <- do.call(rbind.data.frame,d)
# rownames(data) <- tissues
# colnames(data) <- as.character(names(d[[tissues[1]]]))

# DMP/DVP classified as either cis-driven or not cis-driven --
results <-  lapply(tissues, function(tissue) 

readRDS(paste0(inpath, tissue, '/Testing_', model, '_', input, "_", correlation, '_Ancestry.Classified.qval005.rds')))
names(results) <- tissues

# Expand dataframe --
data <- cbind.data.frame(data,
                         sapply(tissues, function(tissue) sum(results[[tissue]]$Class=="Cis-driven", na.rm = T)),
                         sapply(tissues, function(tissue) sum(results[[tissue]]$Class=="Not_cis-driven", na.rm = T)),
                         sapply(tissues, function(tissue) sum(is.na(results[[tissue]]$Class), na.rm = T))
)
colnames(data)[10] <- "Cis-driven"
colnames(data)[11] <- "Not_cis-driven"
colnames(data)[12] <- "Too_Many_mQTL"
data[,12] <- data[,8] + data[,12] #The ones not modelled [,8] had also many mQTLs. And actually the [,12] could not perform anova due to many collinearities, even if we are very restringent

# Ancestry DEG data.frame  ----
dd <- data[,c(2,3,10,11,12)]

# Data for bar plot --
ddp <- cbind.data.frame("not mGene" = dd[,1] - dd[,2], # not sGenes
                        "not cis-driven" = dd$`Not_cis-driven`,
                        "cis-driven" = dd$`Cis-driven`,
                        "not classified" = dd$`Too_Many_mQTL` # ot modelled: No isQTL with MAF 001, no isQTL with var, no isQTL with no dependance
) 
rownames(ddp) <- tissue_info[tissues,'tissue_abbrv']
splic.ddp <- ddp
d2 <- as.data.frame(sapply(1:4, function(i) apply(splic.ddp, 1, function(x) x[i]/sum(x))))
colnames(d2) <- colnames(splic.ddp)
rownames(splic.ddp) <- rownames(ddp)
splic.ddp$tissue <- rownames(ddp)
library(reshape2)
d2$tissue <- tissue_info[tissues,'tissue_abbrv']
d1 <- melt(splic.ddp)
colnames(d1)[3] <- "count"
d11 <- melt(d2)
colnames(d11)[3] <- "proportion"
cis_data <- cbind.data.frame(d1, d11[,3])
colnames(cis_data)[4] <- "proportion"
cis_data$proportion <- 100*cis_data$proportion

cis_data$tissue <- factor(cis_data$tissue, levels = rev(tissue_info[tissues,'tissue_abbrv']), order = T)
cis_data$variable <- factor(cis_data$variable, levels = rev(levels(cis_data$variable)), order = T)
cols <- c(brewer.pal(11,"PRGn")[c(2)], "Red",brewer.pal(11,"PRGn")[c(9)],"light grey")

pdf(paste0(project_path, "/Plots/Cis_driven_", model, "_", correlation, "_", input, ".qval005.pdf"),
    width = 8, height = 4)
names(cols) <- c("not mGene", "not cis-driven", "cis-driven", "not classified")
ggplot(cis_data, mapping = aes(x = proportion, 
                               fill = variable, #actor(var_2, levels = rev(x_levs)), 
                               y = tissue, #factor(var_1, levels = rev(y_levs)),
                               label = count)) +
  geom_bar( stat = "identity") + 
  geom_text(size = 3, position = position_stack(vjust = 0.5)) + # Add number labelsfacet_grid(~trait) +
  ylab("") + xlab(paste(model, "(%)")) +
  labs(fill = "") +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10))
dev.off()  



#From the ones that we can model, get barplots
pdf(paste0(project_path, "/Plots/Cis_driven_", model, "_", correlation, "_", input, "_with_mQTLs.qval005.pdf"),
    width = 8, height = 4)
ggplot(cis_data[cis_data$variable!="not mGene",], mapping = aes(x = count, 
                               fill = variable, #actor(var_2, levels = rev(x_levs)), 
                               y = tissue, #factor(var_1, levels = rev(y_levs)),
                               label = count)) +
  geom_bar( stat = "identity") + 
  geom_text(size = 3, position = position_stack(vjust = 0.5)) + # Add number labelsfacet_grid(~trait) +
  ylab("") + xlab(paste("Number of ",model)) +
  labs(fill = "") +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10))
dev.off()


for(tissue in tissues){
  test <- readRDS(paste0(inpath, tissue, '/Testing_', model, '_', input, "_", correlation, '_Ancestry.Classified.qval005.rds'))
}

d <- lapply(tissues, function(tissue) {
  file <- readRDS(paste0(inpath, tissue, '/Testing_', model, '_', input, "_", correlation, '_Ancestry.Classified.qval005.rds'))
  if(nrow(file)>0){
    print(tissue)
    cbind(file, tissue)
  }
})

data <- do.call(rbind.data.frame,d)
data <- data[!is.na(data$anova.adj.P.Val),]
data <- data[data$anova.adj.P.Val>=0.05,]
table(data$dvp)[table(data$dvp)>2] #cg12052018
table(data$dvp)[table(data$dvp)>1] #others: cg00349687 cg01808640 cg02044222 cg02921257



### dotplot ####
head(dd)
dd$Tissue <- rownames(dd)
# rownames(dd) <- tissue_info[tissues,'tissue_abbrv']
dd$perc_cis_driven_2 <- 100*(dd$`Cis-driven`/(dd$`Cis-driven`+dd$`Not_cis-driven`))
library(reshape2)
# dd$tissue <- tissue_info[tissues,'tissue_abbrv']
d1 <- melt(dd)

tissues_cols <- tissue_info[, 3]
names(tissues_cols) <- tissue_info[, 1]
# tissues_cols <- tissues_cols[tissues]
# names(tissues_cols) <- tissue_info$tissue_abbrv

to_plot <- d1[d1$variable %in% c('perc_cis_driven_2'),]
# to_plot <- to_plot[c(1:4),]
pdf('Plots/cis_driven_DVPs_dotplot.qval005.pdf', width = 1.2, height = 2)
ggplot(data = to_plot,
       aes(x = variable,
           y = value,
           color =  Tissue),
) +geom_boxplot(col = "black",
                fill = "white",
                outlier.shape = NA,
                width = 0.8) +
  geom_jitter(alpha = 1,
              size = 2) +
  xlab("") +
  ylab("% cis-driven DVPs") +
  # ylim(c(65, 85)) +
  ylim(c(0, 100)) +
  scale_x_discrete(breaks=c("perc_cis_driven_2"),
                   labels=c(""))+
  scale_fill_manual(values = tissues_cols) +
  scale_color_manual(values = tissues_cols) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.title.y = ggtext::element_markdown(size=12),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none") 
dev.off()


to_plot <- dd[,c(1,2)]
to_plot$percentage <- to_plot[,2]/to_plot[,1]*100
to_plot$Tissue <- rownames(to_plot)
to_plot$dummy <- "Dummy"
pdf('Plots/cis_driven_DVPs_with_mQTLs.qval005.pdf', width = 1.2, height = 2)
ggplot(data = to_plot,
       aes(x = dummy,
           y = percentage,
           color =  Tissue),
) +geom_boxplot(col = "black",
                fill = "white",
                outlier.shape = NA,
                width = 0.8) +
  geom_jitter(alpha = 1,
              size = 2) +
  xlab("") +
  ylab("% DVPs with mQTLs") +
  ylim(c(0, 100)) +
  scale_x_discrete(breaks=c("Dummy"),
                   labels=c(""))+
  scale_fill_manual(values = tissues_cols) +
  scale_color_manual(values = tissues_cols) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.title.y = ggtext::element_markdown(size=12),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none") 
dev.off()
