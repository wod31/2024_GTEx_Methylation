#!/usr/bin/env Rscript
# @Author: Winona Oliveros
# @E-mail: winona.oliveros@bsc.es
# @Description: Code to do supplementary plots of additive effects (OLD)
# @software version: R=4.2.2


# Libraries ----
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(viridis)
library(reshape2)
library(ComplexHeatmap)
library(scales)

# GO enrichment --
library(WebGestaltR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(DOSE)


# ---- Data ----- ####
# Individual Traits ----
traits_cols <- c("Age" = "#56B4E9",
                 "Ancestry" = "#E69F00",
                 "Sex" =  "#009E73",
                 "BMI" = "#CC79A7")
traits <- names(traits_cols)

# Tissues ----
### read methylation results ####
first_dir <- "~/marenostrum/"
project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Age", "Ancestry", "BMI", "Sex")

tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Jose/00_Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

sex_tissues <- c('Ovary','Prostate','Testis')

tissues <- tissue_info$tissue_ID
# tissues cols
tissues_cols <- tissue_info$colcodes 
names(tissues_cols) <- tissue_info$tissue_abbrv

get_DMPs <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "SEX2"){
    NA
  }else{
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.rds"))[[trait]]
    #rownames(model[[trait]][model[[trait]]$adj.P.Val<0.05,])
    model
  }
}
DMPs_Res <- lapply(c('EURv1','SEX2','AGE','BMI'), function(trait) lapply(tissues, function(tissue) get_DMPs(tissue, trait)))
names(DMPs_Res) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(DMPs_Res[[trait]]) <- tissues}


# Lists of DEGs ----
get_DMPs <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "Sex"){
    NA
  }else{
    rownames(DMPs_Res[[trait]][[tissue]][DMPs_Res[[trait]][[tissue]]$adj.P.Val<0.05,])
  }
}
DMPs <- lapply(tissues, function(tissue) lapply(c("Ancestry", "Sex", "Age", "BMI"), function(trait) get_DMPs(tissue, trait)))
names(DMPs) <- tissues
for(trait in tissues){names(DMPs[[trait]]) <- c("Ancestry", "Sex", "Age", "BMI")}



# ---- Analysis ----- ####

# 1. Fraction of DEGs with 1, 2, 3 or 4 traits ----
# Data
d <- sapply(1:4, function(i)
  sapply(tissues, function(tissue) 
    sum(table(unlist(DMPs[[tissue]]))==i)
  ))
rownames(d) <- tissue_info$tissue_abbrv


# 1.1 Heatmap ----
library(wesanderson)
n_cols <- wes_palettes$GrandBudapest1[c(1,2,3,4)]

names(n_cols) <- c(1:4)
# Number of samples
metadata <- lapply(tissues, function(tissue) {
  metadata_ind <- readRDS(paste0(project_path, "Tissues/",tissue, "/metadata.rds"))
  print("Reading Admixture results")
  admixture_ancestry <- read.table('~/marenostrum_scratch/GTEx/v8/genotype_data/admixture_inferred_ancestry.txt')
  colnames(admixture_ancestry) <- c('SUBJID','AFRv1','EURv1','inferred_ancestry','AFRv2','EURv2')
  metadata_ind <- merge(metadata_ind, admixture_ancestry[,c("SUBJID","EURv1")], by='SUBJID')
  metadata_ind
})
names(metadata) <- tissues
n_samples <- sapply(tissues, function(tissue) nrow(mdata[[tissue]]))
names(n_samples) <- tissues
# Row annotation --
row_ha_left <- HeatmapAnnotation("Samples" = anno_barplot(n_samples,
                                                          gp = gpar(fill = tissue_info$colcodes,
                                                                    col = tissue_info$colcodes),
                                                          border=F),
                                 gap = unit(0.25,"cm"),
                                 show_legend = T, 
                                 show_annotation_name = T,
                                 annotation_name_rot = 90,
                                 annotation_name_gp = gpar(fontsize = 10),
                                 which = "row")
# Column annotation --
column_ha_top <- HeatmapAnnotation("Number of variables" = as.character(1:4),
                                   col = list("Number of variables" = n_cols),
                                   show_legend = T, show_annotation_name = F,
                                   simple_anno_size = unit(0.3,"cm"))
counts <- sapply(1:4, function(i)
  sapply(tissues, function(tissue) 
    sum(table(unlist(DMPs[[tissue]]))==i)
  ))
rownames(counts) <- tissue_info$tissue_abbrv
colnames(counts) <-   c("1 demographic trait",
                        "2 demographic traits",
                        "3 demographic traits",
                        "4 demographic traits")
# cell color % of tissue DEGs
head(counts)
ht <- Heatmap(apply(counts, 2,function(x) x/max(x,na.rm=T)),
              #ht <- Heatmap(100*t(apply(counts, 1,function(x) x/sum(x))),
              col = brewer.pal(9,"BuPu"),
              na_col = "white",
              cluster_rows = F,
              cluster_columns = F,
              name = "DE signal",
              row_names_side = "left",
              column_names_side = "top",
              column_names_rot =  45,
              top_annotation = column_ha_top,
              left_annotation = row_ha_left,
              row_names_gp = gpar(fontsize = 9),
              column_names_gp = gpar(fontsize = 9),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(prettyNum(counts[i, j], big.mark = ","), x, y, gp = gpar(fontsize = 10))},
              heatmap_legend_param = list(direction = "horizontal")
)
pdf("",
    width = 4, height = 9)
draw(ht,
     heatmap_legend_side = "bottom")
dev.off()

# Plot --
pdf("",
    width = 3, height = 8)
par(mfrow=c(1,1))
par(oma = c(0,4,0,0), mar = c(4,1,0,1))
barplot(apply(d, 1, function(x) x/sum(x))[,length(tissues):1],
        horiz = T, 
        las = 2,
        xaxt = 'n',
        cex.names = 0.9,
        xlab = "DEGs (%)",
        #col = brewer.pal(11, "Spectral")[c(8:11)],
        #col = brewer.pal(4, "GnBu"),#[c(3,5,8,9)],
        col = n_cols,
        border = NA)
axis(1, at =axTicks(1), labels = 100*axTicks(1))

dev.off()


# 01.Exprs.OverlapBetweenVariables.Heatmap 5x9

# 2. DEGs with 2 traits per tissue ----
pw.traits <- c("Age:Ancestry", "Age:Sex", "Age:BMI", 
               "Ancestry:Sex", "Ancestry:BMI",
               "Sex:BMI")
# List of DEGs per tissue
deg.pw.traits <- lapply(pw.traits, function(pw)
  lapply(tissues, function(tissue) 
    intersect(DMPs[[tissue]][[unlist(strsplit(pw, split = ":"))[[1]]]], DMPs[[tissue]][[unlist(strsplit(pw, split = ":"))[[2]]]])
  )
)
names(deg.pw.traits) <- pw.traits
for(pw in pw.traits){
  names(deg.pw.traits[[pw]]) <- tissues
}
# Number of DEGs
number_of_DEGs <- sapply(pw.traits, function(pw) 
  sapply(tissues, function(tissue)
    length(deg.pw.traits[[pw]][[tissue]])
  )
)
# Total number of DEGs with 2 traits (non unique)
total_DEGs <- sapply(pw.traits, function(pw) 
  length(unlist(deg.pw.traits[[pw]]))) 
pw.traits <- names(sort(total_DEGs, decreasing = T))
total_DEGs <- total_DEGs[pw.traits]

# Total number of DEGs with 2 traits ( unique)
total_DEGs_unique <- sapply(pw.traits, function(pw) 
  length(unique(unlist(deg.pw.traits[[pw]]))) )
#pw.traits <- names(sort(total_DEGs, decreasing = T))
total_DEGs_unique <- total_DEGs_unique[pw.traits]

# Tissue sharing ----
# Numbers
n_DEGs <- cbind.data.frame(sapply(c(1,2,3), function(i)
  sapply(pw.traits, function(pw) 
    sum(table(unlist(deg.pw.traits[[pw]]))==i)
  )),
  sapply(pw.traits, function(pw) 
    sum(table(unlist(deg.pw.traits[[pw]])) > 3)
  ))
n_DEGs <- n_DEGs[names(sort(total_DEGs,decreasing = T)),]
colnames(n_DEGs) <- c("1", "2", "3", ">3")

# Proportions
d <- cbind.data.frame(sapply(c(1,2,3), function(i)
  sapply(pw.traits, function(pw) 
    sum(table(unlist(deg.pw.traits[[pw]]))==i)
  )/total_DEGs_unique
),
sapply(pw.traits, function(pw) 
  sum(table(unlist(deg.pw.traits[[pw]])) > 3)
)/total_DEGs_unique
)
d <- d[names(sort(total_DEGs,decreasing = T)),]
colnames(d) <- c("1", "2", "3", ">3")
tissue_sharing <- d

pdf("",
    width = 8, height = 5)
#par(mfrow=c(1,2))
par(oma = c(4,0,0,0))
barplot(t(as.matrix(tissue_sharing)),
        horiz = F, 
        las = 2,
        cex.names = 0.9,
        ylab = "% of DEGs",
        col = brewer.pal(5, "Greys")[2:5],
        border = NA)
plot.new()
legend("center",
       c("N = 1",
         "N = 2",
         "N = 3",
         "N > 3"),
       pch = 15,
       col = brewer.pal(5, "Greys")[2:5],
       bty = 'n')
dev.off()


# 3. Main plots ----
# Counts: summary bar plot ----
df <- melt(number_of_DEGs)
colnames(df) <- c("Tissue","Variables","DEGs")
df$Variables <- factor(df$Variables, levels = pw.traits, order = T)
df$Tissue <- factor(df$Tissue, levels = tissue_info$tissue_ID, order = T)

# number of DEGs
# coloured by tissue, tissues ordered by sample size
p1 <- ggplot(df) +
  geom_bar(aes(x= Variables, y = DEGs,  fill = Tissue),#, order = DEGs),
           stat = "identity") +
  scale_fill_manual(values = tissue_info$colcodes) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0)) +
  xlab('') + ylab ("Number of DEGs") +
  scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE))
p1

# coloured by tissue, tissues ordered by contribution 
p1 <- ggplot(df) +
  geom_bar(aes(x= Variables, y = DEGs,  fill = Tissue),# group = DEGs),
           stat = "identity") +
  scale_fill_manual(values = tissue_info$colcodes) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0)) +
  xlab('') + ylab ("Number of DEGs") +
  scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE))
p1
