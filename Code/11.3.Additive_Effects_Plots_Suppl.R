#!/usr/bin/env Rscript

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

# Differential expression analyses: results tables ----
# for(tissue in tissues){
#   if(!file.exists(paste0("~/GTEx_v8/Raquel/03_DEA/01.DEA/Tissues/",tissue,"/", tissue,".voom_limma.covariates_and_traits.results.rds"))){print(tissue)}
# }
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
# n_cols <- c("#EA6B70",
#             "#CA3C70",
#             "#8F2390",
#             "#00299C")
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
# plot.new()
# legend("center",
#        c("1 demographic variable",
#          "2 demographic variables",
#          "3 demographic variables",
#          "4 demographic variables"),
#        pch = 15,
#        col = brewer.pal(5, "Greys")[2:5],
#        bty = 'n')
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
# % of DEGs 
# df <- melt( apply(number_of_DEGs, 2, function(x) x/sum(x)) )
# colnames(df) <- c("Tissue","Variables","DEGs")
# df$Variables <- factor(df$Variables, levels = pw.traits, order = T)
# df$Tissue <- factor(df$Tissue, levels = tissue_info$tissue_ID, order = T)
# p2 <- ggplot(df) +
#   geom_bar(aes(x= Variables, y = DEGs,  fill = Tissue),
#            stat = "identity") +
#   scale_fill_manual(values = tissue_info$colcodes) +
#   theme_minimal() +
#   theme(legend.position = "none",
#         axis.text.x = element_text(angle = 90, vjust = 0)) +
#   xlab('') + ylab ("% of DEGs")
# p2
# 
# p2 <- ggplot(df) +
#   geom_bar(aes(x= Variables, y = DEGs,  fill = Tissue, group = DEGs),
#            stat = "identity") +
#   scale_fill_manual(values = tissue_info$colcodes) +
#   theme_minimal() +
#   theme(legend.position = "none",
#         axis.text.x = element_text(angle = 90, vjust = 0)) +
#   xlab('') + ylab ("% of DEGs")
# p2

# 01.OverlapBetweenPairsOfVariables.Bar_plots 8x5
# grid.arrange(p1, p2, ncol =2)
# 
# # 'Grid' plot ----
# # Prepare data --
# # count data ----
# df <- melt(number_of_DEGs)
# colnames(df) <- c("Tissue","Variables","DEGs")
# df$Variables <- factor(df$Variables, levels = pw.traits, order = T)
# df$Tissue <- factor(df$Tissue, levels = tissue_info$tissue_ID, order = T)
# #df$DEGs <- 100*df$DEGs
# df$var1 <- sapply(as.character(df$Variables), function(i) unlist(strsplit(i, split = ":"))[[1]])
# df$var2 <- sapply(as.character(df$Variables), function(i) unlist(strsplit(i, split = ":"))[[2]])
# 
# # mock data ----
# mock.data <- rbind.data.frame(cbind.data.frame("Tissue" = tissue_info$tissue_ID,
#                                                "Variables" = "Age:Age",
#                                                "DEGs" = 0,
#                                                "var1" = "Age",
#                                                "var2" = "Age"),
#                               cbind.data.frame("Tissue" = tissue_info$tissue_ID,
#                                                "Variables" = "Ancestry:Ancestry",
#                                                "DEGs" = 0,
#                                                "var1" = "Ancestry",
#                                                "var2" = "Ancestry"),
#                               cbind.data.frame("Tissue" = tissue_info$tissue_ID,
#                                                "Variables" = "Sex:Sex",
#                                                "DEGs" = 0,
#                                                "var1" = "Sex",
#                                                "var2" = "Sex"
#                               ),
#                               cbind.data.frame("Tissue" = tissue_info$tissue_ID,
#                                                "Variables" = "BMI:BMI",
#                                                "DEGs" = 0,
#                                                "var1" = "BMI",
#                                                "var2" = "BMI"
#                               ),
#                               cbind.data.frame("Tissue" = tissue_info$tissue_ID,
#                                                "Variables" = "Sex:Ancestry",
#                                                "DEGs" = 0,
#                                                "var1" = "Sex",
#                                                "var2" = "Ancestry"
#                               ),
#                               cbind.data.frame("Tissue" = tissue_info$tissue_ID,
#                                                "Variables" = "BMI:Sex",
#                                                "DEGs" = 0,
#                                                "var1" = "BMI",
#                                                "var2" = "Sex"
#                               ),
#                               cbind.data.frame("Tissue" = tissue_info$tissue_ID,
#                                                "Variables" = "Sex:Age",
#                                                "DEGs" = 0,
#                                                "var1" = "Sex",
#                                                "var2" = "Age"
#                               ),
#                               cbind.data.frame("Tissue" = tissue_info$tissue_ID,
#                                                "Variables" = "Ancestry:Age",
#                                                "DEGs" = 0,
#                                                "var1" = "Ancestry",
#                                                "var2" = "Age"
#                               ),
#                               cbind.data.frame("Tissue" = tissue_info$tissue_ID,
#                                                "Variables" = "BMI:Ancestry",
#                                                "DEGs" = 0,
#                                                "var1" = "BMI",
#                                                "var2" = "Ancestry"
#                               ),
#                               cbind.data.frame("Tissue" = tissue_info$tissue_ID,
#                                                "Variables" = "BMI:Age",
#                                                "DEGs" = 0,
#                                                "var1" = "BMI",
#                                                "var2" = "Age"
#                               )
# )
# 
# # Add fdr to color bar border --
# fisher_results <- readRDS("~/GTEx_v8/Raquel/Draft/Analysis/Expression.OverlapBetweenTraits/Data/03.OverlapBetweenVariables.Fisher_results.rds")
# # black border if ENRICHMENT
# df$fdr <- sapply(1:nrow(df), function(row) ifelse(is.na(fisher_results[fisher_results$Tissue_ID == df[row,"Tissue"] & 
#                                                                          fisher_results$Traits == df[row, "Variables"], "FDR"]),
#                                                   "n.s",
#                                                   ifelse(fisher_results[fisher_results$Tissue_ID == df[row,"Tissue"] & 
#                                                                           fisher_results$Traits == df[row, "Variables"], "FDR"] < 0,
#                                                          "n.s",
#                                                          "s")))
# mock.data$fdr <- "n.s"
# data <- rbind.data.frame(df, mock.data)                              
# data$var1 <- factor(data$var1, levels = traits, order = T)
# data$var2 <- factor(data$var2, levels = traits, order = T)
# data$Tissue <- factor(data$Tissue, levels = tissue_info$tissue_ID, order = T)
# # Add x-dummy value
# data$x <- 1
# # Add border color if FDR < 0.05
# data$fdr <- factor(data$fdr, levels = c("s", "n.s"), order = T)
# # 
# data$col <- sapply(1:nrow(data), function(row) 
#   ifelse(data[row,"fdr"]=="n.s",
#          alpha(tissue_info[tissue_info$tissue_ID== data[row, "Tissue"], "colcodes"],0.25),
#          tissue_info[tissue_info$tissue_ID== data[row, "Tissue"], "colcodes"]
#   ))
# 
# # 01.OverlapBetweenVariables.Grid.Number_DEGs.pdf 9x3
# ggplot(data) + 
#   geom_bar(aes(x = x, y = DEGs,  fill = Tissue, group = DEGs),
#            stat = "identity") +
#   scale_fill_manual(values = tissue_info$colcodes) +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   xlab('') + ylab ("% of DEGs") +
#   facet_grid(var1 ~ var2) +
#   coord_flip() +
#   scale_y_continuous(position = "right") +
#   scale_x_continuous(position = "top") +
#   scale_colour_manual(values = c("black", NA))
# 
# # 'Grid' plot
# # 01.OverlapBetweenVariables.Grid_pot.Border (9 x 3)
# # ggplot(data) + 
# #   geom_bar(aes(x = x, y = DEGs,  fill = Tissue, col = fdr, group = DEGs),
# #            stat = "identity") +
# #   scale_fill_manual(values = tissue_info$colcodes) +
# #   theme_minimal() +
# #   theme(legend.position = "none") +
# #   xlab('') + ylab ("Number of DEGs") +
# #   facet_grid(var1 ~ var2) +
# #   coord_flip() +
# #   scale_y_continuous(position = "right") +
# #   scale_colour_manual(values = c("black", NA))
# 
# # Prepare data --
# # Grid sharing data ----
# df <- melt(d)
# colnames(df) <- c("Sharing","value")
# df$Variables <- rep(rownames(d), 4)
# df$var2 <- sapply(as.character(df$Variables), function(i) unlist(strsplit(i, split = ":"))[[1]])
# df$var1 <- sapply(as.character(df$Variables), function(i) unlist(strsplit(i, split = ":"))[[2]])
# df$Variables <- paste0(df$var1, ":", df$var2)
# 
# # mock data ----
# mock.data <- rbind.data.frame(cbind.data.frame("Sharing" = rep(1,6),
#                                                "value" = 0,
#                                                "Variables" = "Age:Age",
#                                                "var1" = "Age",
#                                                "var2" = "Age"),
#                               cbind.data.frame("Sharing" = rep(1,6),
#                                                "value" = 0,
#                                                "Variables" = "Ancestry:Ancestry",
#                                                "var1" = "Ancestry",
#                                                "var2" = "Ancestry"),
#                               cbind.data.frame("Sharing" = rep(1,6),
#                                                "value" = 0,
#                                                "Variables" = "Sex:Sex",
#                                                "var1" = "Sex",
#                                                "var2" = "Sex"),
#                               cbind.data.frame("Sharing" = rep(1,6),
#                                                "value" = 0,
#                                                "Variables" = "BMI:BMI",
#                                                "var1" = "BMI",
#                                                "var2" = "BMI"),
#                               cbind.data.frame("Sharing" = rep(1,6),
#                                                "value" = 0,
#                                                "Variables" = "Age:Ancestry",
#                                                "var1" = "Age",
#                                                "var2" = "Ancestry"),
#                               cbind.data.frame("Sharing" = rep(1,6),
#                                                "value" = 0,
#                                                "Variables" = "Age:Sex",
#                                                "var1" = "Age",
#                                                "var2" = "Sex"),
#                               cbind.data.frame("Sharing" = rep(1,6),
#                                                "value" = 0,
#                                                "Variables" = "Age:BMI",
#                                                "var1" = "Age",
#                                                "var2" = "BMI"),
#                               cbind.data.frame("Sharing" = rep(1,6),
#                                                "value" = 0,
#                                                "Variables" = "Ancestry:Sex",
#                                                "var1" = "Ancestry",
#                                                "var2" = "Sex"),
#                               cbind.data.frame("Sharing" = rep(1,6),
#                                                "value" = 0,
#                                                "Variables" = "Ancestry:BMI",
#                                                "var1" = "Ancestry",
#                                                "var2" = "BMI"),
#                               cbind.data.frame("Sharing" = rep(1,6),
#                                                "value" = 0,
#                                                "Variables" = "Sex:BMI",
#                                                "var1" = "Sex",
#                                                "var2" = "BMI")
#                               
# )
# 
# data <- rbind.data.frame(df, mock.data)                              
# data$var1 <- factor(data$var1, levels = traits, order = T)
# data$var2 <- factor(data$var2, levels = traits, order = T)
# data$Sharing <- factor(data$Sharing, levels = rev(c("1", "2", "3", ">3")), order = T)
# # Add x-dummy value
# data$x <- 1
# # 01.OverlapBetweenVariables.Grid.Number_DEGs.pdf 9x3
# ggplot(data) + 
#   geom_bar(aes(x = x, y = value,  fill = Sharing),
#            stat = "identity") +
#   scale_fill_manual(values = rev(brewer.pal(4, "Greys"))) +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   xlab('') + ylab ("% of DEGs") +
#   facet_grid(var1 ~ var2) +
#   coord_flip() +
#   scale_y_continuous(position = "left") +
#   scale_x_continuous(position = "bottom") 

# 01.OverlapBetweenVariables.Grid.Tissue_sharing.pdf 9x3