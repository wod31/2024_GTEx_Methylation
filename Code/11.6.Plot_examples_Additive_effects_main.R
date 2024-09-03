#!/usr/bin/env Rscript
set.seed(1)

# Libraries ----
library(ggplot2)
library(ggtext)
library(ggpubr)

# function to plot N in box plots ---
get_box_stats <- function(y, upper_limit = max(y) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "N =", length(y), "\n"#,
      #"Mean =", round(mean(y), 2), "\n",
      #"Median =", round(median(y), 2), "\n"
    )
  ))
}

# ---- Data ---- ####

# Individual Traits ----
traits_cols <- c("Age" = "#56B4E9",
                 "Ancestry" = "#E69F00",
                 "Sex" = "#009E73",
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

# Metadata ----
metadata <- lapply(tissues, function(tissue) {
  metadata_ind <- readRDS(paste0(project_path, "Tissues/",tissue, "/metadata.rds"))
  print("Reading Admixture results")
  admixture_ancestry <- read.table('~/marenostrum_scratch/GTEx/v8/genotype_data/admixture_inferred_ancestry.txt')
  colnames(admixture_ancestry) <- c('SUBJID','AFRv1','EURv1','inferred_ancestry','AFRv2','EURv2')
  metadata_ind <- merge(metadata_ind, admixture_ancestry[,c("SUBJID","EURv1")], by='SUBJID')
  metadata_ind
})
names(metadata) <- tissues

for(tissue in sex_tissues[c(1)]){
  metadata[[tissue]]$SEX <- "2"
  metadata[[tissue]]$SEX <- as.factor(metadata[[tissue]]$SEX)
  metadata[[tissue]] <- metadata[[tissue]][, colnames(metadata$MuscleSkeletal)[c(1:8,10:21)]]
}
for(tissue in sex_tissues[c(2,3)]){
  metadata[[tissue]]$SEX <- "1"
  metadata[[tissue]]$SEX <- as.factor(metadata[[tissue]]$SEX)
  metadata[[tissue]] <- metadata[[tissue]][, colnames(metadata$MuscleSkeletal)[c(1:8,10:21)]]
}


# 1. Hier.part:expression ----
# All genes expressed in tissue 
hier.part.exprs <- lapply(tissues, function(tissue)
  readRDS(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/",tissue,"/","hier.part.5peers.continous.rds")))
names(hier.part.exprs) <- tissues

# 2. Differential expression results ----
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

# Functions ----
pseudo_categorize_age <- function(age){
  if(age < 30){
    return("[20-30)")
  }else if(age < 40){
    return("[30-40)")
  }else if(age < 50){
    return("[40-50)")
  }else if(age < 60){
    return("[50-60)")
  }else{
    return("[60-70]")
  }
}
pseudo_categorize_bmi <- function(bmi){
  if(bmi < 25){
    return("Normal")
  }else if(bmi < 30){
    return("Overweight")
  }else{
    return("Obese")
  }
}

pseudo_categorize_ancestry <- function(ancestry){
  if(ancestry < 0.25){
    return("AFR")
  }else if(ancestry > 0.75){
    return("EUR")
  }else{
    return("ADX")
  }
}

get.gene.data <- function(gene_name){
  #  Expresidualssion data ----
  # Gene TPM ----
  df <- cbind.data.frame(t(readRDS(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/", tissue, "/", "data.rds"))[gene,]))
  
  # Gene residuals ----
  #residuals <- readRDS(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/", tissue, "/", "methylation_residuals.continous.rds"))
  
  #identical(colnames(tpm), colnames(res))
  
  #df <- cbind.data.frame(t(tpm[gene,]))
  colnames(df) <- "Beta"
  df$residuals<- readRDS(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/", tissue, "/", "methylation_residuals.continous.rds"))[gene,]
  df$Age_int <- sapply(rownames(df), function(i) metadata[[tissue]][metadata[[tissue]]$SUBJID==i,"AGE"])
  df$BMI_int <- sapply(rownames(df), function(i) metadata[[tissue]][metadata[[tissue]]$SUBJID==i,"BMI"])
  #df$Age <- sapply(df$Age_int, function(i) ifelse(i<45, "[20-45)", "[45-70]"))
  #df$Age <- factor(df$Age, 
  #                 levels = c("[20-45)", "[45-70]"),
  #                 order = T)
  df$Age <- sapply(df$Age_int, function(age) pseudo_categorize_age(age))
  df$Age <- factor(df$Age, 
                   levels = c("[20-30)",
                              "[30-40)",
                              "[40-50)",
                              "[50-60)",
                              "[60-70]"),
                   order = T)
  df$Ancestry_int <- sapply(rownames(df), function(i) metadata[[tissue]][metadata[[tissue]]$SUBJID==i,"EURv1"])
  df$Ancestry <- sapply(df$Ancestry_int, function(ancestry) pseudo_categorize_ancestry(ancestry))
  df$Ancestry <- factor(df$Ancestry, 
                        levels = c('AFR','ADX','EUR'))

  df$Sex <- sapply(rownames(df), function(i) metadata[[tissue]][metadata[[tissue]]$SUBJID==i,"SEX"])
  df$Sex <- gsub("1", "Male", df$Sex)
  df$Sex <- gsub("2", "Female", df$Sex)
  df$Sex <- factor(df$Sex, levels = c("Male", "Female"), order = T)
  df$BMI <- sapply(df$BMI_int, function(bmi) pseudo_categorize_bmi(bmi))
  df$BMI <- factor(df$BMI,
                   levels = c("Normal", "Overweight", "Obese"),
                   order = T)
  df$Tissue <- rep(tissue_info[tissue_info$tissue_ID==tissue, "tissue_abbrv"], nrow(df))
  df$Gene <- rep(gene, nrow(df))
  head(df)
  return(df)
}

# ColonTransverse - CDH1 cg05356648 Age & Ancestry // cg21473782 AATK // cg14188106 TNXB----
tissue <- "ColonTransverse"
gene_name <- "cg14188106" 
gene <- gene_name
#sapply(traits, function(trait) gene %in% DMPs[[tissue]][[trait]])

# # Gene TPM ----
# tpm <- readRDS(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/", tissue, "/", "data.rds"))
# 
# # Gene residuals ----
# residuals <- readRDS(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/", tissue, "/", "methylation_residuals.continous.rds"))

# get gene data ----
data <- get.gene.data(gene_name)
head(data)

# log2(TPM) ----
# option 1 ----
# plot 1 ----
cols <- c(traits_cols["Age"], traits_cols["Age"],traits_cols["Age"],traits_cols["Age"],traits_cols["Age"])
names(cols) <- c("[20-30)",
                 "[30-40)",
                 "[40-50)",
                 "[50-60)",
                 "[60-70]")
p1 <- ggplot(data = data,
             aes(x = Age,
                 y = as.numeric(Beta),
                 fill = Age),
) +
  geom_violin(col = "black") +
  geom_boxplot(col = "black",
               fill = "white",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  geom_jitter(col = "black", 
              alpha = 0.1,
              size = 0.8) +
  xlab("") +
  ylab("Beta") +
  scale_fill_manual(values = cols) +
  stat_summary(fun.data = get_box_stats, geom = "text",
               hjust = 0.5, vjust = 0.9, size = 2) +
  labs(title=gene_name) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.ticks = element_line(size = 0.1),
        axis.text.x = element_text(size = 10),
        axis.title.y = ggtext::element_markdown(size=10),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5,
                                  size = 10),
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none") 

cols <- c(traits_cols["Ancestry"], traits_cols["Ancestry"], traits_cols["Ancestry"])
names(cols) <- c("AFR", "ADX", "EUR")
p2 <- ggplot(data = data,
             aes(x = Ancestry,
                 y = as.numeric(Beta),
                 fill = Ancestry),
) +
  geom_violin(col = "black") +
  geom_boxplot(col = "black",
               fill = "white",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  geom_jitter(col = "black", 
              alpha = 0.1,
              size = 0.8) +
  xlab("") +
  ylab("Beta") +
  scale_fill_manual(values = cols) +
  stat_summary(fun.data = get_box_stats, geom = "text",
               hjust = 0.5, vjust = 0.9, size = 2) +
  labs(title=gene_name) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.ticks = element_line(size = 0.1),
        axis.text.x = element_text(size = 10),
        axis.title.y = ggtext::element_markdown(size=10),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5,
                                  size = 10),
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none") 
p2
# plot 2 ----
d <- data.frame("variable" = c("Ancestry", "Age"),
                "value"  = as.numeric(hier.part.exprs[[tissue]][gene_name, c("EURv1_abs", "AGE_abs")]))
variables_col <- c(traits_cols["Ancestry"], traits_cols["Age"])
names(variables_col) <- c("Ancestry", "Age")
d$variable <- factor(d$variable, levels = rev(c("Ancestry", "Age")), order = T)
p3 <- ggplot(d, aes(x=1, y=value, fill = variable)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = variables_col) + 
  coord_flip() +
  xlab("") + ylab("Expression variation explained (%)") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.ticks.x =  element_line(size = 0.5),
        axis.ticks.y =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = ggtext::element_markdown(size=10),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5,
                                  size = 10),
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none") 
p3
pdf("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/cg14188106_Colon_additive.pdf",
    width = 3, height = 6)
ggarrange(p1, p2, p3, nrow = 3, heights =  c(4,4,1))
dev.off()

#  Expression residuals ----
# option 1 ----
# plot 1 ----
cols <- c(traits_cols["Age"], traits_cols["Age"],traits_cols["Age"],traits_cols["Age"],traits_cols["Age"])
names(cols) <- c("[20-30)",
                 "[30-40)",
                 "[40-50)",
                 "[50-60)",
                 "[60-70]")
p1 <- ggplot(data = data,
             aes(x = Age,
                 y = residuals,
                 fill = Age),
) +
  geom_violin(col = "black") +
  geom_boxplot(col = "black",
               fill = "white",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  geom_jitter(col = "black", 
              alpha = 0.1,
              size = 0.8) +
  xlab("") +
  ylab("Methylation residuals") +
  scale_fill_manual(values = cols) +
  stat_summary(fun.data = get_box_stats, geom = "text",
               hjust = 0.5, vjust = 0.9, size = 2) +
  labs(title=gene_name) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.ticks = element_line(size = 0.1),
        axis.text.x = element_text(size = 10),
        axis.title.y = ggtext::element_markdown(size=10),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5,
                                  size = 10),
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none") 

cols <- c(traits_cols["Ancestry"], traits_cols["Ancestry"], traits_cols["Ancestry"])
names(cols) <- c("AFR", "ADX", "EUR")
p2 <- ggplot(data = data,
             aes(x = Ancestry,
                 y = residuals,
                 fill = Ancestry),
) +
  geom_violin(col = "black") +
  geom_boxplot(col = "black",
               fill = "white",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  geom_jitter(col = "black", 
              alpha = 0.1,
              size = 0.8) +
  xlab("") +
  ylab("Methylation residuals") +
  scale_fill_manual(values = cols) +
  stat_summary(fun.data = get_box_stats, geom = "text",
               hjust = 0.5, vjust = 0.9, size = 2) +
  labs(title=gene_name) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.ticks = element_line(size = 0.1),
        axis.text.x = element_text(size = 10),
        axis.title.y = ggtext::element_markdown(size=10),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5,
                                  size = 10),
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none") 
p2
# plot 2 ----
d <- data.frame("variable" = c("Ancestry", "Age"),
                "value"  = as.numeric(hier.part.exprs[[tissue]][gene_name, c("EURv1_abs", "AGE_abs")]))
variables_col <- c(traits_cols["Ancestry"], traits_cols["Age"])
names(variables_col) <- c("Ancestry", "Age")
d$variable <- factor(d$variable, levels = rev(c("Ancestry", "Age")), order = T)
p3 <- ggplot(d, aes(x=1, y=value, fill = variable)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = variables_col) + 
  coord_flip() +
  xlab("") + ylab("Expression variation explained (%)") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.ticks.x =  element_line(size = 0.5),
        axis.ticks.y =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = ggtext::element_markdown(size=10),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5,
                                  size = 10),
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none") 
p3
# pdf("~/GTEx_v8/Raquel/Draft/Analysis/Expression.OverlapBetweenTraits/figures/04.Sex_BMI.AdiposeSubcutaneous.EGFL6.option1.Expression_residuals.pdf",
#     width = 3, height = 8)
ggarrange(p1, p2, p3, nrow = 3, heights =  c(4,4,1))
#dev.off()

# option 2 ----
# cols <- c(traits_cols["Age"], traits_cols["Age"],traits_cols["Age"],traits_cols["Age"],traits_cols["Age"])
# names(cols) <- c("[20-30)",
#                  "[30-40)",
#                  "[40-50)",
#                  "[50-60)",
#                  "[60-70]")
# p4 <- ggplot(data = data,
#              aes(x = Age,
#                  y = residuals,
#                  fill = Age),
# ) +
#   geom_violin(col = "black") +
#   geom_boxplot(col = "black",
#                fill = "white",
#                outlier.shape = NA,
#                notch = T,
#                width = 0.25) +
#   geom_jitter(col = "black", 
#               alpha = 0.1,
#               size = 0.8) +
#   xlab("") +
#   ylab("Expression residuals") +
#   scale_fill_manual(values = cols) +
#   facet_grid(~Ancestry) +
#   stat_summary(fun.data = get_box_stats, geom = "text",
#                hjust = 0.5, vjust = 0.9, size = 2) +
#   labs(title=gene_name) +
#   theme(panel.background = element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#         axis.ticks = element_line(size = 0.1),
#         axis.text.x = element_text(size = 10),
#         axis.title.y = ggtext::element_markdown(size=10),
#         axis.title = element_text(size = 10),
#         plot.title = element_text(hjust = 0.5,
#                                   size = 10),
#         strip.background = element_rect(fill="#CC79A7"),
#         legend.position = "none") 
# 
# cols <- c(traits_cols["BMI"], traits_cols["BMI"], traits_cols["BMI"])
# names(cols) <- c("Normal", "Overweight", "Obese")
# p5 <- ggplot(data = data,
#              aes(x = BMI,
#                  y = residuals,
#                  fill = BMI),
# ) +
#   geom_violin(col = "black") +
#   geom_boxplot(col = "black",
#                fill = "white",
#                outlier.shape = NA,
#                notch = T,
#                width = 0.25) +
#   geom_jitter(col = "black", 
#               alpha = 0.1,
#               size = 0.8) +
#   xlab("") +
#   ylab("Expression residuals") +
#   scale_fill_manual(values = cols) +
#   facet_grid(~Sex) +
#   stat_summary(fun.data = get_box_stats, geom = "text",
#                hjust = 0.5, vjust = 0.9, size = 2) +
#   labs(title=gene_name) +
#   theme(panel.background = element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#         axis.ticks = element_line(size = 0.1),
#         axis.text.x = element_text(size = 10),
#         axis.title.y = ggtext::element_markdown(size=10),
#         axis.title = element_text(size = 10),
#         plot.title = element_text(hjust = 0.5,
#                                   size = 10),
#         strip.background = element_rect(fill="#009E73"),
#         legend.position = "none") 
# 
# pdf("~/GTEx_v8/Raquel/Draft/Analysis/Expression.OverlapBetweenTraits/figures/04.Sex_BMI.AdiposeSubcutaneous.EGFL6.option2.Expression_residuals.pdf",
#     width = 5, height = 6)
# ggarrange(p4, p5, p3, nrow = 3, heights =  c(4,4,1))
# dev.off()

