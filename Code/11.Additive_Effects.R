#!/usr/bin/env Rscript
# @Author: Winona Oliveros
# @E-mail: winona.oliveros@bsc.es
# @Description: Code to get additive effects between demographic traits
# @software version: R=4.2.2

##########################################
#### Additive effects ####################
##########################################

#!/usr/bin/env Rscript

# Libraries ----
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)
library(circlize)

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

# Differential methylation analyses: results tables ----

get_DMPs <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "SEX2"){
    NA
  }else{
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/DML_results_5_PEERs_continous.rds"))[[trait]]
    model
  }
}
DMPs_Res <- lapply(c('EURv1','SEX2','AGE','BMI'), function(trait) lapply(tissues, function(tissue) get_DMPs(tissue, trait)))
names(DMPs_Res) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(DMPs_Res[[trait]]) <- tissues}


# Lists of DMPs ----
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

# 1. Overlap between pairwise combinations of traits ----

tested_cpgs <- lapply(tissues, function(tissue) rownames(DMPs_Res[['Age']][[tissue]]))
names(tested_cpgs) <- tissues

# Pairwise combination of traits --
pw.traits <- c("Age:Ancestry",
               "Ancestry:Sex",
               "Age:Sex",
               "Sex:BMI",
               "Ancestry:BMI",
               "Age:BMI"
               
)

# List of DMPs per tissue
deg.pw.traits <- lapply(pw.traits, function(pw)
  lapply(tissues, function(tissue) 
    intersect(DMPs[[tissue]][[unlist(strsplit(pw, split = ":"))[[1]]]], DMPs[[tissue]][[unlist(strsplit(pw, split = ":"))[[2]]]])
  )
)
names(deg.pw.traits) <- pw.traits
for(pw in pw.traits){
  names(deg.pw.traits[[pw]]) <- tissues
}

length(unique(unlist(lapply(names(deg.pw.traits), function(pw) unique(unlist(deg.pw.traits[[pw]]))))))

# Number of DMPs
number_of_DEGs <- sapply(pw.traits, function(pw) 
  sapply(tissues, function(tissue)
    length(deg.pw.traits[[pw]][[tissue]])
  )
)

total_unique_pw_DEGs <- sapply(pw.traits, function(pw) 
  length(unique(unlist(lapply(tissues, function(tissue) deg.pw.traits[[pw]][[tissue]]))))
)
total_pw_DEGs <- sapply(pw.traits, function(pw) 
  length(unlist(lapply(tissues, function(tissue) deg.pw.traits[[pw]][[tissue]])))
)
total_unique_pw_DEGs
total_pw_DEGs

# Number of tissues with genes with additive effects ----
apply(number_of_DEGs, 2, function(x) sum(x>0, na.rm=T))


rownames(number_of_DEGs) <- tissue_info$tissue_abbrv
tissues_cols <- tissue_info$colcodes 
names(tissues_cols) <- tissue_info$tissue_abbrv
row_ha_left <- HeatmapAnnotation("Tissues" = rownames(number_of_DEGs),
                                 col = list("Tissues" = tissues_cols),
                                 show_legend = F, 
                                 show_annotation_name = F,
                                 simple_anno_size = unit(0.3,"cm"),
                                 which = "row")

colnames(number_of_DEGs) <- c("Ancestry-age DMPs",
                              "Ancestry-sex DMPs",
                              "Sex-age DMPs",
                              "Sex-BMI DMPs",
                              "Ancestry-BMI DMPs",
                              "Age-BMI DMPs")
column_ha_top <- HeatmapAnnotation("Number of DMPs" = anno_barplot(total_unique_pw_DEGs,
                                                                   gp = gpar(fill = "light gray",
                                                                             col = "light gray"),
                                                                   border = F),
                                   annotation_name_side = "left",
                                   annotation_name_rot = 90,
                                   annotation_name_gp = gpar(fontsize = 8),
                                   height = unit(2.5, "cm"),
                                   gap = unit(0.25,"cm"))
column_ha_bottom <- HeatmapAnnotation("variable1" = sapply(pw.traits, function(i) unlist(strsplit(i, split = ":"))[[1]]),
                                      "variable2" = sapply(pw.traits, function(i) unlist(strsplit(i, split = ":"))[[2]]),
                                      col = list("variable1" = traits_cols,
                                                 "variable2" = traits_cols),
                                      show_legend = F, 
                                      show_annotation_name = F,
                                      simple_anno_size = unit(0.3,"cm"))


ht <- Heatmap(apply(number_of_DEGs, 2, function(x) x/sum(x)),
              col = brewer.pal(9, "BuPu"),
              na_col = "white",
              cluster_rows = F,
              cluster_columns = F,
              name = "DE signal",
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(prettyNum(number_of_DEGs[i, j], big.mark = ","), x, y, gp = gpar(fontsize = 10))},
              row_names_side = "left",
              column_names_side = "bottom",
              left_annotation = row_ha_left,
              row_names_gp = gpar(fontsize = 10),
              bottom_annotation = column_ha_bottom,
              column_names_gp = gpar(fontsize = 10),
              top_annotation = column_ha_top,
              heatmap_legend_param = list(direction = "horizontal"))


# Fisher's exact test to test if there are more DMPs with 2 traits than expected ----
overlap.enrichment.fun <- function(tissue, trait.1, trait.2){
  if(tissue %in% sex_tissues & (trait.1 == "Sex" | trait.2 == "Sex")){
    return(list("overlap" = NA,
                "odds.ratio" = NA,
                "p.value" = NA))
  }else{
 
    x11 <- length(intersect(DMPs[[tissue]][[trait.1]], DMPs[[tissue]][[trait.2]]))
    x12 <- length(DMPs[[tissue]][[trait.1]][! DMPs[[tissue]][[trait.1]] %in% DMPs[[tissue]][[trait.2]]])
    x21 <- length(DMPs[[tissue]][[trait.2]][! DMPs[[tissue]][[trait.2]] %in% DMPs[[tissue]][[trait.1]]])
    x22 <- length(tested_cpgs[[tissue]][! tested_cpgs[[tissue]] %in% unique(c(DMPs[[tissue]][[trait.1]], DMPs[[tissue]][[trait.2]]))])
    m <- matrix(c(x11,x12,x21,x22),2,2,byrow = T)
    rownames(m) <- c(paste0(trait.1,".deg"), paste0(trait.1, ".not_deg"))
    colnames(m) <- c(paste0(trait.2,".deg"), paste0(trait.2, ".not_deg"))
    #sum(m) == length(exprs.genes[[tissue]])
    f <- fisher.test(m)
    f$estimate
    return(list("overlap" = x11,
                "odds.ratio" = f$estimate,
                "p.value" = f$p.value,
                "counts.matrix" = m,
                "lower_CI" = f$conf.int[1],
                "upper_CI" = f$conf.int[2]))
  }
}


# Enrichment analysis  --
fisher.results <- lapply(pw.traits, function(pw)
  lapply(tissues, function(tissue) 
    overlap.enrichment.fun(tissue, unlist(strsplit(pw, split = ":"))[[1]], unlist(strsplit(pw, split = ":"))[[2]])
  ))
names(fisher.results) <- pw.traits
for(pw in pw.traits){names(fisher.results[[pw]]) <- tissues}

# Enrichment statistics --
p.value <- sapply(pw.traits, function(pw)
  sapply(tissues, function(tissue) 
    fisher.results[[pw]][[tissue]][["p.value"]]
  )
)

# Multiple testing correction by tissue --
fdr <- t(apply(p.value, 1, function(x) p.adjust(x, method = "BH")))
fdr[fdr>=0.05] <- NA

# multiple testing correction across tissues and number of pairwise combinations of traits
adj_p_values <- p.adjust(unlist(p.value), method = "BH")
adjPVal_matrix <- matrix(adj_p_values, 
                         nrow = length(tissues), ncol = length(pw.traits),
                         byrow = F)

# Number of DMPs in overlap --
overlap <- sapply(pw.traits, function(pw)
  sapply(tissues, function(tissue) 
    fisher.results[[pw]][[tissue]][["overlap"]]
  )
)
rownames(overlap) <- tissues

# Odds ratio --
odds.ratio <- sapply(pw.traits, function(pw)
  sapply(tissues, function(tissue) 
    fisher.results[[pw]][[tissue]][["odds.ratio"]]
  )
)
rownames(odds.ratio) <- tissues
or <- odds.ratio
or[is.na(fdr)] <- NA

log2_OR <- log2(or)
range(log2_OR, na.rm = T)
col_fun <- colorRamp2(seq(0,3,length.out=9),
                      (brewer.pal(9, "Reds")))


colnames(log2_OR) <- c("Ancestry-age DMPs",
                       "Ancestry-sex DMPs",
                       "Sex-age DMPs",
                       "Sex-BMI DMPs",
                       "Ancestry-BMI DMPs",
                       "Age-BMI DMPs")
rownames(log2_OR) <- tissue_info$tissue_abbrv
tissues_cols <- tissue_info$colcodes 
names(tissues_cols) <- tissue_info$tissue_abbrv
row_ha_left <- HeatmapAnnotation("Tissues" = rownames(log2_OR),
                                 col = list("Tissues" = tissues_cols),
                                 show_legend = F, 
                                 show_annotation_name = F,
                                 simple_anno_size = unit(0.3,"cm"),
                                 which = "row")

m <- as.matrix(rbind(apply(or, 2, function(x) sum(x<1, na.rm=T)),
                     apply(or, 2, function(x) sum(x>1, na.rm=T))))

column_ha_top_ht2 <- HeatmapAnnotation("Number of tissues" = anno_barplot(t(m),
                                                                          gp = gpar(fill = brewer.pal(11, "RdBu")[c(10,2)],
                                                                                    col = brewer.pal(11, "RdBu")[c(10,2)]),
                                                                          border = F),
                                       annotation_name_side = "right",
                                       annotation_name_rot = 90,
                                       annotation_name_gp = gpar(fontsize = 8),
                                       height = unit(2.5, "cm"),
                                       gap = unit(0.25,"cm"))

pw.traits[1] <- "Ancestry:Age"
pw.traits[3] <- "Sex:Age"
column_ha_bottom <- HeatmapAnnotation("variable1" = sapply(pw.traits, function(i) unlist(strsplit(i, split = ":"))[[1]]),
                                      "variable2" = sapply(pw.traits, function(i) unlist(strsplit(i, split = ":"))[[2]]),
                                      col = list("variable1" = traits_cols,
                                                 "variable2" = traits_cols),
                                      show_legend = F, 
                                      show_annotation_name = F,
                                      simple_anno_size = unit(0.3,"cm"))

ht2 <- Heatmap(log2_OR,
               col = col_fun,
               na_col = "white",
               cluster_rows = F,
               cluster_columns = F,
               name = "Odds ratio (log2)",
               #cell_fun = function(j, i, x, y, width, height, fill) {
               #  grid.text(round(-log10(fdr[i, j]), 2), x, y, gp = gpar(fontsize = 10))},
               row_names_side = "left",
               column_names_side = "bottom",
               #left_annotation = row_ha_left,
               row_names_gp = gpar(fontsize = 10),
               bottom_annotation = column_ha_bottom,
               column_names_gp = gpar(fontsize = 10),
               top_annotation = column_ha_top_ht2,
               heatmap_legend_param = list(direction = "horizontal"))


log10_fdr_sign <- -log10(fdr)
for(i in 1:nrow(log10_fdr_sign)){
  for(j in 1:ncol(log10_fdr_sign)){
    if(!is.na(log2_OR[i,j]) & log2_OR[i,j]<0){
      log10_fdr_sign[i, j] <- -log10_fdr_sign[i, j]
    }
  }
}

col_fun <- colorRamp2(seq(0,70,length.out=9),
                      (brewer.pal(9, "Reds")))
rownames(log10_fdr_sign) <- tissue_info$tissue_abbrv
ht <- Heatmap(log10_fdr_sign,
              col = col_fun,
              na_col = "white",
              cluster_rows = F,
              cluster_columns = F,
              name = "-log10(FDR)",
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(round(log2_OR[i, j], 2), x, y, gp = gpar(fontsize = 10))},
              row_names_side = "left",
              column_names_side = "bottom",
              left_annotation = row_ha_left,
              row_names_gp = gpar(fontsize = 10),
              bottom_annotation = column_ha_bottom,
              column_names_gp = gpar(fontsize = 10),
              top_annotation = column_ha_top,
              heatmap_legend_param = list(direction = "horizontal"))


#######################################

# Build summary data.table --
# Tissue, pairwise:traits, overlap, odds.ratio, FDR 
data <- do.call(rbind.data.frame,
                lapply(tissues, function(tissue) 
                  data.frame("Tissue" = rep(tissue_info[tissue_info$tissue_ID==tissue, "tissue_abbrv"], length(pw.traits)),
                             "Tissue_ID" = rep(tissue,  length(pw.traits)),
                             "Traits" = pw.traits,
                             "Overlap" = overlap[tissue,], # all tissue and pairwise trait combinations
                             "OR" = or[tissue, ], # only fdr < 0.05
                             "fdr" = fdr[tissue,] # only fdr < 0.05
                  )))
data$FDR <- -log10(data$fdr)
data$FDR <- ifelse(data$OR < 1, (-1)*data$FDR, data$FDR) # if depletion, blue (negative)
data$OddsRatio <- log2(data$OR)
data$Overlap <- ifelse(data$Overlap=="0","",data$Overlap)
data$Overlap <- as.numeric(data$Overlap)
data$Tissue <- factor(data$Tissue, levels = rev(tissue_info$tissue_abbrv), order = T)
data <- data[!is.na(data$fdr),]
data$lower_CI <- sapply(1:nrow(data), function(row) 
  fisher.results[[ data[row, "Traits"] ]][[ data[row, "Tissue_ID"] ]][['lower_CI']])
data$upper_CI <- sapply(1:nrow(data), function(row) 
  fisher.results[[ data[row, "Traits"] ]][[ data[row, "Tissue_ID"] ]][['upper_CI']])
data <- data[data$OR>1 & data$Overlap >= 10,]
# Limit color scale (FDR) --
data$FDR_bounded <- ifelse(is.na(data$FDR), NA,
                           ifelse(data$FDR < -5, -5,
                                  ifelse(data$FDR > 5, 5, data$FDR))) # if depletion, blue (negative)
data$y_dummy <- paste0(data$Tissue, "_", data$Traits)

data$Traits <- gsub("Age:Ancestry", "Ancestry:Age", data$Traits)
data$Traits <- factor(data$Traits, levels = gsub("Age:Ancestry", "Ancestry:Age", pw.traits), order = T)
data <- data[order(data$Traits),]

data$y_dummy <- paste0(data$Tissue, "_", data$Traits)
data$y_dummy <- factor(data$y_dummy, levels = rev(paste0(data$Tissue, "_", data$Traits)), order = T)
#data <- data[data$Overlap>=10 & data$OR>1,]
ggplot(data) + 
  geom_point(aes(x=OR, y = factor(y_dummy),
                 col = FDR_bounded,
                 size = as.numeric(Overlap))) +
  geom_errorbar(aes(y = factor(y_dummy),
                    xmin=lower_CI,
                    xmax=upper_CI,
                    col = FDR_bounded)) +
  theme_minimal() +
  #scale_color_manual(values = brewer.pal(11,"RdBu")) +
  scale_color_gradient2(low=brewer.pal(9,"Reds")[[2]], 
                        mid=brewer.pal(9,"Reds")[[4]],
                        high=brewer.pal(9,"Reds")[[8]], midpoint = 3.5) +
  geom_vline(xintercept = 1, lty = 2) +
  labs(colour="-log10(FDR)",
       size = "Number of DMPs") + 
  ylab("") + xlab("Odds ratio") +
  facet_grid(Traits~.,
             scales = "free_y",
             space = "free",
             drop = T,
             switch = "y"
  )  +
  theme(legend.position="bottom")



data$y_dummy <- factor(data$y_dummy, levels = rev(paste0(data$Tissue, "_", data$Traits)), order = T)
ggplot(data) + 
  geom_point(aes(x=OR, y = factor(y_dummy),
                 col = Tissue,
                 size = -log10(fdr))) +
  # geom_errorbar(aes(y = factor(y_dummy),
  #                   xmin=lower_CI,
  #                   xmax=upper_CI,
  #                   col = Tissue)) +
  theme_minimal() +
  scale_color_manual(values = tissue_cols) +
  geom_vline(xintercept = 1, lty = 2) +
  scale_x_continuous(limits=c(0,6))



# Limit color scale (FDR) --
data$FDR_bounded <- ifelse(is.na(data$FDR), NA,
                           ifelse(data$FDR < -3, -3,
                                  ifelse(data$FDR > 3, 3, data$FDR))) # if depletion, blue (negative)
data$Traits <- factor(data$Traits, levels = pw.traits, order = T)



saveRDS(data, "~/marenostrum/Projects/GTEx_v8/Methylation/Data/03.OverlapBetweenVariables.Fisher_results.rds")
