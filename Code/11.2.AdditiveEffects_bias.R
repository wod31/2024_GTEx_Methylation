#!/usr/bin/env Rscript
# @Author: Winona Oliveros
# @E-mail: winona.oliveros@bsc.es
# @Description: Code to get additive effects between demographic traits, test if there is a bias
# @software version: R=4.2.2

#!/usr/bin/env Rscript

# Libraries ----
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(scales)

# ---- Data ----- ####
# Demographic traits ----
traits_cols <- c("Ancestry" = "#E69F00",
                 "Age" = "#56B4E9",
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

head(metadata$WholeBlood)
mdata <- do.call(rbind.data.frame, lapply(tissues, function(tissue) metadata[[tissue]][,c("SUBJID", "EURv1", "AGE", "SEX", "BMI")]))
mdata <- mdata[!duplicated(mdata),]
boxplot((mdata$EURv1~ mdata$SEX))
wilcox.test((mdata$EURv1 ~ mdata$SEX))
boxplot(mdata$AGE ~ mdata$EURv1)
cor(mdata$AGE,mdata$EURv1)
boxplot(mdata$BMI ~ mdata$EURv1)
cor(mdata$BMI,mdata$EURv1)


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

# Pairwise combination of traits --
pw_traits <- c("Ancestry:Sex", 
               "Ancestry:Age",
               "Ancestry:BMI",
               "Sex:Age", 
               "Sex:BMI",
               "Age:BMI"
)

# Are covariates correlated ----
correlated_covariates <- list()

correlated_covariates[["Ancestry:Sex"]] <- sapply(tissues, function(tissue) 
{if(tissue %in% sex_tissues){
  NA
}else{
  wilcox.test(((metadata[[tissue]][,"EURv1"] ~
                              metadata[[tissue]][,"SEX"])), exact = FALSE)$p.value   
}
}
)
correlated_covariates[["Ancestry:Age"]] <- sapply(tissues, function(tissue) 
{
  cor.test(metadata[[tissue]][,"AGE"],
                metadata[[tissue]][,"EURv1"])$p.value   
}
)
correlated_covariates[["Ancestry:BMI"]] <- sapply(tissues, function(tissue) 
{
  cor.test(metadata[[tissue]][,"BMI"],
                metadata[[tissue]][,"EURv1"])$p.value   
}
)
correlated_covariates[["Sex:Age"]] <- sapply(tissues, function(tissue) 
{if(tissue %in% sex_tissues){
  NA
}else{
  wilcox.test(metadata[[tissue]][,"AGE"] ~
                metadata[[tissue]][,"SEX"], exact = FALSE)$p.value   
}
}
)
correlated_covariates[["Sex:BMI"]] <- sapply(tissues, function(tissue) 
{if(tissue %in% sex_tissues){
  NA
}else{
  wilcox.test(metadata[[tissue]][,"BMI"] ~
                metadata[[tissue]][,"SEX"], exact = FALSE)$p.value   
}
}
)
correlated_covariates[["Age:BMI"]] <- sapply(tissues, function(tissue) 
{
  cor.test(metadata[[tissue]][,"AGE"],
           metadata[[tissue]][,"BMI"])$p.value   
}
)

correlation_matrix <- do.call(cbind.data.frame, correlated_covariates)
correlation_matrix <- -log10(correlation_matrix)
correlation_matrix[correlation_matrix < -log10(0.05)] <- NA
pdf('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/Correlation_dem_traits_suppl_additive.pdf', width = 4, height = 3)
Heatmap(as.matrix(correlation_matrix),
        col = brewer.pal(9, "Reds"),
        na_col ="white",
        cluster_rows = F,
        cluster_columns = F)
dev.off()

# Are there tissue:pairwise combinations with more DMPs with 2 traits than expected in a particular direction ----
# Tested only if at least 20 genes overlap

# Is there a bias towards a particular direction of change --
Xsq_fun <- function(tissue, trait1, trait2){
  #print(paste0(tissue, ": ", trait1, "-", trait2))
  if(tissue %in% sex_tissues & (trait1 == "Sex" | trait2 == "Sex")){
    return(list("P-value" = NA, 
                "O/E" = rep(NA, 4),
                "counts" = rep(NA, 4)))
  }else{
    # Do we observe a higher than expected overlap of DMPs in a particular direction of change?
    # a numeric vector representing the observed proportions
    # a vector of probabilities (of the same length of the observed proportions) representing the expected proportions
    trait1.up <- rownames(DMPs_Res[[trait1]][[tissue]][DMPs_Res[[trait1]][[tissue]]$adj.P.Val < 0.05 &
                                                         DMPs_Res[[trait1]][[tissue]]$logFC > 0,])
    trait1.down <- rownames(DMPs_Res[[trait1]][[tissue]][DMPs_Res[[trait1]][[tissue]]$adj.P.Val < 0.05 &
                                                           DMPs_Res[[trait1]][[tissue]]$logFC < 0,])
    
    trait2.up <- rownames(DMPs_Res[[trait2]][[tissue]][DMPs_Res[[trait2]][[tissue]]$adj.P.Val < 0.05 &
                                                        DMPs_Res[[trait2]][[tissue]]$logFC > 0,])
    trait2.down <- rownames(DMPs_Res[[trait2]][[tissue]][DMPs_Res[[trait2]][[tissue]]$adj.P.Val < 0.05 &
                                                           DMPs_Res[[trait2]][[tissue]]$logFC < 0,])
    
    # Observed counts
    counts <- c(sum(trait1.up %in% trait2.up), # upup
                sum(trait1.down %in% trait2.up), # downup
                sum(trait1.up %in% trait2.down), # updown
                sum(trait1.down %in% trait2.down) # downdown
    )
    # Expected proportions
    trait1.up.p <- length(trait1.up)/length(c(trait1.up, trait1.down))
    trait1.down.p <- length(trait1.down)/length(c(trait1.up, trait1.down))
    trait2.up.p <- length(trait2.up)/length(c(trait2.up, trait2.down))
    trait2.down.p <- length(trait2.down)/length(c(trait2.up, trait2.down))
    expected_prob <- c(trait1.up.p * trait2.up.p, trait1.down.p * trait2.up.p, trait1.up.p * trait2.down.p, trait1.down.p * trait2.down.p)
    expected_counts <- round(sum(counts) * expected_prob)
    
    # Return results
    if(sum(counts) < 20){
      print(paste0(tissue, ": Fewer than 20 genes DM with ", trait1, " and ", trait2))
      return(list("P-value" = NA,
                  "O/E" = rep(NA, 4),
                  "counts" = counts))
      break
    }else{
      if(min(expected_counts) < 5){
        print(paste0(tissue, ": Number of observations is not enough for Chi-Square Test\nUsing Monte Carlo simulations"))
        Xsq <- chisq.test(counts, 
                          p = expected_prob,
                          simulate.p.value = T)
      }else{
        Xsq <- chisq.test(counts,
                          p = expected_prob)
      }
      oe <- Xsq$observed/round(Xsq$expected)
      #Xsq$residuals  # Pearson residuals
      #Xsq$stdres     # standardized residuals
      return(list("P-value" = Xsq$p.value,
                  "O/E" = oe,
                  "counts" = counts))
    }
  }
}

# Enrichment analysis ----
Xsq_results <- lapply(pw_traits, function(pw)
  lapply(tissues, function(tissue) 
    Xsq_fun(tissue,
            unlist(strsplit(pw, split = ":"))[[1]],
            unlist(strsplit(pw, split = ":"))[[2]])
  ))
names(Xsq_results) <- pw_traits
for(pw in pw_traits){names(Xsq_results[[pw]]) <- tissues}

# Enrichment statistics --
# P-values
p_values <- lapply(pw_traits, function(pw)
  sapply(tissues, function(tissue)
    Xsq_results[[pw]][[tissue]][["P-value"]]
  )
)
names(p_values) <- pw_traits

# multiple testing correction across tissues and number of pairwise combinations of traits
adj_p_values <- p.adjust(unlist(p_values), method = "BH")
adjPVal_matrix <- matrix(adj_p_values, 
                         nrow = length(tissues), ncol = length(pw_traits),
                         byrow = F)
fdr <- -log10(adjPVal_matrix)
fdr[fdr <= -log10(0.05)] <- 0 # if not tested NA (grey in plot); if FDR >= 0.05, 0 (white in plot)
#range(fdr, na.rm = T)
colnames(fdr) <- pw_traits
rownames(fdr) <- tissue_info$tissue_abbrv
apply(fdr, 2, function(x) sum(x>0, na.rm = T))
range(fdr, na.rm=T)

# tissues and combinations with fdr < 0.05
col_fun <- colorRamp2(seq(-log10(0.05),200,length.out=9),
                      brewer.pal(9, "Reds"))
colnames(fdr) <- c("Ancestry-sex DMPs",
                   "Ancestry-age DMPs",
                   "Ancestry-BMI DMPs",
                   "Sex-age DMPs",
                   "Sex-BMI DMPs",
                   "Age-BMI DMPs")
# row annotation
row_ha_left <-  HeatmapAnnotation("Tissues" = rownames(fdr),
                                  col = list("Tissues" = tissues_cols),
                                  show_legend = F, 
                                  show_annotation_name = T,
                                  simple_anno_size = unit(0.3,"cm"),
                                  annotation_name_gp =  gpar(fontsize = 9),
                                  which = "row")
# column annotation
m_top_bar <- rbind.data.frame(apply(fdr, 2, function(x) sum(x>0, na.rm=T)))
column_bottom_anno <- HeatmapAnnotation("variable1" = sapply(pw_traits, function(i) unlist(strsplit(i, split = ":"))[[1]]),
                                        "variable2" = sapply(pw_traits, function(i) unlist(strsplit(i, split = ":"))[[2]]),
                                        col = list("variable1" = traits_cols,
                                                   "variable2" = traits_cols),
                                        show_legend = F, 
                                        show_annotation_name = F,
                                        simple_anno_size = unit(0.3,"cm"),
                                        annotation_name_gp =  gpar(fontsize = 9)
)
column_top_anno <- HeatmapAnnotation("Number of tissues" =  anno_barplot(t(m_top_bar),
                                                                         gp = gpar(fill = brewer.pal(11, "RdBu")[2],
                                                                                   col = brewer.pal(11, "RdBu")[2]),
                                                                         border = F),
                                     annotation_name_side = "right",
                                     annotation_name_rot = 90,
                                     annotation_name_gp = gpar(fontsize = 8),   
                                     show_legend = T, 
                                     show_annotation_name = T,
                                     height = unit(2.5, "cm"),
                                     gap = unit(0.25,"cm"),
                                     simple_anno_size = unit(0.3,"cm")
)

# main figure ----
ht <- Heatmap(fdr,
              col = col_fun,
              na_col = "white",
              cluster_rows = F,
              cluster_columns = F,
              name = "-log10(FDR)",
              row_names_side = "left",
              column_names_side = "bottom",
              column_names_rot =  90,
              top_annotation = column_top_anno,
              bottom_annotation = column_bottom_anno,
              #left_annotation = row_ha_left,
              row_names_gp = gpar(fontsize = 9),
              column_names_gp = gpar(fontsize = 9),
              heatmap_legend_param = list(direction = "horizontal")
)

draw(ht, 
     heatmap_legend_side = "bottom")



# supplementary tables ----
# observed versus expected ratio
odds_ratio <- lapply(pw_traits, function(pw)
  log2(t(sapply(tissues, function(tissue) 
    Xsq_results[[pw]][[tissue]][["O/E"]]
  )))
)
names(odds_ratio) <- pw_traits  

# number of DMPs 
number_of_DEGs <- lapply(pw_traits, function(pw)
  t(sapply(tissues, function(tissue) 
    Xsq_results[[pw]][[tissue]][["counts"]]
  ))
)
names(number_of_DEGs) <- pw_traits  
for(pw in pw_traits){
  rownames(number_of_DEGs[[pw]]) <- tissue_info$tissue_abbrv
}

# Supplementary tables
column_names_pw_traits <- list(
  "Ancestry:Sex" = c("EUR - female",
                     "AFR - female",
                     "EUR - male",
                     "AFR - male"),
  "Ancestry:Age" = c("EUR - old",
                     "AFR - old",
                     "EUR - young",
                     "AFR - young"),
  "Ancestry:BMI" = c("EUR - high BMI",
                     "AFR - high BMI",
                     "EUR - low BMI",
                     "AFR - low BMI"),
  "Sex:Age" = c("female - old",
                "male - old",
                "female - young",
                "male - young"),
  "Sex:BMI" = c("female - high BMI",
                "male - high BMI",
                "female - low BMI",
                "male - low BMI"),
  "Age:BMI" = c("old - high BMI",
                "young - high BMI",
                "old - low BMI",
                "young - low BMI")
)
names(column_names_pw_traits) <- pw_traits
for(pw in pw_traits){
  colnames(number_of_DEGs[[pw]]) <- column_names_pw_traits[[pw]]
}
for(pw in pw_traits){
  colnames(odds_ratio[[pw]]) <- column_names_pw_traits[[pw]]
}


outpath <- "~/marenostrum/Projects/GTEx_v8/Methylation/Data/"
#dir.create(outpath, recursive = T)
colnames(fdr) <- pw_traits
rownames(correlation_matrix) <- tissue_info$tissue_abbrv
correlation_matrix <- as.data.frame(correlation_matrix)
summary_tables <- list()
for(pw in pw_traits){
  my_df <- cbind.data.frame(number_of_DEGs[[pw]][,c(4, 2, 3,1)],
                            odds_ratio[[pw]][,c(4, 2, 3,1)])
  my_df$p_value <- p_values[[pw]]
  my_df[,"FDR (-log10)"] <- fdr[,pw]
  my_df$correlation_p_value <- sapply(tissue_info$tissue_abbrv, function(tissue) ifelse(is.na(my_df[tissue, "FDR (-log10)"]),
                                                                                        NA, 
                                                                                        correlation_matrix[tissue, pw]))
  write.table(my_df,
              paste0(outpath,
                     gsub(":", "_", pw),
                     ".biased_additive_effects.chi_square_results.tab"),
              col.names = T, row.names = T,
              quote = F,
              sep = "\t")  
  summary_tables[[pw]] <- my_df
}


sapply(pw_traits, function(pw) 
  rownames(summary_tables[[pw]][!is.na(summary_tables[[pw]][,"correlation_p_value"]),])
)
correlation_matrix[is.na(fdr)] <- NA
Heatmap(as.matrix(correlation_matrix),
        col = brewer.pal(9, "Reds"),
        na_col ="white",
        cluster_rows = F,
        cluster_columns = F)

#dir.create(outpath, recursive = T)

#index <- c(4, 2, 3, 1)
ordered_columns <- list(
  "Ancestry:Sex" = c("EUR - female",
                     "AFR - female",
                     "EUR - male",
                     "AFR - male"),
  "Ancestry:Age" = c("EUR - old",
                     "AFR - old",
                     "EUR - young",
                     "AFR - young"),
  "Ancestry:BMI" = c("EUR - high BMI",
                     "AFR - high BMI",
                     "EUR - low BMI",
                     "AFR - low BMI"),
  "Sex:Age" = c("female - old",
                "male - old",
                "female - young",
                "male - young"),
  "Sex:BMI" = c("female - high BMI",
                "male - high BMI",
                "female - low BMI",
                "male - low BMI"),
  "Age:BMI" = c("old - high BMI",
                "young - high BMI",
                "old - low BMI",
                "young - low BMI")
)
pw_colnames <- as.list(colnames(fdr))
names(pw_colnames) <- pw_traits

# ancestry-sex DMPs  ----
pw <- "Ancestry:Sex"    
pw_colname <- pw_colnames[[pw]]
# male - young # down-down
# male - old # down-up
# female - young # up - down
# female - old # up-up

# observed vs expected ratio (cell color)
m <- odds_ratio[[pw]][]
colnames(m) <- ordered_columns[[pw]]
rownames(m) <- tissue_info$tissue_abbrv
range(m, na.rm=T)
col_fun <- colorRamp2(seq(-2,2,length.out=11),
                      rev(brewer.pal(11, "RdBu")))

# number of DMPs (cell value)
m_anno <- number_of_DEGs[[pw]][]
m_anno <- matrix(as.character(m_anno),
                 nrow = length(tissues),
                 ncol = 4,
                 byrow = F)
m_anno[m_anno == "0"] <- ""
colnames(m_anno) <- ordered_columns[[pw]]
rownames(m_anno) <- tissue_info$tissue_abbrv

# fdr (row annotation)
m_fdr <- fdr[rownames(fdr[!is.na(fdr[,pw_colname]) & fdr[,pw_colname]>0,]),pw_colname] # NA: not test; 0: not significant
range(m_fdr, na.rm = T)
col_fun_fdr <-  colorRamp2(c(0, seq(-log10(0.05), max(m_fdr, na.rm=T), length.out=9)),
                           c("white", brewer.pal(9, "Reds")))

m <- m[rownames(fdr[!is.na(fdr[,pw_colname]) & fdr[,pw_colname]>0,]),]
m_anno <- m_anno[rownames(fdr[!is.na(fdr[,pw_colname]) & fdr[,pw_colname]>0,]),]

# row annotation
row_ha_left <-  HeatmapAnnotation("Tissues" = rownames(m),
                                  "FDR" = m_fdr,
                                  col = list("Tissues" = tissues_cols,
                                             "FDR" = col_fun_fdr),
                                  show_legend = F, 
                                  show_annotation_name = T,
                                  simple_anno_size = unit(0.3,"cm"),
                                  annotation_name_gp =  gpar(fontsize = 9),
                                  which = "row")
# column annotation
cols <- c("EUR" = alpha(traits_cols["Ancestry"], 0.5),
          "AFR" = traits_cols["Ancestry"],
          "male" = alpha(traits_cols["Sex"], 0.5),
          "female" = traits_cols["Sex"])
names(cols) <- c("EUR", "AFR", "male", "female")

m_top_bar <- rbind.data.frame(apply(m, 2, function(x) sum(x<0, na.rm=T)),
                              apply(m, 2, function(x) sum(x>0, na.rm=T))
)
column_top_anno <- HeatmapAnnotation("Sex" = sapply(colnames(m), function(i) unlist(strsplit(i, split = " - "))[[2]]),
                                     "Ancestry" = sapply(colnames(m), function(i) unlist(strsplit(i, split = " - "))[[1]]), 
                                     "number of tissues" =  anno_barplot(t(m_top_bar),
                                                                         gp = gpar(fill = brewer.pal(11, "RdBu")[c(10,2)],
                                                                                   col = brewer.pal(11, "RdBu")[c(10,2)])),
                                     col = list("Sex" = cols,
                                                "Ancestry" = cols),
                                     show_legend = T, 
                                     show_annotation_name = T,
                                     simple_anno_size = unit(0.3,"cm"),
                                     annotation_name_gp =  gpar(fontsize = 9)
)
# heatmap
ht <- Heatmap(m,
              col = col_fun,
              na_col = "white",
              cluster_rows = F,
              cluster_columns = F,
              name = "log2(O/E)",
              row_names_side = "left",
              column_names_side = "top",
              column_names_rot =  45,
              top_annotation = column_top_anno,
              left_annotation = row_ha_left,
              row_names_gp = gpar(fontsize = 9),
              column_names_gp = gpar(fontsize = 9),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(m_anno[i, j], x, y, gp = gpar(fontsize = 10))},
              heatmap_legend_param = list(direction = "horizontal")
              
)

draw(ht,
     heatmap_legend_side = "bottom",
     column_title = pw_colnames[[pw]],
     column_title_gp = gpar(fontsize = 10))
# dev.off()


# sex-age DEGs  ----
pw <- "Sex:Age"    
pw_colname <- pw_colnames[[pw]]
# male - young # down-down
# male - old # down-up
# female - young # up - down
# female - old # up-up

# observed vs expected ratio (cell color)
m <- odds_ratio[[pw]][]
colnames(m) <- ordered_columns[[pw]]
rownames(m) <- tissue_info$tissue_abbrv
range(m, na.rm=T)
col_fun <- colorRamp2(seq(-2,2,length.out=11),
                      rev(brewer.pal(11, "RdBu")))

# number of DEGs (cell value)
m_anno <- number_of_DEGs[[pw]][]
m_anno <- matrix(as.character(m_anno),
                 nrow = length(tissues),
                 ncol = 4,
                 byrow = F)
m_anno[m_anno == "0"] <- ""
colnames(m_anno) <- ordered_columns[[pw]]
rownames(m_anno) <- tissue_info$tissue_abbrv

# fdr (row annotation)
m_fdr <- fdr[rownames(fdr[!is.na(fdr[,pw_colname]) & fdr[,pw_colname]>0,]),pw_colname] # NA: not test; 0: not significant
range(m_fdr, na.rm = T)
col_fun_fdr <-  colorRamp2(c(0, seq(-log10(0.05), max(m_fdr, na.rm=T), length.out=9)),
                           c("white", brewer.pal(9, "Reds")))

m <- m[rownames(fdr[!is.na(fdr[,pw_colname]) & fdr[,pw_colname]>0,]),]
m_anno <- m_anno[rownames(fdr[!is.na(fdr[,pw_colname]) & fdr[,pw_colname]>0,]),]

# row annotation
row_ha_left <-  HeatmapAnnotation("Tissues" = rownames(m),
                                  "FDR" = m_fdr,
                                  col = list("Tissues" = tissues_cols,
                                             "FDR" = col_fun_fdr),
                                  show_legend = F, 
                                  show_annotation_name = T,
                                  simple_anno_size = unit(0.3,"cm"),
                                  annotation_name_gp =  gpar(fontsize = 9),
                                  which = "row")
# column annotation
cols <- c("male" = alpha(traits_cols["Sex"], 0.5),
          "female" = traits_cols["Sex"],
          "young" = alpha(traits_cols["Age"], 0.5),
          "old" = traits_cols["Age"])
names(cols) <- c("male", "female", "young", "old")

m_top_bar <- rbind.data.frame(apply(m, 2, function(x) sum(x<0, na.rm=T)),
                              apply(m, 2, function(x) sum(x>0, na.rm=T))
)
column_top_anno <- HeatmapAnnotation("Age" = sapply(colnames(m), function(i) unlist(strsplit(i, split = " - "))[[2]]),
                                     "Sex" = sapply(colnames(m), function(i) unlist(strsplit(i, split = " - "))[[1]]), 
                                     "number of tissues" =  anno_barplot(t(m_top_bar),
                                                                         gp = gpar(fill = brewer.pal(11, "RdBu")[c(10,2)],
                                                                                   col = brewer.pal(11, "RdBu")[c(10,2)])),
                                     col = list("Sex" = cols,
                                                "Age" = cols),
                                     show_legend = T, 
                                     show_annotation_name = T,
                                     simple_anno_size = unit(0.3,"cm"),
                                     annotation_name_gp =  gpar(fontsize = 9)
)
# heatmap
ht <- Heatmap(m,
              col = col_fun,
              na_col = "white",
              cluster_rows = F,
              cluster_columns = F,
              name = "log2(O/E)",
              row_names_side = "left",
              column_names_side = "top",
              column_names_rot =  45,
              top_annotation = column_top_anno,
              left_annotation = row_ha_left,
              row_names_gp = gpar(fontsize = 9),
              column_names_gp = gpar(fontsize = 9),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(m_anno[i, j], x, y, gp = gpar(fontsize = 10))},
              heatmap_legend_param = list(direction = "horizontal")
              
)

# pdf(paste0(outpath, pw, "_DEGs.chi_square_heatmap.pdf"),
#     width = 3, height = 5)
draw(ht,
     heatmap_legend_side = "bottom",
     column_title = pw_colnames[[pw]],
     column_title_gp = gpar(fontsize = 10))
# dev.off()
