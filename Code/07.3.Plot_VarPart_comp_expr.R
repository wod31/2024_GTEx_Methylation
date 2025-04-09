##### plot hier.part results ##### overlap with expression
#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Plot variance partition results on DNA methylation data per tissue on demographic traits + overlap gene expression
# @software version: R=4.2.2


rm(list=ls())
suppressPackageStartupMessages(library(ComplexHeatmap))
library(RColorBrewer)
suppressPackageStartupMessages(library(circlize))

first_dir <- "~/marenostrum/"

project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "WholeBlood", "KidneyCortex", "Testis","MuscleSkeletal")
names <- c("Age", "Ancestry", "BMI", "Sex")
data <- matrix(nrow = length(tissues), ncol=4, dimnames = list(tissues, names))

tissue_info <- readRDS(paste0(project_path, "Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

sex_tissues <- c('Ovary','Prostate','Testis')

### hier.part
# proportion of total tissue methylation variation explained by each trait --
get_tissue_methylation_variation_explained <- function(tissue){
  print(tissue)
  
  hier_res <- readRDS(paste0(project_path, "varPart/",tissue, "_var_part.rds"))
  
  if(tissue %in% sex_tissues){
    tissue_methylation_variation_explained <- sapply(traits[-2], function(trait) sum(hier_res[,(trait)]))/sum(sapply(traits[-2], function(trait) sum(hier_res[,(trait)])))
    tissue_methylation_variation_explained <- c(tissue_methylation_variation_explained[c(1)], 0, tissue_methylation_variation_explained[c(2,3)])
    names(tissue_methylation_variation_explained)[2] <- "Sex"
  }else{
    tissue_methylation_variation_explained <- sapply(traits, function(trait) sum(hier_res[,(trait)]))/sum(sapply(traits, function(trait) sum(hier_res[,(trait)])))
  }
  return(tissue_methylation_variation_explained)  
}


traits <- c('EURv1','SEX','AGE','BMI')
tissue_methylation_variation_explained <- do.call(rbind.data.frame,
                                                 lapply(tissues, function(tissue) get_tissue_methylation_variation_explained(tissue)))
colnames(tissue_methylation_variation_explained) <- c("Ancestry","Sex","Age","BMI")
rownames(tissue_methylation_variation_explained) <- tissues    

traits_cols <- c('#C49122','#4B8C61','#70A0DF','#A76595')
names(traits_cols) <- c("Ancestry","Sex","Age","BMI")

# bar plot
plot_path <- paste0(project_path, "Plots/")
pdf(paste0(plot_path, "Figure_S2B.tissue_methylation_variation_explained_new.meth.pdf"),
    width = 3, height = 4)
# width = 3, height = 9)
barplot(t(tissue_methylation_variation_explained[length(tissues):1,]),
        horiz = T,
        border = NA,
        col = traits_cols,
        xlab = "Tissue methylation \n variation explained (%)",
        las = 2, 
        cex.names = 0.9,
        xaxt = 'n')#,
        #yaxt = 'n')
axis(1, at = axTicks(1))
dev.off()

#### get genes from dna methylation 
gene_annotation <- read.delim("~/marenostrum/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/data/public/gencode.v26.GRCh38.genes.bed", header=F)[,c(6,7)]

colnames(gene_annotation) <- c("gene", "symbol")

#These were genes duplicated, I changed their names to their correct one
gene_annotation$symbol[gene_annotation$gene=="ENSG00000253972.5"] <- "MAL2-AS1" #I found this a posteriori
gene_annotation$symbol[gene_annotation$gene=="ENSG00000283992.1"] <- "SLURP2" #Insted of LYNX1
gene_annotation$symbol[gene_annotation$gene=="ENSG00000235271.5"] <- "GC22P027299"
gene_annotation$symbol[gene_annotation$gene=="ENSG00000229694.6"] <- "C9orf73"
gene_annotation$symbol[gene_annotation$gene=="ENSG00000228741.2"] <- "GC13P024553" 

data_path <- "~/marenostrum/Projects/GTEx_v8/Methylation/Data/"
annotation <- read.delim(paste0(data_path, "Methylation_Epic_gene_promoter_enhancer_processed.txt"), sep = '\t', header = T)
Sys.time()

get_genes <- function(tissue) {
  hier_res <- readRDS(paste0(project_path, "varPart/",tissue, "_var_part.rds"))
  cpgs <- rownames(hier_res)
  genes <- annotation$UCSC_RefGene_Name[annotation$IlmnID %in% cpgs]
  
  ensembl <- gene_annotation$gene[gene_annotation$symbol %in% genes]
  return(ensembl)
}

genes_meth <- lapply(tissues[tissues != 'KidneyCortex'], function(tissue) get_genes(tissue))
names(genes_meth) <- tissues[tissues != 'KidneyCortex']


### expression ####
get_tissue_expression_variation_explained <- function(tissue){
  print(tissue)

  hier_res <- readRDS(paste0(project_path, "varPart/",tissue, "_var_part_expression.rds"))
  hier_res <- hier_res[genes_meth[[tissue]],]
  
  if(tissue %in% sex_tissues){
    tissue_expression_variation_explained <- sapply(traits[-2], function(trait) sum(hier_res[,(trait)]))/sum(sapply(traits[-2], function(trait) sum(hier_res[,(trait)])))
    tissue_expression_variation_explained <- c(tissue_expression_variation_explained[c(1)], 0, tissue_expression_variation_explained[c(2,3)])
    names(tissue_expression_variation_explained)[2] <- "Sex"
  }else{
    tissue_expression_variation_explained <- sapply(traits, function(trait) sum(hier_res[,(trait)]))/sum(sapply(traits, function(trait) sum(hier_res[,(trait)])))
  }
  return(tissue_expression_variation_explained)  
}


traits <- c('EURv1','Sex','Age','BMI')
tissue_expression_variation_explained <- do.call(rbind.data.frame,
                                                 lapply(tissues[tissues != 'KidneyCortex'], function(tissue) get_tissue_expression_variation_explained(tissue)))
colnames(tissue_expression_variation_explained) <- c("Ancestry","Sex","Age","BMI")
rownames(tissue_expression_variation_explained) <- tissues[tissues != 'KidneyCortex']    

traits_cols <- c('#C49122','#4B8C61','#70A0DF','#A76595')
names(traits_cols) <- c("Ancestry","Sex","Age","BMI")

# bar plot
plot_path <- paste0(project_path, "Plots/")
pdf(paste0(plot_path, "Figure_S2B.tissue_expression_variation_explained_new.expr.pdf"),
    width = 3, height = 4)
# width = 3, height = 9)
barplot(t(tissue_expression_variation_explained[length(tissues):1,]),
        horiz = T,
        border = NA,
        col = traits_cols,
        xlab = "Tissue gene expression \n variation explained (%)",
        las = 2, 
        cex.names = 0.9,
        xaxt = 'n')#,
#yaxt = 'n')
axis(1, at = axTicks(1))
#axis(2)#, labels = rownames(tissue_expression_variation_explained))
dev.off()


