#########################################
#### correlation DVG and DVP #####
#########################################
#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez & Winona Oliveros
# @E-mail: jose.ramirez1@bsc.es & winona.oliveros@bsc.es
# @Description: Code to correlate gene expression and DNA methylation at the variability level
# @software version: R=4.2.2

#Loading libraries
library(dplyr)

#Set path 
setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")

#Correlate all DVPs to all DVGs:

# Parsing
library(optparse)
parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-t", "--tissue"), type="character",
                     dest="tissue",
                     help="Tissue")
options=parse_args(parser)
tissue=options$tissue
# tissue <- "Lung"

print(tissue)


print("Reading annotation")
#Get a table with probe and a gene
Sys.time()
data_path <- "/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/"
annotation <- read.delim(paste0(data_path, "Methylation_Epic_gene_promoter_enhancer_processed.txt"), sep = '\t', header = T)
Sys.time()


#From ensembl id to gene symbol
gene_annotation <- read.delim("/gpfs/projects/bsc83/MN4/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/data/public/gencode.v26.GRCh38.genes.bed", header=F)[,c(6,7)]
colnames(gene_annotation) <- c("gene", "symbol")

#These were genes duplicated, I changed their names to their correct one
gene_annotation$symbol[gene_annotation$gene=="ENSG00000253972.5"] <- "MAL2-AS1" #I found this a posteriori
gene_annotation$symbol[gene_annotation$gene=="ENSG00000283992.1"] <- "SLURP2" #Insted of LYNX1
gene_annotation$symbol[gene_annotation$gene=="ENSG00000235271.5"] <- "GC22P027299"
gene_annotation$symbol[gene_annotation$gene=="ENSG00000229694.6"] <- "C9orf73"
gene_annotation$symbol[gene_annotation$gene=="ENSG00000228741.2"] <- "GC13P024553" 


print("Running")
print(tissue)
Sys.time()

### read DEG ####

sexual_tissues <- c("Prostate", "Testis", "Ovary")
if(tissue %in% sexual_tissues){
  print("Sexual tissue")
  traits_meth <- c("EURv1", "AGE")
  traits_expr <- c("Ancestry", "Age")
  names(traits_meth) <- traits_expr
  names(traits_expr) <- traits_meth
  
} else{
  traits_meth <- c("EURv1", "SEX2", "AGE")
  traits_expr <- c("Ancestry", "Sex", "Age")
  names(traits_meth) <- traits_expr
  names(traits_expr) <- traits_meth
  
}

for (trait in traits_expr) {
  #res <- results_DML[[trait]]
  res <- readRDS(paste0("Tissues/",tissue, "/DVP_",trait,".rds"))
  
  signif <- res[res$P.Value<0.05,]
  # table(signif$logFC>0)
  #Reading methylationresiduals
  beta <- readRDS(paste0("Tissues/", tissue, '/',traits_meth[trait],"_methylation_residuals.continous.rds"))
  Sys.time() 
  
  #Reading expression residuals in the lung:
  expression <- readRDS(paste0("Tissues/", tissue, '/',traits_meth[trait],"_expression_residuals.continous.rds"))
  
  rownames(expression) <- sapply(rownames(expression), function(gene) gene_annotation$symbol[gene_annotation$gene==gene])
  #From sample id to donor id
  colnames(expression) <- sapply(colnames(expression), function(id) paste0(strsplit(id, "-")[[1]][-3], collapse="-"))
  #Subset expression data to match the donors in DNA methylation data
  expression <- expression[,colnames(expression) %in% colnames(beta)]
  beta <- beta[,colnames(beta) %in% colnames(expression)] 
  identical(colnames(beta), colnames(expression)) #same order
  
  Sys.time() 
  print("Starting to compute correlations")
  
  #Compute correlations
  
  #Now correlate beta values for each position and expression for its associated gene and compute Pearson's correlation
  correlation_function <- function(gene, probe){
    if(!gene %in% rownames(expression)){ return(NA) }
    exp <- expression[rownames(expression)==gene,]
    if(!probe %in% rownames(beta)){ return(NA) }
    bet <- beta[rownames(beta)==probe,]
    test <- cor.test(exp, as.numeric(bet))
    to_return <- list("p.val"= test$p.value, "cor"=test$estimate, "gene"=gene, "probe"=probe)
    return(to_return)
  }
  
  #Running a different analysis per promoter, enhancer and gene body
  annotation_promoter <- annotation[annotation$Type=="Promoter_Associated",]
  annotation_enhancer <- annotation[annotation$Type=="Enhancer_Associated",]
  annotation_gene_body <- annotation[annotation$Type=="Gene_Associated",]
  

  expr <- readRDS(paste0("Tissues/",tissue, "/DVG_",trait,".rds"))
  deg <-  rownames(expr[expr$Adj.P.Value<0.05,])
  
  deg_symbol <- gene_annotation$symbol[gene_annotation$gene %in% deg]
  
  annotation_promoter <- annotation_promoter[annotation_promoter$IlmnID %in% rownames(signif) & annotation_promoter$UCSC_RefGene_Name %in% deg_symbol,]
  annotation_enhancer <- annotation_enhancer[annotation_enhancer$IlmnID %in% rownames(signif) & annotation_enhancer$UCSC_RefGene_Name %in% deg_symbol,]
  annotation_gene_body <- annotation_gene_body[annotation_gene_body$IlmnID %in% rownames(signif) & annotation_gene_body$UCSC_RefGene_Name %in% deg_symbol,]
  #annotation_other <- annotation_other[annotation_other$IlmnID %in% rownames(signif),]
  
  if (nrow(annotation_promoter)>0) {
    annotation_promoter_df <- annotation_promoter$UCSC_RefGene_Name
    names(annotation_promoter_df) <- annotation_promoter$IlmnID
    annotation_promoter_df <- stack(annotation_promoter_df) #From list to data frame
    annotation_promoter_df <- annotation_promoter_df %>% distinct() #Keep only non duplicated rows
    annotation_promoter_df$ind <- as.character(annotation_promoter_df$ind)
  } else {
    annotation_promoter_df <- NA
  }
  if (nrow(annotation_enhancer)>0) {
    annotation_enhancer_df <- annotation_enhancer$UCSC_RefGene_Name
    names(annotation_enhancer_df) <- annotation_enhancer$IlmnID
    annotation_enhancer_df <- stack(annotation_enhancer_df) #From list to data frame
    annotation_enhancer_df <- annotation_enhancer_df %>% distinct() #Keep only non duplicated rows
    annotation_enhancer_df$ind <- as.character(annotation_enhancer_df$ind)
    
  } else {
    annotation_enhancer_df <- NA
  }
  if (nrow(annotation_gene_body)>0) {
    annotation_gene_body_df <- annotation_gene_body$UCSC_RefGene_Name
    names(annotation_gene_body_df) <- annotation_gene_body$IlmnID
    annotation_gene_body_df <- stack(annotation_gene_body_df) #From list to data frame
    annotation_gene_body_df <- annotation_gene_body_df %>% distinct() #Keep only non duplicated rows
    annotation_gene_body_df$ind <- as.character(annotation_gene_body_df$ind)
  } else {
    annotation_gene_body_df <- NA
  }
  
  if (is.na(annotation_promoter_df)[1]) {
    output_promoters <- list('prom'=list("p.val"= 1, "cor"=0, "gene"=NA, "probe"=NA))
    #names(output_promoters) <- c('prom')
  } else{
    output_promoters <- mapply(function(x,y) correlation_function(x, y), annotation_promoter_df[,1], annotation_promoter_df[,2] , SIMPLIFY = FALSE)
  }
  if (is.na(annotation_enhancer_df)[1]) {
    output_enhancers <- list('enh'=list("p.val"= 1, "cor"=0, "gene"=NA, "probe"=NA))
    #names(output_enhancers) <- c(NA)
  } else{
    output_enhancers <- mapply(function(x,y) correlation_function(x, y), annotation_enhancer_df[,1], annotation_enhancer_df[,2] , SIMPLIFY = FALSE)
  }
  if (is.na(annotation_gene_body_df)[1]) {
    output_gene_body <- list('gene'=list("p.val"= 1, "cor"=0, "gene"=NA, "probe"=NA))
    #names(output_gene_body) <- c(NA)
  } else{
    output_gene_body <- mapply(function(x,y) correlation_function(x, y), annotation_gene_body_df[,1], annotation_gene_body_df[,2] , SIMPLIFY = FALSE)
  }
  #output_other <- mapply(function(x,y) correlation_function(x, y), annotation_other_df[,1], annotation_other_df[,2] )
  
  cleaning_output <- function(output){
    
    output_2 <- as.data.frame(t(do.call(cbind,output)))
    
    output_2$cor <- unlist(output_2$cor)
    output_2$p.val <- unlist(output_2$p.val)
    output_2$gene <- unlist(output_2$gene)
    output_2$probe <- unlist(output_2$probe)
    
    output_2 <- output_2[!is.na(output_2$p.val),]   #Removing NAs
    
    output_2$p.adj <- p.adjust(output_2$p.val, method = "BH")
    return(output_2)
  }
  output_promoters <- cleaning_output(output_promoters)
  output_promoters$class <- "promoter"
  output_enhancers <- cleaning_output(output_enhancers)
  output_enhancers$class <- "enhancer"
  output_gene_body <- cleaning_output(output_gene_body)
  output_gene_body$class <- "gene_body"
  # output_other <- cleaning_output(output_other)
  # output_other$class <- "other"
  
  print(paste0("Finished computing correlations for", tissue))
  Sys.time()
  
  output <- rbind(output_promoters, output_enhancers, output_gene_body)
  # saveRDS(output, paste0("tissues/Lung/Correlations.rds"))
  
  saveRDS(output, paste0("Tissues/", tissue, '/',trait,"_Correlations_probes_genes_DVG_DVP.pnominal.rds"))
  
}
#---------------------------------------------------










