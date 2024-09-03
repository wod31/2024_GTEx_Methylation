##### plot hier.part results ##### overlap with expression
rm(list=ls())
suppressPackageStartupMessages(library(ComplexHeatmap))
library(RColorBrewer)
suppressPackageStartupMessages(library(circlize))

first_dir <- "~/marenostrum/"

project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "WholeBlood", "KidneyCortex", "Testis","MuscleSkeletal")
names <- c("Age", "Ancestry", "BMI", "Sex")
data <- matrix(nrow = length(tissues), ncol=4, dimnames = list(tissues, names))

#tissue_info <- readRDS(paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Jose/00_Data/Tissue_info_whole.rds"))
tissue_info <- readRDS(paste0(project_path, "Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

sex_tissues <- c('Ovary','Prostate','Testis')

### hier.part
# proportion of total tissue expression variation explained by each trait --
get_tissue_expression_variation_explained <- function(tissue){
  print(tissue)
  
  # if (tissue %in% c('Testis')) {
  #   hier_res <- readRDS(paste0(project_path, "Tissues/",tissue, "/hier.part.5peers.continous.rds"))
  #   tissue_expression_variation_explained <- sapply(traits[c(3:4)], function(trait) sum(hier_res[,paste0(trait,'_abs')]))/sum(sapply(traits[c(3:4)], function(trait) sum(hier_res[,paste0(trait,'_abs')])))
  #   tissue_expression_variation_explained <- c(0, 0, tissue_expression_variation_explained[c(1,2)])
  #   names(tissue_expression_variation_explained)[1] <- "Ancestry"
  #   names(tissue_expression_variation_explained)[2] <- "Sex"
  #   return(tissue_expression_variation_explained)  
  # } else {
  hier_res <- readRDS(paste0(project_path, "varPart/",tissue, "_var_part.rds"))
  
  if(tissue %in% sex_tissues){
    tissue_expression_variation_explained <- sapply(traits[-2], function(trait) sum(hier_res[,(trait)]))/sum(sapply(traits[-2], function(trait) sum(hier_res[,(trait)])))
    tissue_expression_variation_explained <- c(tissue_expression_variation_explained[c(1)], 0, tissue_expression_variation_explained[c(2,3)])
    names(tissue_expression_variation_explained)[2] <- "Sex"
  }else{
    tissue_expression_variation_explained <- sapply(traits, function(trait) sum(hier_res[,(trait)]))/sum(sapply(traits, function(trait) sum(hier_res[,(trait)])))
  }
  return(tissue_expression_variation_explained)  
}


traits <- c('EURv1','SEX','AGE','BMI')
tissue_expression_variation_explained <- do.call(rbind.data.frame,
                                                 lapply(tissues, function(tissue) get_tissue_expression_variation_explained(tissue)))
colnames(tissue_expression_variation_explained) <- c("Ancestry","Sex","Age","BMI")
rownames(tissue_expression_variation_explained) <- tissues    

traits_cols <- c('#C49122','#4B8C61','#70A0DF','#A76595')
names(traits_cols) <- c("Ancestry","Sex","Age","BMI")

# bar plot
plot_path <- paste0(project_path, "Plots/")
pdf(paste0(plot_path, "Figure_S2B.tissue_expression_variation_explained_new.meth.pdf"),
    width = 3, height = 4)
# width = 3, height = 9)
barplot(t(tissue_expression_variation_explained[length(tissues):1,]),
        horiz = T,
        border = NA,
        col = traits_cols,
        xlab = "Tissue methylation \n variation explained (%)",
        las = 2, 
        cex.names = 0.9,
        xaxt = 'n')#,
        #yaxt = 'n')
axis(1, at = axTicks(1))
#axis(2)#, labels = rownames(tissue_expression_variation_explained))
dev.off()

#### get genes from dna methylation 
gene_annotation <- read.delim("~/marenostrum/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/data/public/gencode.v26.GRCh38.genes.bed", header=F)[,c(6,7)]
#gene_annotation <- read.delim("~/marenostrum/MN4/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/data/public/gencode.v26.GRCh38.genes.bed", header=F)[,c(6,7)]

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
  
  # if (tissue %in% c('Testis')) {
  #   hier_res <- readRDS(paste0(project_path, "Tissues/",tissue, "/hier.part.5peers.continous.rds"))
  #   tissue_expression_variation_explained <- sapply(traits[c(3:4)], function(trait) sum(hier_res[,paste0(trait,'_abs')]))/sum(sapply(traits[c(3:4)], function(trait) sum(hier_res[,paste0(trait,'_abs')])))
  #   tissue_expression_variation_explained <- c(0, 0, tissue_expression_variation_explained[c(1,2)])
  #   names(tissue_expression_variation_explained)[1] <- "Ancestry"
  #   names(tissue_expression_variation_explained)[2] <- "Sex"
  #   return(tissue_expression_variation_explained)  
  # } else {
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
        xlab = "Tissue methylation \n variation explained (%)",
        las = 2, 
        cex.names = 0.9,
        xaxt = 'n')#,
#yaxt = 'n')
axis(1, at = axTicks(1))
#axis(2)#, labels = rownames(tissue_expression_variation_explained))
dev.off()


#### prop explained per tissue dotplot 
DEGs_per_trait_and_tissue <- unlist(lapply(traits, function(trait) sapply(tissues, function(tissue)  ifelse(tissue %in% sex_tissues & trait == "Sex", NA, length(DEGs[[trait]][[tissue]])))))
DEGs_per_tissue <- sapply(tissues, function(tissue) length(unique(unlist(lapply(traits, function(trait) DEGs[[trait]][[tissue]])))[!is.na(unique(unlist(lapply(traits, function(trait) DEGs[[trait]][[tissue]]))))]))
proporton_of_DEGs_per_trait <- 100*DEGs_per_trait_and_tissue/rep(DEGs_per_tissue,4)

get_tissue_mean_variation <- function(tissue){
  print(tissue)
  
  hier_res <- readRDS(paste0(project_path, "varPart/",tissue, "_var_part.rds"))
  
  if(tissue %in% sex_tissues){
    tissue_expression_variation_explained <- sapply(traits_prev[-2], function(trait) mean(hier_res[,(trait)]))
    tissue_expression_variation_explained <- c(tissue_expression_variation_explained[c(1)], 0, tissue_expression_variation_explained[c(2,3)])
    names(tissue_expression_variation_explained)[2] <- "Sex"
  }else{
    tissue_expression_variation_explained <- sapply(traits_prev, function(trait) mean(hier_res[,(trait)]))
  }
  return(tissue_expression_variation_explained)  
}

traits_prev <- c('EURv1','SEX','AGE','BMI')
tissue_expression_variation_explained <- do.call(rbind.data.frame,
                                                 lapply(tissues, function(tissue) get_tissue_mean_variation(tissue)))
colnames(tissue_expression_variation_explained) <- traits
rownames(tissue_expression_variation_explained) <- tissues  

# average gene expression variation explained per trait and tissue --
data <- cbind.data.frame("Tissue" = rep(tissues, 4),
                         "Trait" = unlist(lapply(traits, function(trait) rep(trait, length(tissues)))),
                         "DEGs" = unlist(lapply(traits, function(trait) sapply(tissues, function(tissue)  ifelse(tissue %in% sex_tissues & trait == "Sex", NA, length(DEGs[[trait]][[tissue]]))))),
                         "proportion_of_DEGs" <- proporton_of_DEGs_per_trait,
                         "R2" = unlist(tissue_expression_variation_explained)
)
data$Tissue <- factor(data$Tissue, levels = rev(tissues), order = T)
data$Trait <- factor(data$Trait, levels = traits, order = T)

pdf(paste0(plot_path, "Figure_1C.gene_expression_variation_explained_new.pdf"),
    width = 3, height = 4)
ggplot(data,
       aes(x = Tissue,
           y = R2,
           col = Trait,
           size = proportion_of_DEGs)) +
  geom_point(alpha = 0.5) +
  coord_flip() +
  theme_bw() +
  scale_color_manual(values = traits_cols) +
  ylab("Mean CpG methylation \n variation explained (%)") +
  xlab("") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12))
dev.off()

pdf(paste0(plot_path, "Figure_1C.gene_expression_variation_explained.legend.pdf"),
    width = 5, height = 9)
ggplot(data,
       aes(x = Tissue,
           y = R2,
           col = Trait,
           size = proportion_of_DEGs)) +
  geom_point(alpha = 0.5) +
  coord_flip() +
  theme_bw() +
  scale_color_manual(values = traits_cols) +
  ylab("Mean CpG methylation \n variation explained (%)") +
  xlab("") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12))
dev.off()

#### get CpGs for which a single trait explains more than 10% of their inter-individual DNa methylation variation #####
get_tissue_mean_variation <- function(tissue){
  print(tissue)
  
  hier_res <- readRDS(paste0(project_path, "varPart/",tissue, "_var_part.rds"))
  
  if(tissue %in% sex_tissues){
    tissue_expression_variation_explained <- hier_res[,traits_prev[-2]]
  }else{
    tissue_expression_variation_explained <- hier_res[,traits_prev]
  }
  return(tissue_expression_variation_explained)  
}

traits_prev <- c('EURv1','SEX','AGE','BMI')
tissue_expression_variation_explained_all <- lapply(tissues, function(tissue) get_tissue_mean_variation(tissue))
names(tissue_expression_variation_explained_all) <- tissues

#install.packages('yarrr')
library(yarrr)
library(reshape2)
library(scales)

# 5.1.1 Prepare expression data ----
# Function to parse tissue data
melt.data <- function(tissue, data){
  print(tissue)
  df <- data[[tissue]]
  df$cpg <- rownames(df)
  if(tissue %in% sex_tissues){
    d <- reshape2::melt(df[,c((traits_prev[-2]),"cpg")],
                        id.vars = "cpg",
                        variable.name = 'Trait',
                        value.name = 'R2')
  }else{
    d <- reshape2::melt(df[,c((traits_prev),"cpg")],
                        id.vars = "cpg",
                        variable.name = 'Trait',
                        value.name = 'R2')
    
  }
  min(d$R2)
  #d$Trait <- gsub(pattern = "_abs",replacement = "", d$Trait)
  d$Tissue <- rep(tissue_info[tissue_info$tissue_ID==tissue, "tissue_abbrv"], nrow(d))
  return(d)
}

meth.data <- do.call(rbind.data.frame,
                     lapply(tissues, function(tissue) 
                       melt.data(tissue, tissue_expression_variation_explained_all)
                     ))
meth.data <- meth.data[!is.na(meth.data$R2),]

# total unique number of CpGs per tissue

meth.data <- merge(meth.data, tissue_info[,c("tissue_ID", "tissue_abbrv")], by.x='Tissue', by.y='tissue_abbrv')

dma_10 <- list()
traits <- c("Ancestry", "Sex", "Age", 'BMI')
names(traits) <- traits_prev
for(tissue in tissues){
  dma_10[[tissue]] <- list()
  for(trait in traits_prev){
    if(tissue %in% sex_tissues & trait=="SEX"){
      dma_10[[tissue]][[trait]] <-  NA
    }else{
      trait_deg <- 
        dma_10[[tissue]][[trait]] <-  meth.data[meth.data$cpg %in% DEGs[[traits[trait]]][[tissue]] & meth.data$R2>0.1 & meth.data$tissue_ID == tissue & meth.data$Trait == trait,]  
    }
  }
}

get_DMPs <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "SEX"){
    NA
  }else{
    unique(dma_10[[tissue]][[trait]]$cpg)
  }
}

DMPs <- lapply(traits_prev, function(trait) lapply(tissues, function(tissue) get_DMPs(tissue, trait)))
DMPs_ns <- lapply(traits_prev, function(trait) lapply(tissues[!tissues %in% sex_tissues], function(tissue) get_DMPs(tissue, trait)))

names(DMPs) <- traits_prev
for(trait in traits_prev){names(DMPs[[trait]]) <- tissues}
names(DMPs_ns) <- traits_prev
for(trait in traits_prev){names(DMPs_ns[[trait]]) <- tissues[!tissues %in% sex_tissues]}

genes <- unique(unlist(lapply(traits_prev, function(trait) unlist(DMPs[[trait]]))))
length(genes)

### get CpGs per trait and perform enrichment #####
cpgs_traits <- lapply(traits_prev, function(trait) unlist(DMPs[[trait]]))
names(cpgs_traits) <- traits_prev
cpgs_traits_2 <- lapply(traits_prev, function(trait) unlist(DMPs_ns[[trait]]))
names(cpgs_traits_2) <- traits_prev
cpgs_traits_all <- lapply(traits_prev, function(trait) unique(cpgs_traits[[trait]],cpgs_traits_2[[trait]]))
names(cpgs_traits_all) <- traits_prev

#### enrichment + annotation 
annot <- read.delim('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Methylation_Epic_gene_promoter_enhancer_processed.txt')
library(missMethyl)

betas <- read.csv("~/marenostrum_scratch/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")$Name
res <- list()
for (trait in traits_prev) {
  print(trait)
  res_1 <- gometh(cpgs_traits[[trait]], all.cpg=(betas),
                  collection="GO", array.type="EPIC")
  res[[trait]] <- res_1#[res_1$ONTOLOGY=="BP",]
  print(table(res[[trait]]$FDR<0.05))
}

### plot examples and number of explained
traits <- c("Ancestry","Sex","Age","BMI")
driver_traits <- (sapply(traits_prev, function(x) length(cpgs_traits[[x]])))
names(driver_traits) <- traits
#driver_traits['Ancestry'] <- 3

# top annotation bar plot
plot_path <- '~/marenostrum/Projects/GTEx_v8/Methylation/Plots/'
pdf(paste0(plot_path, "Figure_2C.barplot_n10_variation_explained_new.nosex.pdf"),
    width = 4, height = 4)
barplot(driver_traits, col = traits_cols, names.arg=c("Ancestry","Sex","Age","BMI"), ylab = 'Number of CpGs R2 > 0.1')
dev.off()
plot.new()



### plots like in cell genomics
get_box_stats <- function(y, upper_limit = max(y) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "n =", length(y)
    )
  ))
}

admixture_ancestry <- read.table('~/marenostrum_scratch/MN4/bsc83/bsc83535/GTEx/v8/genotype_data/admixture_inferred_ancestry.txt')
colnames(admixture_ancestry) <- c('SUBJID','AFRv1','EURv1','inferred_ancestry','AFRv2','EURv2')
admixture_ancestry$Ancestry <- 'ADX'
admixture_ancestry$Ancestry[admixture_ancestry$EURv1>=0.5] <- 'EUR'
admixture_ancestry$Ancestry[admixture_ancestry$EURv1<0.5] <- 'AFR'

metadata <- lapply(tissues, function(tissue) readRDS(paste0(project_path, "Tissues/",tissue, "/metadata.rds")))
names(metadata) <- tissues


get_gene_data <- function(tissue, df, gene_name){
  #df <- cbind.data.frame(t(tpm[gene_name,]))
  colnames(df) <- "Beta"
  df$Ancestry <- sapply(rownames(df), function(i) admixture_ancestry[admixture_ancestry$SUBJID==i,"Ancestry"])
  df$Ancestry <- gsub("EUR", "EA", df$Ancestry)
  df$Ancestry <- gsub("AFR", "AA", df$Ancestry)
  df$Ancestry <- factor(df$Ancestry, levels = c("EA", "AA"), order = T)
  df$Sex <- sapply(rownames(df), function(i) metadata[[tissue]][metadata[[tissue]]$SUBJID==i,"SEX"])
  df$Sex <- gsub("1", "Male", df$Sex)
  df$Sex <- gsub("2", "Female", df$Sex)
  df$Sex <- factor(df$Sex, levels = c("Male", "Female"), order = T)
  df$Age_int <- sapply(rownames(df), function(i) metadata[[tissue]][metadata[[tissue]]$SUBJID==i,"AGE"])
  df$Age <- sapply(df$Age_int, function(i) ifelse(i<45, "[20-45)", "[45-70]"))
  df$Age <- factor(df$Age, 
                   levels = c("[20-45)", "[45-70]"),
                   order = T)
  # df$BMI_int <- sapply(rownames(df), function(i) mdata[[tissue]][mdata[[tissue]]$Sample==i,"BMI"])
  # df$BMI <- sapply(df$BMI_int, function(bmi) pseudo_categorize_bmi(bmi))
  # df$BMI <- factor(df$BMI,
  #                  levels = c("Normal", "Overweight", "Obese"),
  #                  order = T)
  df$Tissue <- rep(tissue_info[tissue_info$tissue_ID==tissue, "tissue_abbrv"], nrow(df))
  df$cpg <- rep(gene_name, nrow(df))
  return(df)
}

names(traits_prev) <- traits
names(traits) <- traits_prev
get_gene_plot <- function(trait, tissue, gene_name){
  print(paste0(trait, " - ", tissue, " - ", gene_name))
  # Gene TPM --
  tpm <- cbind.data.frame(t(readRDS(paste0(project_path,"Tissues/", tissue, "/data.rds"))[gene_name,]))
  # tpm: matrix with gene TPM expression data for spleen samples used in this study. T
  
  # get gene data --
  data <- get_gene_data(tissue, tpm, gene_name)
  cols <- rep(traits_cols[trait], length(levels(data[,trait])))
  names(cols) <- as.character(levels(data[,trait]))
  
  # plot 1 --
  p1 <- ggplot(data = data,
               aes(x = eval(parse(text=trait)),
                   y = as.numeric(Beta),
                   fill =  eval(parse(text=trait))),
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
                 hjust = 0.5, vjust = 0.9, size = 3) +
    labs(title=paste0(gene_name, " (", tissue_info[tissue_info$tissue_ID==tissue, "tissue_abbrv"], ")")) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.title.y = ggtext::element_markdown(size=14),
          plot.title = element_text(hjust = 0.5,
                                    size = 15),
          strip.background = element_rect(fill="#B3B3B3"),
          legend.position = "none") 
  
  # plot 2 --
  d <- data.frame("variable" = trait,
                  "value"  = dma_10[[tissue]][[traits_prev[trait]]][dma_10[[tissue]][[traits_prev[trait]]]$cpg == gene_name, "R2"])
  variable_col <- traits_cols[trait]
  names(variable_col) <- trait
  
  p2 <- ggplot(d, aes(x=variable, y=value)) + 
    geom_bar(aes(fill = as.factor(variable)), stat = "identity") +
    scale_fill_manual(values = variable_col) + ylim(c(0,1))+
    coord_flip() + 
    xlab("") + ylab("DNA methylation variation explained (%)") +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5,
                                    size = 15),
          strip.background = element_rect(fill="#B3B3B3"),
          legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y =  element_blank()) 
  p <- ggarrange(p1, p2, nrow = 2, heights =  c(9,2))
  return(p)
}

gene_examples <- list(#c("Ancestry", "Ovary", "cg19649313"),
  #c("Ancestry", "MuscleSkeletal", "cg10306532"),
  c('Ancestry',"Lung",'cg01878807'),
  c('Ancestry',"Lung",'cg05171021')
  #c("Ancestry", "Lung", "cg16999677"),
  #c("Sex", "MuscleSkeletal", "cg19765154"),
  #c("Sex", "MuscleSkeletal", "cg02989351"),
  #c("Age", "Ovary", "cg05708550"))
  #c("Age", "ColonTransverse", "cg25940946"),
  #c("Age", "Testis", "cg13233461")
)
library(ggpubr)
plot_path <- '~/marenostrum/Projects/GTEx_v8/Methylation/Plots/'
for(i in 1:length(gene_examples)){
  p <- get_gene_plot(gene_examples[[i]][1], gene_examples[[i]][2], gene_examples[[i]][3])  
  ggexport(p,filename = paste0(plot_path, "Figure_2E.", gene_examples[[i]][1], "_", gene_examples[[i]][2], "_", gene_examples[[i]][3], ".v2.pdf"),
           width = 3, height = 4)
}

