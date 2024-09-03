#### overlap varpart methylation and expression ####
### read expression data ###
library(variancePartition)
library(tidyr)
library(dplyr)
first_dir <- "~/"
setwd(paste0(first_dir, "marenostrum/Projects/GTEx_v8/Methylation/"))

# -------------- #
print(Sys.time())
#-------------- #

tissues <- list.dirs("Tissues/", full.names = F)[-1]
tissues <- tissues[-grep('Old',tissues)]

#### reading varPart values ####
beta <- readRDS(paste0("var_part_expression.rds"))
beta$gene <- gsub('\\..*','',rownames(beta))

library(reshape2)
betas_locations_m <- melt(data = beta[,c("gene", "Donor", "Tissue")], id.vars = c('gene'),
                          variable.name = 'Variable', value.name = 'Proportion')
head(betas_locations_m)

### read methylation 
chuncks <- c(1:16)

#### reading varPart values ####
beta_m <- lapply(chuncks, function(chnk) readRDS(paste0('varPart/', chnk, "_chunck_var_part.rds")))
beta_df <- do.call("rbind",beta_m)

## add anotation
Sys.time()
data_path <- "~/marenostrum/Projects/GTEx_v8/Methylation/Data/"
annotation <- read.delim(paste0(data_path, "Methylation_Epic_gene_promoter_enhancer_processed.txt"), sep = '\t', header = T)
Sys.time()
beta_df$cpg <- rownames(beta_df)
# beta_df$Type <- 'Other'
# beta_df$Type[rownames(beta_df) %in% annotation$IlmnID[annotation$Type=="Promoter_Associated"]] <- "Promoter_Associated"
# beta_df$Type[rownames(beta_df) %in% annotation$IlmnID[annotation$Type=="Enhancer_Associated"]] <- "Enhancer_Associated"
# beta_df$Type[rownames(beta_df) %in% annotation$IlmnID[annotation$Type=="Gene_Associated"]] <-"Gene_Associated"
# table(beta_df$Type)

beta_df <- merge(beta_df, annotation[,c("IlmnID","UCSC_RefGene_Name","Type")], all.x=TRUE, by.x='cpg', by.y='IlmnID')
beta_df <- beta_df %>% distinct()

head(beta_df)

### genes_high tissue sp 
genes <- betas_locations_m$gene[betas_locations_m$Variable=='Donor' & betas_locations_m$Proportion>0.5]
genes_from_meth <- unique(beta_df$UCSC_RefGene_Name[beta_df$SUBJID>0.5])

gene_annotation <- read.delim("~/marenostrum/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/data/public/gencode.v26.GRCh38.genes.bed", header=F)[,c(6,7)]
#gene_annotation <- read.delim("~/marenostrum/MN4/bsc83/Projects/GTEx_v8/Jose/04_Smoking/github/analysis/data/public/gencode.v26.GRCh38.genes.bed", header=F)[,c(6,7)]

colnames(gene_annotation) <- c("gene", "symbol")

#These were genes duplicated, I changed their names to their correct one
gene_annotation$symbol[gene_annotation$gene=="ENSG00000253972.5"] <- "MAL2-AS1" #I found this a posteriori
gene_annotation$symbol[gene_annotation$gene=="ENSG00000283992.1"] <- "SLURP2" #Insted of LYNX1
gene_annotation$symbol[gene_annotation$gene=="ENSG00000235271.5"] <- "GC22P027299"
gene_annotation$symbol[gene_annotation$gene=="ENSG00000229694.6"] <- "C9orf73"
gene_annotation$symbol[gene_annotation$gene=="ENSG00000228741.2"] <- "GC13P024553" 

genes_symbol <- gene_annotation$symbol[gsub('\\..*','',gene_annotation$gene) %in% genes]

sum(genes_symbol %in% genes_from_meth)

### get varpart values for highly variable genes at meth level  #####
head(betas_locations_m)
gene_annotation$gene <- gsub('\\..*','',gene_annotation$gene)
betas_locations_m <- merge(betas_locations_m, gene_annotation)
betas_locations_m$var <- 'Non-ind-variable'
betas_locations_m$var[betas_locations_m$symbol %in% genes_from_meth] <- 'Ind-variable'

#png('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/tissue_sharing_ancestry_violin.png', width = 900, height = 700, res = 300)
colors <- c('#A39A92','#9C731C')
library(ggpubr)
p <- (ggplot(betas_locations_m[betas_locations_m$Variable=='Donor',], aes(var, Proportion)) + #geom_jitter(aes(Type, number, col=Type), alpha=0.5) +
        theme_bw() + #geom_violin(aes(x = type, y = number, fill=type),alpha=0.8) +
        geom_boxplot(aes(var, Proportion),col = "black",
                     fill = "white",
                     outlier.shape = NA,
                     notch = T,
                     width = 0.25) +
        geom_violin(aes(var, Proportion,fill=var),col = "black")+
        scale_fill_manual(values = rev(colors)) +
        scale_color_manual(values = rev(colors)) + xlab('') + ylab('Varience explained') +
        stat_compare_means(label.y = 0.6,size=2)+
        scale_x_discrete(breaks=c("Non-ind-variable",'Ind-variable'),
                         labels=c("Non-ind-\nvariable",'Ind-\nvariable'))+
        #ggrepel::geom_text_repel(aes(type, number,label=ifelse(label=='Yes',as.character(UCSC_RefGene_Name),'')), max.overlaps = Inf, )+
        theme(legend.position = "none",
              axis.text = element_text(size = 9),
              text = element_text(size = 6),
              axis.title = element_text(size = 9)))
#dev.off()
ggsave('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/expression_ind_variable_meth.pdf', p, dpi = 300, width = 2, height = 2, units = 'in')

mu <- betas_locations_m %>% 
  group_by(Variable , var) %>%
  summarise(grp.mean = mean(Proportion), grp.median=median(Proportion))
mu

#### tissue specific ####
#betas_locations_m <- merge(betas_locations_m, gene_annotation)
genes_from_meth_t <- unique(beta_df$UCSC_RefGene_Name[beta_df$Tissue>0.5])
betas_locations_m$var_t <- 'Non-tissue-variable'
betas_locations_m$var_t[betas_locations_m$symbol %in% genes_from_meth_t] <- 'Tissue-variable'

#png('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/tissue_sharing_ancestry_violin.png', width = 900, height = 700, res = 300)
colors <- c('#A39A92','#9C731C')
library(ggpubr)
p <- (ggplot(betas_locations_m[betas_locations_m$Variable=='Tissue',], aes(var_t, Proportion)) + #geom_jitter(aes(Type, number, col=Type), alpha=0.5) +
        theme_bw() + #geom_violin(aes(x = type, y = number, fill=type),alpha=0.8) +
        geom_boxplot(aes(var_t, Proportion),col = "black",
                     fill = "white",
                     outlier.shape = NA,
                     notch = T,
                     width = 0.25) +
        geom_violin(aes(var_t, Proportion,fill=var_t),col = "black")+
        scale_fill_manual(values = (colors)) +
        scale_color_manual(values = (colors)) + xlab('') + ylab('Varience explained') +
        stat_compare_means(label.y = 0.6,size=2)+
        scale_x_discrete(breaks=c("Non-tissue-variable",'Tissue-variable'),
                         labels=c("Non-tissue-\nvariable",'Tissue-\nvariable'))+
        #ggrepel::geom_text_repel(aes(type, number,label=ifelse(label=='Yes',as.character(UCSC_RefGene_Name),'')), max.overlaps = Inf, )+
        theme(legend.position = "none",
              axis.text = element_text(size = 9),
              text = element_text(size = 6),
              axis.title = element_text(size = 9)))
#dev.off()
ggsave('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/expression_tissue_variable_meth.pdf', p, dpi = 300, width = 2, height = 2, units = 'in')

mu <- betas_locations_m %>% 
    group_by(Variable , var_t) %>%
  summarise(grp.mean = mean(Proportion), grp.median=median(Proportion))
mu

### now do it per region #####
head(betas_locations_m)
head(beta_df)

genes_prom <- unique(beta_df$UCSC_RefGene_Name[beta_df$SUBJID>0.5 & beta_df$Type == 'Promoter_Associated'])
genes_enh <- unique(beta_df$UCSC_RefGene_Name[beta_df$SUBJID>0.5 & beta_df$Type == 'Enhancer_Associated'])
genes_gene <- unique(beta_df$UCSC_RefGene_Name[beta_df$SUBJID>0.5 & beta_df$Type == 'Gene_Associated'])

betas_locations_m$type <- 'Non-ind-variable'
betas_locations_m$type[betas_locations_m$symbol %in% genes_gene] <- 'Ind-variable-Gene'
betas_locations_m$type[betas_locations_m$symbol %in% genes_enh] <- 'Ind-variable-Enh'
betas_locations_m$type[betas_locations_m$symbol %in% genes_prom] <- 'Ind-variable-Prom'

colors <- c('#A39A92','#9C731C')
library(ggpubr)
for (type_m in c('Ind-variable-Gene','Ind-variable-Enh','Ind-variable-Prom')) {
  p <- (ggplot(betas_locations_m[betas_locations_m$Variable=='Donor' & betas_locations_m$type %in% c(type_m,'Non-ind-variable'),], aes(type, Proportion)) + #geom_jitter(aes(Type, number, col=Type), alpha=0.5) +
          theme_bw() + #geom_violin(aes(x = type, y = number, fill=type),alpha=0.8) +
          geom_boxplot(aes(type, Proportion),col = "black",
                       fill = "white",
                       outlier.shape = NA,
                       notch = T,
                       width = 0.25) +
          geom_violin(aes(type, Proportion,fill=type),col = "black")+
          scale_fill_manual(values = rev(colors)) +
          scale_color_manual(values = rev(colors)) + xlab('') + ylab('Varience explained') +
          stat_compare_means(label.y = 0.6,size=2)+
          scale_x_discrete(breaks=c("Non-ind-variable",type_m),
                           labels=c("Non-ind-\nvariable",type_m))+
          #ggrepel::geom_text_repel(aes(type, number,label=ifelse(label=='Yes',as.character(UCSC_RefGene_Name),'')), max.overlaps = Inf, )+
          theme(legend.position = "none",
                axis.text = element_text(size = 9),
                text = element_text(size = 6),
                axis.title = element_text(size = 9)))
  #dev.off()
  ggsave(paste0('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/expression_',type_m,'_meth.pdf'), p, dpi = 300, width = 2, height = 2, units = 'in')
  
  mu <- betas_locations_m[betas_locations_m$type %in% c(type_m,'Non-ind-variable'),] %>% 
    group_by(Variable , type) %>%
    summarise(grp.mean = mean(Proportion), grp.median=median(Proportion), n=n())
  print(mu)
  
  
}

#### plot corr gene meth #####
#### reading metadata ####
### example gene ind variable and not ####
# ENSG00000198502 cg00440797
# ENSG00000157470 FAM81A cg00000905 
tissues <- list.dirs("Tissues/", full.names = F)[-1]
tissues <- tissues[-grep('Old',tissues)]
tissues <- c('Lung')

metadata <- lapply(tissues, function(tissue) readRDS(paste0("Tissues/", tissue, "/metadata.rds")))
names(metadata) <- tissues

#### reading beta values ####
beta <- lapply(tissues, function(tissue) readRDS(paste0("Tissues/", tissue, "/data.rds"))[c('cg00440797','cg00000905'),])
names(beta) <- tissues

for (tissue in tissues) {
  metadata[[tissue]][,'Tissue'] <- tissue
}

### reading admixture and joining metadata #### 
metadata_df <- do.call("rbind",lapply(tissues, function(tissue) metadata[[tissue]][,c("SUBJID", "Tissue")]))
head(metadata_df)
rownames(metadata_df) <- paste0(metadata_df$SUBJID,':',metadata_df$Tissue)

### create label with both subj-id pairs ####
for (tissue in tissues) {
  colnames(beta[[tissue]]) <- paste0(colnames(beta[[tissue]]),':',tissue)
}

### merging all data together ####
beta_df <- do.call("cbind",lapply(tissues, function(tissue) beta[[tissue]]))
beta_df <- as.data.frame(t(beta_df))
beta_df$Tissue <- gsub('.*:','',rownames(beta_df))
head(beta_df)

metadata_df$label <- rownames(metadata_df)
beta_df$label <- rownames(beta_df)

### expression ###
metadata <- lapply(tissues, function(tissue) readRDS(paste0("Tissues/", tissue, "/metadata_expression.rds")))
names(metadata) <- tissues
for (tissue in tissues) {
  metadata[[tissue]][,'Tissue'] <- tissue
}

metadata_df_expr <- do.call("rbind",lapply(tissues, function(tissue) metadata[[tissue]][,c("Donor", "Tissue")]))
head(metadata_df_expr)
rownames(metadata_df_expr) <- paste0(metadata_df_expr$Donor,':',metadata_df_expr$Tissue)
metadata_df_expr$label <- rownames(metadata_df_expr)

### subset samples ####
metadata_df <- metadata_df[metadata_df$label %in% metadata_df_expr$label,]
metadata_df_expr <- metadata_df_expr[metadata_df$label,]
beta_df <- beta_df[metadata_df$label,]

## read expression ####
#### reading beta values ####
counts <- lapply(tissues, function(tissue) readRDS(paste0("Tissues/", tissue, "/counts.rds"))[c('ENSG00000198502.5','ENSG00000157470.11'),])
names(counts) <- tissues

# for (tissue in tissues) {
#   counts[[tissue]][,'Tissue'] <- tissue
# }

for (tissue in tissues) {
  colnames(counts[[tissue]]) <- paste0(gsub('-....$','',colnames(counts[[tissue]])),':',tissue)
}

### merging all data together ####
counts_df <- do.call("cbind",lapply(tissues, function(tissue) counts[[tissue]]))
counts_df <- as.data.frame(t(counts_df))
counts_df$Tissue <- gsub('.*:','',rownames(counts_df))
head(counts_df)
counts_df$label <- rownames(counts_df)

##3 subset melt and merge ####
counts_df <- counts_df[metadata_df_expr$label,]
colnames(counts_df) <- c('Expr','Tissue','label')
colnames(beta_df) <- c('Meth','Tissue','label')

all <- merge(beta_df, counts_df, by = 'label')
all <- all %>% distinct()
all$Meth <- as.numeric(all$Meth)
all$Expr <- as.numeric(all$Expr)
all$log10Expr <- log10(all$Expr) 

library("ggpubr")
pdf('Plots/corr_expr_meth_hla_DRB5.pdf', height = 3, width = 3)
ggscatter(all, x = "log10Expr", y = "Meth", color = "Tissue.x",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "log10(Expression)", ylab = "Methylation", label.y = 0.9)
dev.off()

### plot genes ####
## read data ##
genes <- read.delim('~/marenostrum/Projects/GTEx_v8/Methylation/Data/to_plot_genes_varpart.txt')

get_box_stats <- function(y, upper_limit = max(y) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "n =", length(y)
    )
  ))
}

p1 <- ggplot(data = genes,
             aes(x = Gene,
                 y = as.numeric(counts),
                 fill =  Gene),
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
  ylab("Counts") +
  #scale_fill_manual(values = cols) +
  stat_summary(fun.data = get_box_stats, geom = "text",
               hjust = 0.5, vjust = 0.9, size = 3) +
  #labs(title=paste0(gene_name, " (", tissue_info[tissue_info$tissue_ID==tissue, "tissue_abbrv"], ")")) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.title.y = element_text(size=14),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none") 
ggexport(p1,filename = paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Plots/Figure_1SE.example_genes_lung.counts.pdf"),
         width = 3, height = 4)

#### plot top expression ####
#### name genes in immune pathways #####
mhc <- read.table('Data/GO_MHC_Complex_assembly.txt', sep = '\t', row.names = NULL)
peptide <- read.table('Data/GO_antigen_processing.txt', sep = '\t', row.names = NULL)
beta <- merge(beta, gene_annotation)
beta$immune <- 'No'
beta$immune[beta$symbol %in% c(mhc$V2, peptide$V2)] <- 'Immune'

beta$name <- 'No'
beta$name[beta$Donor > 0.1 & beta$immune=='Immune'] <- 'Immune'

library(ggrepel)
pdf(paste0("Plots/","Imune_genes_ind_variable.dotplot.pdf"), width = 6, height = 5)
ggplot(beta, aes(x=Tissue, y=Donor)) +
  geom_point(alpha=0.6, colour = '#BFC0C0') +
  geom_point(data = beta[beta$immune == 'Immune',], aes(x=Tissue, y=Donor),alpha=0.8, colour='#058ED9') +
  geom_text_repel(aes(label=ifelse(name=='Immune',as.character(symbol),'')), max.overlaps = Inf) +
  theme_classic()
dev.off()

p <- (ggplot(beta, aes(immune, Donor)) + #geom_jitter(aes(Type, number, col=Type), alpha=0.5) +
        theme_bw() + #geom_violin(aes(x = type, y = number, fill=type),alpha=0.8) +
        geom_boxplot(aes(immune, Donor),col = "black",
                     fill = "white",
                     outlier.shape = NA,
                     notch = T,
                     width = 0.25) +
        geom_violin(aes(immune, Donor,fill=immune),col = "black")+
        scale_fill_manual(values = rev(colors)) +
        scale_color_manual(values = rev(colors)) + xlab('') + ylab('Varience explained') +
        stat_compare_means(label.y = 0.6,size=2)+
        # scale_x_discrete(breaks=c("Non-ind-variable",type_m),
        #                  labels=c("Non-ind-\nvariable",type_m))+
        #ggrepel::geom_text_repel(aes(type, number,label=ifelse(label=='Yes',as.character(UCSC_RefGene_Name),'')), max.overlaps = Inf, )+
        theme(legend.position = "none",
              axis.text = element_text(size = 9),
              text = element_text(size = 6),
              axis.title = element_text(size = 9)))
#dev.off()
ggsave(paste0('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/expression_','immune_donor.pdf'), p, dpi = 300, width = 2, height = 2, units = 'in')

mu <- beta %>% 
  group_by(immune) %>%
  summarise(grp.mean = mean(Donor), grp.median=median(Donor), n=n())
print(mu)
