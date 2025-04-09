#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Plot fisher enrichments of DMPs per chromatin state and tissue; for sex-age specific DMPs
# @software version: R=4.2.2

#### plot enrichments
library(ggplot2)
library(ggh4x)

### read enrichments
hallmark <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Hallmarks_enrichments_DMP.hyper.rds')
hallmark_hypo<- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Hallmarks_enrichments_DMP.hypo.rds')
KEGG <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/KEGG_enrichments_DMP.rds')
KEGG <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/KEGG_enrichments_DMP.hyper.rds')
KEGG_hypo <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/KEGG_enrichments_DMP.hypo.rds')
GO <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/GO_enrichments_DMP_peers5.hyper.continous.rds')
GO <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/GO_enrichments_DMP_peers5.continous.rds')
GO_hypo <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/GO_enrichments_DMP_peers5.hypo.continous.rds')

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Age", "Ancestry", "BMI", "Sex")
traits_to_use <- c('EURv1','SEX2','AGE','BMI')

### transform data for plotting

##all final table
final_table_all <- data.frame(name="X",term="X", GeneRatio=1, FDR=1, trait=1, tissue=1)
for (tissue in tissues) {
  print(tissue)
  probes <- data.frame(name="X",term="X",GeneRatio=1, FDR=1, trait=1, tissue=1)
  for (trait in traits_to_use) {
    res <- as.data.frame(GO[[tissue]][[trait]])
    if (nrow(res)<1) {
      next}
    res$name <- rownames(res)
    ###for GO 
    res$term <- res$TERM
    res$tissue <- tissue
    res$trait <- trait
    res$GeneRatio <- res$DE/res$N
    probes <- rbind(probes, res[,c('name','term','GeneRatio','FDR','trait','tissue')])
  }
  probes <- probes[-1,]
  final_table_all <- rbind(final_table_all, probes)
}

final_table_all <- final_table_all[-1,]

head(final_table_all)

### all traits hyper
final_table_hyper <- data.frame(name="X",term="X", GeneRatio=1, FDR=1, trait=1, tissue=1)
for (tissue in tissues) {
  print(tissue)
  probes <- data.frame(name="X",term="X",GeneRatio=1, FDR=1, trait=1, tissue=1)
  for (trait in traits_to_use) {
    res <- as.data.frame(GO[[tissue]][[trait]])
    if (nrow(res)<1) {
      next}
    res$name <- rownames(res)
    ###for GO 
    res$term <- res$TERM
    res$tissue <- tissue
    res$trait <- trait
    res$GeneRatio <- res$DE/res$N
    probes <- rbind(probes, res[,c('name','term','GeneRatio','FDR','trait','tissue')])
  }
  probes <- probes[-1,]
  final_table_hyper <- rbind(final_table_hyper, probes)
}

final_table_hyper <- final_table_hyper[-1,]

head(final_table_hyper)
### all traits hypo
final_table_hypo <- data.frame(name="X",term="X",GeneRatio=1, FDR=1, trait=1, tissue=1)
for (tissue in tissues) {
  print(tissue)
  probes <- data.frame(name="X",term="X",GeneRatio=1, FDR=1, trait=1, tissue=1)
  for (trait in traits_to_use) {
    res <- as.data.frame(GO_hypo[[tissue]][[trait]])
    if (nrow(res)<1) {
      next}
    res$name <- rownames(res)
    ###for GO 
    res$term <- res$TERM
    res$tissue <- tissue
    res$trait <- trait
    res$GeneRatio <- res$DE/res$N
    probes <- rbind(probes, res[,c('name','term','GeneRatio','FDR','trait','tissue')])
  }
  probes <- probes[-1,]
  final_table_hypo <- rbind(final_table_hypo, probes)
}

final_table_hypo <- final_table_hypo[-1,]

head(final_table_hypo)

#
# #### Bubble plot ####
#
# pdf('/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Tissues/')
counts <- as.data.frame(table(final_table_all$term,final_table_all$trait))
counts$Var3 <- 'all'
table(counts$Var2)

traits_cols <- c('#70A0DF','#C49122')
names(traits_cols) <- c('AGE', 'EURv1')

pal <- rev(c('#0061FF', '#1872FF', '#3183FF', '#4994FF', '#61A5FF', '#7AB7FF', '#92C8FF', '#AAD9FF', '#C3EAFF', '#DBFBFF'))

df <- final_table_all[final_table_all$term %in% counts$Var1[counts$Var2=='AGE' & counts$Freq>3] & final_table_all$trait=='AGE',]

pdf('marenostrum/Projects/GTEx_v8/Methylation/Plots/Enrichment_GO_age_shared.pdf', width = 6, height = 5)
ggplot(df, aes(y=term, x=tissue, size=-log10(FDR), fill=GeneRatio)) +
  geom_point(alpha=0.8, shape=21) +
  scale_size(range = c(0, 15), name="-log10(adjPvalue)") +
  scale_fill_gradientn(colours = pal) + 
  #scale_fill_gradient2(low="yellow", high="red", name='GeneRatio', midpoint = 0.5)+
  #theme_ipsum() +
  theme(legend.position="right") +
  #coord_flip()+
  ylab("") +
  xlab("") +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        panel.grid = element_line(colour = 'light grey'), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
dev.off()


### sharing pathways 
final_table_hypo$dir <- 'Hypomethylated'
final_table_hyper$dir <- 'Hypermethylated'
all <- rbind(final_table_hypo, final_table_hyper)
head(all)

counts <- as.data.frame(table(all$term,all$dir,all$trait))

traits_cols <- c('#70A0DF','#4B8C61')
names(traits_cols) <- c('AGE','SEX2')
strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols))

library(ggplot2)
library(ggh4x)
pdf('~/marenostrum/Projects/GTEx_v8/Methylation/Plots/Enrichment_Hyper_Hypo_barplot.pdf', 
    width = 6, height = 4)
ggplot(counts[counts$Freq>0,], aes(x = Var2, fill = as.factor(Freq))) +
  geom_bar() + ylab('Nº of Terms') + xlab('') +
  geom_text(aes(label=after_stat(count), y = after_stat(count+15)), stat='count', position='stack') +
  labs(fill='Nº of Tissues') + scale_fill_grey(start = 0.9, end = 0) + theme_bw()+
  facet_wrap2(~ Var3, strip = strip, nrow = 1) 
dev.off()

#all
counts <- as.data.frame(table(final_table_all$term,final_table_all$trait))
counts$Var3 <- 'all'
table(counts$Var2)

traits_cols <- c('#70A0DF','#C49122')
names(traits_cols) <- c('AGE', 'EURv1')
strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols))

counts$Var2 <- factor(counts$Var2, levels = c('AGE','EURv1'))

ggplot(counts[counts$Freq>0,], aes(x = Var3, fill = as.factor(Freq))) +
  geom_bar() + ylab('Nº of Terms') + xlab('') +
  geom_text(aes(label=after_stat(count), y = after_stat(count+15)), stat='count', position='stack') +
  labs(fill='Nº of Tissues') + scale_fill_grey(start = 0.9, end = 0) + theme_bw() + 
  facet_wrap2(~ Var2, strip = strip, nrow = 1) 


##### plot example of enrichments shared in ageing #####
head(all)
head(counts)
shared <- all[all$term %in% counts$Var1[counts$Freq == 4],]
hypo_shared <- all[all$term %in% counts$Var1[counts$Freq == 3 & counts$Var2 == 'Hypomethylated'],]
merged <- rbind(shared, hypo_shared)

ggplot(merged, aes(y=tissue, x=term, size=-log10(FDR), fill=dir)) +
  geom_point(alpha=0.8, shape=21) +
  scale_size(range = c(0, 21), name="-log10(adjPvalue)") + 
  #scale_fill_gradient2(low="yellow", high="red", name='GeneRatio', midpoint = 0.5)+
  #theme_ipsum() +
  theme(legend.position="right") +
  ylab("") +
  xlab("") +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        panel.grid = element_line(colour = 'light grey'), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

