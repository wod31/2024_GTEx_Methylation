library(ComplexHeatmap)
library(RColorBrewer)
suppressPackageStartupMessages(library(circlize))

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

tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Methylation/Data/Tissue_info_whole.rds"))

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
  admixture_ancestry <- read.table('~/marenostrum/Projects/GTEx_v8/Methylation/Data/admixture_inferred_ancestry.txt')
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


### plot rho of correlated probes by location ############
### first tissue sep ####
### read Correlations
get_corr <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "SEX"){
    NA
  }else{
    model <- readRDS(paste0(project_path, "Tissues/",tissue, "/",trait,'_Correlations_probes_genes.rds'))
    #rownames(model[[trait]][model[[trait]]$adj.P.Val<0.05,])
    model
  }
}
DMPs_cor <- lapply(c('EURv1','SEX','AGE','BMI'), function(trait) lapply(tissues, function(tissue) get_corr(tissue, trait)))
names(DMPs_cor) <- c("Ancestry", "Sex", "Age", "BMI")
for(trait in c("Ancestry", "Sex", "Age", "BMI")){names(DMPs_cor[[trait]]) <- tissues}

all_cor <- unlist(DMPs_cor, recursive = FALSE)
all_cor_df <- do.call("rbind", all_cor)
all_cor_df[,c('Trait','Tissue','Gene')] <- str_split_fixed(rownames(all_cor_df), '\\.', 3)
write.table(all_cor_df, '~/marenostrum/Projects/GTEx_v8/Methylation/Data/Correlations_all_traits_tissues.txt', sep = '\t', quote = F, 
            row.names = F, col.names = T)
all_cor_df <- read.delim('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Correlations_all_traits_tissues.txt', sep = '\t')

### density
names(tissues_cols) <- tissue_info$tissue_ID
all_cor_df <- all_cor_df[!is.na(all_cor_df$class),]
all_cor_df$class <- factor(all_cor_df$class, levels = c("promoter", "enhancer", "gene_body"))
all_cor_df$class <- droplevels(all_cor_df$class)
ggplot(all_cor_df[all_cor_df$p.adj<0.05,], aes(x=abs(cor), color=Tissue))+
  geom_density()+facet_grid(class ~ Trait,scales = "free_y") +
  scale_color_manual(values=tissues_cols)+
  theme(axis.title.y = element_text(margin = margin(r = 2), size = 11),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.text.x = element_text(size = 9, colour = "black", angle = 90, vjust = 0.5)) + theme_classic()


#### numbers of correlations #####
traits_cols <- c('#C49122','#4B8C61','#70A0DF','#A76595')
names(traits_cols) <- c("Ancestry", "Sex", "Age", "BMI")
strip <- strip_themed(background_x = elem_list_rect(fill = traits_cols))
all_cor_df$Trait <- factor(all_cor_df$Trait, levels = c("Ancestry", "Sex", "Age", "BMI"))

ggplot(all_cor_df, aes(y = Tissue, fill=class)) +
  geom_bar(position = 'fill', alpha=0.8) + 
  geom_text(aes(label=after_stat(count), x = after_stat(count+2)), stat='count', position='fill', size=3,hjust=1.6, angle=20) +
  theme_bw() + ylab('Tissue') + xlab('Proportion correlated probes') +
  facet_wrap2(~ Trait, strip = strip, nrow = 1) + theme_classic() + scale_x_continuous(breaks=seq(0, 1, 0.5))

ggplot(all_cor_df, aes(y = Tissue, fill=class)) +
  geom_bar(position = 'stack', alpha=0.8) + 
  #geom_text(aes(label=after_stat(count), x = after_stat(count+1.5)), stat='count', position='stack', size=3,hjust=1, angle=45) +
  theme_bw() + ylab('Tissue') + xlab('Proportion correlated probes') +
  facet_wrap2(~ Trait, strip = strip, nrow = 1) + theme_classic() #+ scale_x_continuous(breaks=seq(0, 1, 0.5))


#### median rho per tissue separated by trait and location ######
grouped <- all_cor_df %>%
  group_by(Trait, Tissue, class)%>% 
  summarise(Mean=mean(abs(cor)), Max=max(abs(cor)), Min=min(abs(cor)), Median=median(abs(cor)), Std=sd(abs(cor)))
head(grouped)

grouped$Trait <- factor(grouped$Trait, levels = c("Ancestry", "Sex", "Age", "BMI"))
ggplot(grouped, aes(class, Median, color=Tissue)) + geom_jitter() + 
  geom_boxplot(outlier.shape = NA,
               notch = T,
               width = 0.25, color = 'grey', alpha=0.5) + theme_classic() + xlab("") + ylab("abs(rho)") +
  facet_wrap2(~ Trait, strip = strip, nrow = 1) +
  scale_color_manual(values=tissues_cols)

### tissue sharing #####
grouped <- all_cor_df %>%
  group_by(Trait, gene, probe)%>% 
  summarise(N=n()) 

cpgs_gene <- read.delim('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Methylation_Epic_gene_promoter_enhancer_processed.txt')

### enrichment correlated positions
library(missMethyl)
my_fisher <- function(trait){
  #             DS      Not DS
  # type
  # other_types
  print(trait)
  res <- grouped[grouped$Trait == trait,]
  chrom_tissue <- cpgs_gene
  
  if (nrow(res) < 1 ) {
    return(NA)
  }
  tryCatch(
    {res_2 <- gometh(unique(res$probe), all.cpg=unique(chrom_tissue$IlmnID),
                     collection="GO", array.type="EPIC")
    res_2 <- res_2#[res$ONTOLOGY=="BP",]
    print(table(res_2$FDR<0.05))
    return(res_2)
    },  error=function(cond) {
      message("Error")
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    })
}
# Two-tailed Fisher test
#families <- as.vector(unique(shared_cpgs$region_chromhmm))
fisher_results <- lapply(c("Ancestry" ,"Sex" , "Age"  , "BMI"), function(trait)  my_fisher(trait))
names(fisher_results) <- c("Ancestry" ,"Sex" , "Age"  , "BMI")

grouped_2 <- all_cor_df %>%
  group_by(Trait,Tissue, probe)%>% 
  summarise(N=n()) 

cpgs_gene <- read.delim('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Methylation_Epic_gene_promoter_enhancer_processed.txt')

### enrichment correlated positions
library(missMethyl)
my_fisher <- function(trait, tissue){
  #             DS      Not DS
  # type
  # other_types
  print(trait)
  res <- grouped_2[grouped_2$Trait == trait & grouped_2$Tissue == tissue,]
  chrom_tissue <- cpgs_gene
  
  if (nrow(res) < 1 ) {
    return(NA)
  }
  tryCatch(
    {res_2 <- gometh(unique(res$probe), all.cpg=unique(chrom_tissue$IlmnID),
                     collection="GO", array.type="EPIC")
    res_2 <- res_2#[res$ONTOLOGY=="BP",]
    print(table(res_2$FDR<0.05))
    return(res_2)
    },  error=function(cond) {
      message("Error")
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    })
}
# Two-tailed Fisher test
#families <- as.vector(unique(shared_cpgs$region_chromhmm))
fisher_results <- lapply(c("Ancestry" ,"Sex" , "Age"  , "BMI"), function(trait)  lapply(tissues, function(tissue) my_fisher(trait, tissue)))
names(fisher_results) <- c("Ancestry" ,"Sex" , "Age"  , "BMI")

for (trait in c("Ancestry" ,"Sex" , "Age"  , "BMI")) {
  names(fisher_results[[trait]]) <- tissues
}
