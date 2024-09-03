#!/usr/bin/env Rscript

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
# first_dir <- "~/Documents/mn4/"
first_dir <- "/gpfs/"
project_path <- paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Methylation/")

start_time <- Sys.time()
print("Reading data")
# Sys.time()
# data <- readRDS(paste0(project_path,'Tissues/', tissue, "/data.rds")) #From whole compressed data in 5.6G to compressed 1.4G/1.1Gb only in Lung (the highest number of samples)
# Sys.time() #12 minutes to load 15 Gb
# beta <- data

print("Reading metadata")
metadata <- readRDS(paste0(project_path, 'Tissues/', tissue, "/metadata.rds"))
# metadata$Ancestry <- as.factor(metadata$Ancestry)
# if(sum(metadata$Ancestry=="AMR")<5){
#   print("Not enough AMRs")
#   metadata <- metadata[metadata$Ancestry!="AMR",]
#   #beta <- beta[,metadata$SUBJID]
# }
# 
# print("Filtering ASN")
# metadata <- metadata[metadata$Ancestry!="ASN",]
# metadata$Ancestry <- droplevels(metadata$Ancestry)
# #beta <- beta[,colnames(beta) %in% metadata$SUBJID]

print("Reading Admixture results")
admixture_ancestry <- read.table(paste0(project_path,'Data/admixture_inferred_ancestry.txt'))
colnames(admixture_ancestry) <- c('SUBJID','AFRv1','EURv1','inferred_ancestry','AFRv2','EURv2')
metadata <- merge(metadata, admixture_ancestry[,c("SUBJID","EURv1")], by='SUBJID')

metadata$SEX <- as.factor(metadata$SEX)
if(length(levels(metadata$SEX))==1){
  print("Sexual tissue")
  metadata <- metadata[,-which(names(metadata) == "SEX")]
  individual_variables <- c("EURv1", "AGE", "BMI")
} else{
  individual_variables <- c("EURv1", "SEX", "AGE", "BMI")
}
metadata$DTHHRDY <- as.factor(metadata$DTHHRDY)
rownames(metadata) <- metadata$SUBJID
metadata$SUBJID <- NULL
# if(sum(metadata$Ancestry=="AFR")<5){
#   print("Not enough AFRs")
#   metadata <- metadata[, !colnames(metadata) %in% c("Ancestry")]
#   individual_variables <- individual_variables[!individual_variables %in% "Ancestry"]
# }

### make sure order is the same
#beta <- beta[,rownames(metadata)]

#probes <- rownames(beta)

metadata_2 <- metadata[,c(individual_variables)]

# Model cpg as:
# meth_res ~ ieVariants_in_eGene + traits
# The genotypes of the sVariants are modelled as a factor, where: 0|0 -> 0, 0|1 -> 1; 1|1 ->2

print("metadata is prepared")

print("Reading mQTL data")
# eGpGs and independent mQTLs ----
inpath_mqtls <- "/gpfs/scratch/bsc83/MN4/bsc83/bsc83535/GTEx/v9/mQTLs/"
#inpath_eqtls <- "~/GTEx_v8_data/cisQTLs/"
mgene_data <- as.data.frame(data.table::fread(paste0(inpath_mqtls,tissue,".mQTLs.conditional.txt.gz")))
mgene_data <- mgene_data[mgene_data$V7<0.05 & abs(mgene_data$V3)<250000,]
mCpGs <- unique(gsub(':.*','',mgene_data$V1[mgene_data$V7<0.05 & abs(mgene_data$V3)<250000]))

## read data from p-nominal cpgs ####3
mgene_data_extra <- as.data.frame(data.table::fread(paste0(inpath_mqtls,tissue,".regular.perm.fdr.txt")))
mgene_data_extra <- mgene_data_extra[!mgene_data_extra$cpg_id %in% mCpGs,]
mgene_data_extra <- mgene_data_extra[mgene_data_extra$qval<0.05,]
mCpGs <- c(mCpGs,mgene_data_extra$cpg_id)

#### DMAnalysis
dea_res <- readRDS(paste0(project_path,'Tissues/', tissue,"/DML_results_5_PEERs_continous.rds"))[["EURv1"]]

# Save summary data ----
d <- list()
d[["mQTLs:CpG"]] <- length(mCpGs)
d[["Ancestry:DMP"]] <- sum(dea_res$adj.P.Val < 0.05)
d[["Ancestry:DMP:eGpG"]] <- sum(rownames(dea_res[dea_res$adj.P.Val < 0.05,])  %in% mCpGs)

# Print to std.out summary data ----
print(paste0("No. of mQTLs:CpG: ", length(mCpGs)))
print(paste0("No. of Ancestry:DMP:  ", sum(dea_res$adj.P.Val < 0.05)))
print(paste0("Ancestry:DMP:eGpG:  ", sum(rownames(dea_res[dea_res$adj.P.Val < 0.05,])  %in% mCpGs)))

# Select differentially methylated cpgs with at least 1 independent mQTL ----
dea_res <- dea_res[dea_res$adj.P.Val < 0.05,] # significant CpGs
dea_res <- dea_res[rownames(dea_res) %in% mCpGs,] # i-mCpG

# Subset mQTLs in eCpGs differentially methylated with at least 1 independent mQTL ----
mgene_data <- mgene_data[gsub(':.*','',mgene_data$V1) %in% rownames(dea_res),]
mgene_data_extra <- mgene_data_extra[mgene_data_extra$cpg_id %in% rownames(dea_res),]

# eVariants in eCpGs differentially methylated with at leats 1 independent mQTL ----
mVariants <- c(unique(mgene_data$V2),mgene_data_extra$variant_id) # imQTLs
print(paste0("No. of mVariants (i-mQTL) in ancestry DMP: ", length(mVariants)))

# Genotyped eVariants ----
# Read the genotypes for the sVariants which are also independent cis-eQtls 
# Job run in cluster
#vcftools --gzvcf ${vcf_file} --snps <(zcat /gpfs/projects/bsc83/Data/GTEx_v8/cisQTLs/GTEx_Analysis_v8_sQTL_independent/${tissue_id}.v8.independent_sqtls.txt.gz  | cut -f7 | awk '{if(NR>1)print}')  --recode --stdout | gzip -c > ${pathOut}/${tissue}.isQTLs.Donor_genotypes.vcf.gz
#vcf-query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%SAMPLE=%GTR]\n' ${pathOut}/${tissue}.isQTLs.Donor_genotypes.vcf.gz | gzip -c > ${pathOut}/${tissue}.isQTLs.Donor_genotypes.txt.gz
gt_path <- "/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/Fst/imQTL_GT/"
variants_genotype <- as.data.frame(data.table::fread(paste0(gt_path, tissue, ".imQTLs.Donor_genotypes.txt.gz")))  

# Set proper colnames --
library(stringr)
colnames.variants <- c('CHROM','POS','variant_id','REF','ALT')
colnames.individuals <- str_split_fixed(as.character(variants_genotype[1,-c(1:5)]),'=',2)[,1]
colnames(variants_genotype) <- c(colnames.variants,colnames.individuals)

# Retain selected tissue samples --
variants_genotype <- variants_genotype[,c(colnames.variants, rownames(metadata_2))]
# Donors in same order as in metadata
if(!identical(rownames(metadata_2), colnames(variants_genotype)[6:ncol(variants_genotype)])){
  print("Donors not in the same order in metadata and genotype tables")
  quit()
}

if(length(mVariants) > 1){
  # Subet genotype file and only keep eVariants in eGenes differentially expressed --
  variants_genotype <- variants_genotype[variants_genotype$variant_id %in% mVariants,] 
}

# Customize the genotypes files --
vg.tissue_mod <- cbind(variants_genotype[1:5], apply(variants_genotype[6:ncol(variants_genotype)],2, FUN=function(x){str_split_fixed(x,'=',2)[,2]}))
vg.tissue_wide <- cbind(vg.tissue_mod[1:5], apply(vg.tissue_mod[6:ncol(vg.tissue_mod)],2, FUN=function(x){ifelse(x=='0/0','0',
                                                                                                                 ifelse(x=='0/1','1',
                                                                                                                        ifelse(x=='1/1','2',NA)))}))
vg.tissue_long <- reshape2::melt(vg.tissue_wide,
                                 id.vars = colnames.variants,
                                 variable.name = 'Individual_ID',
                                 value.name = 'genotype')
if(length(mVariants) == 1){
  # Subet genotype file and only keep eVariants in mCpGs differentially methylated --
  variants_genotype <- variants_genotype[variants_genotype$variant_id %in% mVariants,]
  vg.tissue_mod <- vg.tissue_mod[vg.tissue_mod$variant_id %in% mVariants,]
  vg.tissue_wide <- vg.tissue_wide[vg.tissue_wide$variant_id %in% mVariants,]
  vg.tissue_long <- vg.tissue_long[vg.tissue_long$variant_id %in% mVariants,]
}

# Subset tissue samples ----
vg.tissue_long <- vg.tissue_long[vg.tissue_long$Individual_ID %in% rownames(metadata_2),]

# eVariants in VCF -> eVariants in eGenes with at leats 1 independent eQTL with MAF > 0.001 ----
mVariants_lost <- mVariants[!mVariants %in% unique(vg.tissue_long$variant_id)] # MAF < 0.001?
mVariants <- mVariants[mVariants %in% unique(vg.tissue_long$variant_id)]
print(paste0("No. of eVariants (i-mQTL) with MAF > 0.01 in ancestry DMP: ", length(mVariants)))
print(paste0("Minimum MAF value of genotyped eVariants: ",  min(mgene_data[mgene_data$V2 %in% mVariants, "V6"])))

# Add info to summary data ----
d[["imQTL:MAF_b_01"]] <- length(mVariants_lost)
d[["imQTL:MAF_01"]] <- length(mVariants)
d[["minMAF"]] <-  min(mgene_data[mgene_data$V2 %in% mVariants, "V6"])

# Subset eQTLs in eGenes differentially spliced with at least 1 genotyped independent eQTL with MAF > 0.001 ----
mgene_data <- mgene_data[mgene_data$V2 %in% mVariants,]
mgene_data_extra <- mgene_data_extra[mgene_data_extra$variant_id %in% mVariants,]
imGenes.de.maf01 <- c(unique(gsub(':.*','',mgene_data$V1[mgene_data$V7<0.05])),mgene_data_extra$cpg_id[mgene_data_extra$qval<0.05])

# Subset genes with ieQTL with MAF 001 (in vcf file) ----
dea_res <- dea_res[rownames(dea_res) %in% imGenes.de.maf01,]
d[["Ancestry:DMP-imQTL:MAF01"]] <- length(imGenes.de.maf01)


# Create 'dictionary' (named list)  ----
gene_variants.list <- sapply(imGenes.de.maf01, function(x){
  if (x %in% unique(gsub(':.*','',mgene_data$V1[mgene_data$V7<0.05]))) {
    variant <- unique(mgene_data[gsub(':.*','',mgene_data$V1)==x,]$V2)
  } else {
  variant <- unique(mgene_data_extra[mgene_data_extra$cpg_id==x,]$variant_id)
  }
  return(variant)
}, simplify=F)

print('Finished parsing mQTLs')

##### run model #######
############### Linear model for each event #########################

## To build the model:
# 1. For each gene, recover its eVariants (ieQTLS with MAF 001) ----
# 2. Compare 2 models:
# residuals ~ Age + Sex + BMI + cis-effects (ieQTL(s))
# residuals ~ Age + Sex + BMI + cis-effects (ieQTL(s)) + Ancestry 
# Is Ancestry adding someting that the sQTLs did not capture, a not cis-driven effect?
library(limma)

genes <- rownames(dea_res)

# Recover the expression residuals associated with DEGs ----
meth_residuals <- readRDS(paste0(project_path,'Tissues/', tissue,"/methylation_residuals.continous.rds"))

# Susbet residuals of genes to be modelled ----
if(length(rownames(dea_res))==1){
  meth_residuals <- as.data.frame(t(as.matrix(meth_residuals[rownames(dea_res),])))
  rownames(meth_residuals) <- rownames(dea_res)
}else{
  meth_residuals <- meth_residuals[rownames(dea_res),]    
}

#identical(metadata$Sample, colnames(exprs_residuals))
if(!identical(rownames(metadata_2), colnames(meth_residuals))){
  print("The samples in the metadata and the residual data do not coincide")
  q()
}
if(!identical(rownames(metadata_2), as.character(unique(vg.tissue_long$Individual_ID)))){
  print("The donors in the metadata and the genotyped data do not coincide")
  q()
}

# Is the Ancestry effect on the mCpG cis-driven or not cis-driven ----
# Is the mCpG cis-driven or not cis-driven ----
print('*****************************************************')
print(paste0(tissue, ' has a total of ', length(genes), ' DMPs with cis-independent mQTLs with MAF > 0.01'))
print('Applying glm() function to each of those events')
print('*****************************************************')
cat('\n')
cat('\n')

###### code for cis-driven expressn
#### Function to fit  linear model per gene, w and w/ independent cis-eQTL of the eGene ####
source("/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Scripts/Hier_part_mod.R")
library(tidyverse)
library(dplyr)
library(caret)
lm.cis_models <- function(g){
  
  print(paste0("----  ", g, "  ----"))
  
  # Get gene residual expression across samples ----
  if(is.vector(meth_residuals)){
    res.g <- as.data.frame(meth_residuals)
  }else{
    res.g <- as.data.frame(t(meth_residuals[g,])) 
  }
  #res.g_long <- data.frame(Individual_ID=colnames(res.g), residuals=res.g[1,], row.names=NULL) #residuals wide to long
  if (ncol(res.g)>1) {
    res.g_long <- reshape2::melt(res.g, variable.name="Individual_ID", value.name = "residuals")
  } else {
    res.g_long <- res.g
    res.g_long$Individual_ID <- rownames(res.g_long)
    colnames(res.g_long) <- c('residuals','Individual_ID')
  }
  
  # Merge residuals with metadata -> res.md
  metadata_2$Individual_ID <- rownames(metadata_2)
  res.md <- merge(res.g_long, metadata_2, by = 'Individual_ID')
  #res.md.filtered <- res.md[,-which(colnames(res.md)=='Individual_ID')]
  df <- res.md
  df <- df %>% mutate_if(is.character,as.factor)
  print(head(df))
  
  print(paste0('Fitting the model with cis-eQTLs'))
  # 1. Select gene imQTL
  # 2. Select donors with genotyped imQTL
  # 3. Select traits with variance
  # 4. Select imQTLs with variance
  # 5. Run models
  
  # 1. Retrieve the eVariants (ieQTLs in gene) ----
  variants.tissue_g <- gene_variants.list[[g]]
  print(paste0("The CpG ", g, " has a total of ",length(variants.tissue_g)," mVariants"))
  
  # Retrieve the genotype of the eVariants
  vg_g <- droplevels(vg.tissue_long[vg.tissue_long$variant_id%in%variants.tissue_g,])
  print(head(vg_g))
  
  # 2. Filter out snps if any of the individuals is NA ----
  vg_g$Individual_ID <- as.character(vg_g$Individual_ID)
  snps_na <- unique(vg_g[is.na(vg_g$genotype),]$variant_id) # individuals with missing genotypes
  snps_in <- unique(vg_g$variant_id)[!unique(vg_g$variant_id)%in%snps_na] 
  if(length(snps_in)==0){
    print(paste0(g, " has 0 eVariants without NAs"))
    return(NA)
    break
  }
  vg_g.filtered <- vg_g[vg_g$variant_id%in%snps_in,]
  
  # genotypes long to wide (to have each snps as a column to do the lm())
  vg_g.filtered.wide <- reshape2::dcast(vg_g.filtered[,c('variant_id','Individual_ID','genotype')],
                                        Individual_ID ~ variant_id, value.var="genotype")
  
  # Subset res.md (residuals + metadata) to have the same individuals
  res.md.filtered <- droplevels(res.md[res.md$Individual_ID%in%vg_g.filtered$Individual_ID,])
  #res.md.filtered <- res.md.filtered[,-which(colnames(res.md.filtered)=='Sample_ID')]
  
  # Merge res.md.filtered (residuals + metadata) & variants (vg_g.filtered.wide) by Individual_ID
  df <- merge(res.md.filtered, vg_g.filtered.wide, by = 'Individual_ID')
  df <- df %>% mutate_if(is.character,as.factor)
  # All genotyped individuals are EUR for gene ENSG00000235615.2 !!!!
  # if(length(levels(df$Ancestry))==1){
  #   print(paste0(g, " only has genotyped individuals of one ancestry"))
  #   return(NA)
  #   break
  # }
  variants_def <- colnames(df)[!colnames(df)%in%c(c('PEER1','PEER2','PEER3','PEER4','PEER5'),individual_variables,'Individual_ID','residuals')]
  
  
  # 3. check variance of its (its:individual traits) ----
  its_var <- sapply(individual_variables, function(it) var(as.numeric(df[[it]])), simplify=F)
  its_in <- names(its_var)[its_var>0]
  its_in.list <- its_var[names(its_var)%in%its_in]
  
  # 4. check variance of snps before being included in the 'basal' model ----
  #snps_var <- sapply(variants_def, function(v) var(as.numeric(df[[v]])), simplify=F)
  #snps_in <- names(snps_var)[snps_var>0]
  snps_table <- sapply(variants_def, function(v) table(as.numeric(df[[v]])), simplify=F)
  snps_var <- sapply(variants_def, function(v) sum(table(as.numeric(df[[v]]))>2)==length(table(as.numeric(df[[v]]))) & 
                       length(table(as.numeric(df[[v]])))>1, simplify=F)  # at least 3 donors of each genotype
  snps_in <- names(snps_var)[which(snps_var==T)]
  snps_in.list <- snps_var[names(snps_var)%in%snps_in]
  snps.out <- NA
  
  # 5. Run models ----
  if(length(snps_in)>0){
    # Build models and fit lm ----
    modelA <- paste(c(its_in[its_in!="EURv1"],snps_in),collapse = '+')
    modelB <- paste(c(its_in[its_in!="EURv1"],snps_in,"EURv1"),collapse = '+')
    lm_formula.y <- paste0('residuals')
    lm_formula.modelA <- paste(c(lm_formula.y, modelA),collapse='~')
    lm_formula.modelB <- paste(c(lm_formula.y, modelB),collapse='~')
    multi.fit.modelA <- lm(lm_formula.modelA,
                           data=df)
    multi.fit.modelB <- lm(lm_formula.modelB,
                           data=df)
    #Check if there are linearly dependent sVariants ----
    if(length(snps_in)>1){
      # cor matrix ----
      df.snps <- df[,colnames(df)%in%snps_in]
      df.snps %>% mutate_if(is.factor, as.numeric) -> df.snps
      correlationMatrix <- cor(df.snps)
      
      # findCorrelation from {caret} package (filter cutoff=0.9 -> out) ----
      #print('Using findCorrelation from caret package to filter linearly dependent sVariants: cutoff>0.9 -> out ')
      highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.85, names = T)
      snps.out <- ifelse(length(highlyCorrelated)==0, NA, paste(highlyCorrelated,collapse = ":"))
      
      if(!is.na(snps.out)){
        print(paste0(snps.out, " is linearly dependent and thus excluded"))
        # Exclude correlated
        snps.f.out <- unlist(strsplit(snps.out, split = ":"))
        snps_in <- snps_in[!snps_in%in%snps.f.out]    
        if(length(snps_in)==0){
          print(paste0(g, " has 0 eVariants with variance"))
          return(NA)
          break
        }
        # Update models
        modelA <- paste(c(its_in[its_in!="EURv1"],snps_in),collapse = '+')
        modelB <- paste(c(its_in[its_in!="EURv1"],snps_in,"EURv1"),collapse = '+')
        lm_formula.y <- paste0('residuals')
        lm_formula.modelA <- paste(c(lm_formula.y, modelA),collapse='~')
        lm_formula.modelB <- paste(c(lm_formula.y, modelB),collapse='~')
        multi.fit.modelA <- lm(lm_formula.modelA,
                               data=df)
        multi.fit.modelB <- lm(lm_formula.modelB,
                               data=df)
      }
    }else{
      snps.f.out <- c()
    } 
    # 6. g report ----
    # report data frame
    g.df <- data.frame(Tissue = tissue,
                       ensembl.id = g,
                       indv_total = length(unique(metadata_2$Individual_ID)),
                       # indv_in = length(indv_in),
                       # indv_lost = length(indv_na),
                       EUR_total = sum(res.md$EURv1>0.75),
                       AFR_total = sum(res.md$EURv1<0.75),
                       # EUR_in = sum(res.md[res.md$Individual_ID%in%indv_in,]$EURv1>0.75),
                       # AFR_in = sum(res.md[res.md$Individual_ID%in%indv_in,]$EURv1<0.75),
                       # EUR_lost = table(res.md[res.md$Individual_ID%in%indv_na,]$EURv1>0.75),
                       # AFR_lost = table(res.md[res.md$Individual_ID%in%indv_na,]$URv1<0.75),
                       num_snps = length(variants_def),
                       num_snps_in = length(snps_in),
                       num_snps_lost =  length(variants_def)-length(snps_in),
                       num_snps_filtered_out =  ifelse(is.na(snps.out), NA, length(snps.f.out))
    )
    rownames(g.df) <- g
    
    # 7. hier.part ----
    # if(length(snps_in) <= 8){
    #   hier.part.out <- hier.part.mod(y=df[,"residuals"], x=df[,c(its_in[its_in!="EURv1"],snps_in,"EURv1")], fam = "gaussian", gof = "Rsqu")
    #   if(length(hier.part.out)==1){# & is.na(hier.part.out)){
    #     r2.traits <- sapply(its_in, function(i) NA)
    #     r2.eVariants <- NA
    #   }else{
    #     hier.part.results <- hier.part.out$IJ
    #     r2.traits <- sapply(its_in, function(i) hier.part.results[i,1])
    #     r2.eVariants <- sum(sapply(rownames(hier.part.results)[!rownames(hier.part.results) %in% its_in], function(i) hier.part.results[i,1]))    
    #   }
    # }else{
    #   r2.traits <- sapply(its_in, function(i) NA)
    #   r2.eVariants <- NA
    # }
    # names(r2.eVariants) <- "eVariants"
    
    # 8.output ----
    res <- list("report_summary" = g.df,
                "lm.modelA" = multi.fit.modelA,
                "lm.modelB" = multi.fit.modelB,
                "df" = df[,colnames(df) %in% c("Individual_ID", "residuals", its_in, snps_in)]#,
                # "hier.part" = c(r2.traits, r2.eVariants,
                #                 "n.eVariants" = length(snps_in),
                #                 "EA" = table(res.md[res.md$Individual_ID%in%indv_in,]$EURv1>0.75),
                #                 "AA" = table(res.md[res.md$Individual_ID%in%indv_in,]$EURv1<0.75))
    )
    
    return(res)
    
  }else{
    print(paste0(g, " has 0 snps with variance"))
    return(NA)
    break
  }
  
}


# lm per event ----
g.lm_models <- sapply(genes, function(g) lm.cis_models(g), simplify = F)
names(g.lm_models) <-  genes
# genes that  cannot modelled ----
d[["Ancestry:DMP:NotModelled"]] <- sum(is.na(g.lm_models))
# genes modelled ----
g.lm_models <- g.lm_models[!is.na(g.lm_models)] # SNP with no variance or no SNP not correlated
d[["Ancestry:DMP:Modelled"]] <- length(g.lm_models)

# Compare modelA vs modelB ----
deg <- names(g.lm_models)
# Anova
anova.P.Value <- sapply(deg, function(g) anova(g.lm_models[[g]][["lm.modelA"]],g.lm_models[[g]][["lm.modelB"]], test="F")[2,6])
# Multiple testing correction
anova.adj.P.Val <- p.adjust(anova.P.Value, method = "BH")

# Parse results data ----
results.df <- cbind.data.frame(deg, anova.adj.P.Val)
results.df$Class <- ifelse(results.df$anova.adj.P.Val < 0.05, "Not_cis-driven","Cis-driven")

# Parse hier.part results ----
# if(length(unique(sapply(deg, function(g) length(g.lm_models[[g]][["hier.part"]])))) > 1){
#   print("Some genes do not have trait variance")
#   q()
# }
# 
# hier.part.results <- do.call(rbind.data.frame,
#                              lapply(deg, function(g) g.lm_models[[g]][["hier.part"]]))
# #colnames(hier.part.results) <- c(its, "eVariants",  "n.eVariants", "EA","AA")
# rownames(hier.part.results) <- deg
# 
# if(identical(rownames(hier.part.results), rownames(results.df))){
#   results.df <- cbind.data.frame(results.df, hier.part.results)
# }

# Save ancestry-DEG classified ----
saveRDS(results.df,
        paste0(project_path,'Tissues/',tissue,'/Ancestry_DMP.Classified.qval005.rds'))
saveRDS(g.lm_models,
        paste0(project_path,'Tissues/',tissue,'/Ancestry_DMP.g.lm_models.qval005.rds'))


# Report summary ----
report_df <- do.call(rbind.data.frame,
                     lapply(deg, function(gene)
                       g.lm_models[[gene]][["report_summary"]]))
saveRDS(report_df,
        paste0(project_path,'Tissues/',tissue,'/',tissue,'.imQTLs.Report_summary.qval005.rds'))

# d ----
saveRDS(d,
        paste0(project_path,'Tissues/',tissue,'/',tissue,'.Ancestry_DMP.Classification_summary.qval005.rds'))

# eVariant not in VCF
# saveRDS(eVariants_lost,
#         paste0(project_path, 'Tissues/',tissue, ".mVariants_not_in_vcf.rds"))

# ---------------------- #
end_time <- Sys.time()
# ---------------------- #

# ---------------------- #
print("")
print("# ---- Elapsed time ---- #")
end_time - start_time
print("# ---------------------- #\n")
# ---------------------- #

# ---------------------- #
print("")
print("# ---- Session info ---- #")
sessionInfo()
print("# ---------------------- #\n")
print("")
# ---------------------- #

