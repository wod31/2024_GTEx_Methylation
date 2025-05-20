#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Model if DVPs are genetically driven bu mQTLs (cis-driven)
# @software version: R=4.2.2

#Loading libraries
library(optparse)
library(stringr)
library(limma)
library(missMethyl)
suppressMessages(library(tidyverse))
library(dplyr)
library(caret)

# Parsing
parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-t", "--tissue"), type="character",
                     dest="tissue",
                     help="Tissue")
parser <- add_option(parser, opt_str=c("-i", "--input"), type="character",
                     dest="input", default="M", 
                     help="Type of input, either M values or residuals")
parser <- add_option(parser, opt_str=c("-m", "--model"), type="character",
                     dest="model", default="DVP", 
                     help="What do we want to model? Either DMP or DVP")
options=parse_args(parser)
tissue=options$tissue
input=options$input
model=options$model

# tissue <- "WholeBlood"

# input <- "M"
# input <- "residuals"

# model <- "DMP"
# model <- "DVP"
# correlation <- 0.4
correlation <- 0.8
# correlation <- 0.6

distance <- 250000
# distance <- 100000

print(tissue)
print(correlation)
print(distance)
# first_dir <- "~/Documents/mn4/"
first_dir <- "/gpfs/"
project_path <- paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Methylation/")

start_time <- Sys.time()

print("Reading metadata")
metadata <- readRDS(paste0(project_path, 'Tissues/', tissue, "/metadata.rds"))

if(input=="M"){
  print("Reading data")
  beta <- readRDS(paste0(project_path,'Tissues/', tissue, "/data.rds"))

  print("Computing M values")
  beta <- beta[,metadata$SUBJID]
  probes <- rownames(beta)
  beta <- sapply(beta, as.numeric)
  rownames(beta) <- probes
  M <- log2(beta/(1-beta))
}else if(input=="residuals"){ #We will be calling M the variable for both M values or residuals
  M <- readRDS(paste0(project_path,'Tissues/', tissue,"/methylation_residuals.continous.rds"))
}

print("Reading Admixture results")
admixture_ancestry <- read.table(paste0(project_path, 'Data/admixture_inferred_ancestry.txt'))
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

covariates <- c("TRISCHD", "DTHHRDY", "PEER1", "PEER2", "PEER3", "PEER4", "PEER5")
metadata_2 <- metadata[,c(covariates, individual_variables)]

# The genotypes of the variants (mQTLs) are modelled as a factor, where: 0|0 -> 0, 0|1 -> 1; 1|1 ->2

print("metadata is prepared")

print("Reading mQTL data")
# mCpGs and independent mQTLs ----
inpath_mqtls <- "/gpfs/scratch/bsc83/MN4/bsc83/bsc83535/GTEx/v9/mQTLs/"
mgene_data <- as.data.frame(data.table::fread(paste0(inpath_mqtls,tissue,".mQTLs.conditional.txt.gz")))
mgene_data <- mgene_data[mgene_data$V7<0.05 & abs(mgene_data$V3)<distance,] #Distance from mQTL to CpG
mCpGs <- unique(gsub(':.*','',mgene_data$V1[mgene_data$V7<0.05 & abs(mgene_data$V3)<distance]))

mgene_data_extra <- as.data.frame(data.table::fread(paste0(inpath_mqtls,tissue,".regular.perm.fdr.txt")))
mgene_data_extra <- mgene_data_extra[!mgene_data_extra$cpg_id %in% mCpGs,]
mgene_data_extra <- mgene_data_extra[mgene_data_extra$qval<0.05,]
mCpGs <- c(mCpGs,mgene_data_extra$cpg_id)

#### DM Analysis
if(model=="DVP"){
  dva_res <- readRDS(paste0(project_path,'Tissues/', tissue,"/DVP_Ancestry.rds"))
  colnames(dva_res)[6] <- "adj.P.Val" #Here we only had significant CpGs already
  dva_res <- dva_res[dva_res$adj.P.Val < 0.05,]
  dva_res <- dva_res[rownames(dva_res) != 'cg22122606',]
} else if(model=="DMP"){
  dva_res <- readRDS(paste0(project_path,'Tissues/', tissue,"/DML_results_5_PEERs_continous.rds"))[["EURv1"]]
  dva_res <- dva_res[dva_res$adj.P.Val < 0.05,] # significant CpGs
}

# Save summary data ----
d <- list()
d[["mQTLs:CpG"]] <- length(mCpGs)
d[["Ancestry:DVP"]] <- sum(dva_res$adj.P.Val < 0.05)
d[["Ancestry:DVP:eGpG"]] <- sum(rownames(dva_res[dva_res$adj.P.Val < 0.05,])  %in% mCpGs)

# Print to std.out summary data ----
print(paste0("No. of mQTLs:CpG: ", length(mCpGs)))
print(paste0("No. of Ancestry:DVP:  ", sum(dva_res$adj.P.Val < 0.05)))
print(paste0("Ancestry:DVP:eGpG:  ", sum(rownames(dva_res[dva_res$adj.P.Val < 0.05,])  %in% mCpGs)))

# Select differentially methylated cpgs with at least 1 independent mQTL ----
dva_res <- dva_res[rownames(dva_res) %in% mCpGs,] # i-mCpG

# Subset mQTLs in mCpGs differentially methylated with at least 1 independent mQTL ----
mgene_data <- mgene_data[gsub(':.*','',mgene_data$V1) %in% rownames(dva_res),]

# eVariants in eCpGs differentially methylated with at leats 1 independent mQTL ----
mVariants <- unique(mgene_data$V2) # imQTLs
print(paste0("No. of mVariants (i-mQTL) in ancestry DMP: ", length(mVariants)))

# Genotyped Variants: Filtered in a previous step in the cluster that looked something like this:
#vcftools --gzvcf ${vcf_file} --snps <(zcat /gpfs/projects/bsc83/Data/GTEx_v8/cisQTLs/GTEx_Analysis_v8_sQTL_independent/${tissue_id}.v8.independent_sqtls.txt.gz  | cut -f7 | awk '{if(NR>1)print}')  --recode --stdout | gzip -c > ${pathOut}/${tissue}.isQTLs.Donor_genotypes.vcf.gz
#vcf-query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%SAMPLE=%GTR]\n' ${pathOut}/${tissue}.isQTLs.Donor_genotypes.vcf.gz | gzip -c > ${pathOut}/${tissue}.isQTLs.Donor_genotypes.txt.gz
gt_path <- paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Methylation/Data/Fst/imQTL_GT/")
variants_genotype <- as.data.frame(data.table::fread(paste0(gt_path, tissue, ".imQTLs.Donor_genotypes.txt.gz")))  

# Set proper colnames --
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
  # Subet genotype file and only keep mVariants in mCpGs differentially methylated --
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
vg.tissue_long <- vg.tissue_long[vg.tissue_long$Individual_ID %in% rownames(metadata_2),] #This is already fullfiled, just double checking

# mVariants in VCF -> mVariants in mCpGs with at leats 1 independent mQTL with MAF > 0.01 ----
#In this case all are already MAF > 0.01, we are just double checking
mVariants_lost <- mVariants[!mVariants %in% unique(vg.tissue_long$variant_id)] # MAF < 0.01?
mVariants <- mVariants[mVariants %in% unique(vg.tissue_long$variant_id)]
print(paste0("No. of eVariants (i-mQTL) with MAF > 0.01 in ancestry DMP: ", length(mVariants)))
print(paste0("Minimum MAF value of genotyped eVariants: ",  min(mgene_data[mgene_data$V2 %in% mVariants, "V6"])))

# Add info to summary data ----
d[["imQTL:MAF_b_01"]] <- length(mVariants_lost)
d[["imQTL:MAF_01"]] <- length(mVariants)
d[["minMAF"]] <-  min(mgene_data[mgene_data$V2 %in% mVariants, "V6"])

# Subset mQTLs in mCpGs differentially methylated with at least 1 genotyped independent mQTL with MAF > 0.01 ----
mgene_data <- mgene_data[mgene_data$V2 %in% mVariants,] #double checking
imCpGs.de.maf01 <- unique(gsub(':.*','',mgene_data$V1[mgene_data$V7<0.05]))

# Subset CpGs with imQTL with MAF 001 (in vcf file) ----
dva_res <- dva_res[rownames(dva_res) %in% imCpGs.de.maf01,] #just double checking as V6 always was the MAF we want
d[["Ancestry:DVP-imQTL:MAF001"]] <- length(imCpGs.de.maf01)


# Create 'dictionary' (named list)  ----
gene_variants.list <- sapply(imCpGs.de.maf01, function(x){
  variant <- unique(mgene_data[gsub(':.*','',mgene_data$V1)==x,]$V2)
  return(variant)
}, simplify=F)

print('Finished parsing mQTLs')

##### run model #######
############### Linear model for each event #########################

## To build the model:
# Compare 2 models:
# M_values ~ batch_effects + Age + Sex + BMI + cis-effects (imQTL(s)) 
# M_values ~ batch_effects + Age + Sex + BMI + cis-effects (imQTL(s)) + Ancestry 
# Is Ancestry adding someting that the mQTLs did not capture, a not cis-driven effect?

CpGs <- rownames(dva_res) #the naming is maintained from older scripts, but here CpGs would refer to CpGs

# Susbet CpGs to be modelled ----
if(length(rownames(dva_res))==1){
  meth_M <- as.data.frame(t(as.matrix(M[rownames(dva_res),])))
  rownames(meth_M) <- rownames(dva_res)
}else{
  meth_M <- M[rownames(dva_res),]    
}
if(!identical(rownames(metadata_2), colnames(meth_M))){
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
print(paste0(tissue, ' has a total of ', length(CpGs), ' DVPs with cis-independent mQTLs with MAF > 0.01'))
print('*****************************************************')
cat('\n')


###### code for cis-driven 
save.image(paste0(project_path,'Tissues/',tissue,'/Testing_image_', model, '_', correlation, '.RData'))
# load(paste0(project_path,'Tissues/',tissue,'/Testing_image.RData'))

# g <- CpGs[1]
lm.cis_models <- function(g){
  # g <- "cg11028037" #Why variant 25 gets an error for this cpg in blood?
  print(paste0("----  ", g, "  ----"))
  
  # Edit format ----
  if(is.vector(meth_M)){
    res.g <- as.data.frame(meth_M)
  }else{
    res.g <- as.data.frame(t(meth_M[g,])) 
  }
  
  if (ncol(res.g)>1) {
    res.g_long <- reshape2::melt(res.g, variable.name="Individual_ID", value.name = "residuals", id.vars = NULL)
  } else {
    res.g_long <- res.g
    res.g_long$Individual_ID <- rownames(res.g_long)
    colnames(res.g_long) <- c('residuals','Individual_ID') #Here the naming residuals refer to M values
  }
  
  # Merge M-values with metadata -> res.md
  metadata_2$Individual_ID <- rownames(metadata_2)
  res.md <- merge(res.g_long, metadata_2, by = 'Individual_ID')
  df <- res.md
  df <- df %>% mutate_if(is.character,as.factor)
  
  print(paste0('Fitting the model with cis-mQTLs'))
  # 1. Select gene imQTL
  # 2. Select donors with genotyped imQTL
  # 3. Select traits with variance
  # 4. Select imQTLs with variance
  # 5. Run models
  
  # 1. Retrieve the mVariants (imQTLs in CpG) ----
  variants.tissue_g <- gene_variants.list[[g]]
  print(paste0("The CpG ", g, " has a total of ",length(variants.tissue_g)," mVariants"))
  
  # Retrieve the genotype of the eVariants
  vg_g <- droplevels(vg.tissue_long[vg.tissue_long$variant_id%in%variants.tissue_g,])
  
  # 2. Filter out snps if any of the individuals is NA ----
  vg_g$Individual_ID <- as.character(vg_g$Individual_ID)
  snps_na <- unique(vg_g[is.na(vg_g$genotype),]$variant_id) # individuals with missing genotypes
  snps_in <- unique(vg_g$variant_id)[!unique(vg_g$variant_id)%in%snps_na] 
  vg_g.filtered <- vg_g[vg_g$variant_id%in%snps_in,]
  
  if(length(snps_in)==0){
    print(paste0(g, " has 0 eVariants with all genotypes"))
    return(NA)
    break
  }
  
  # genotypes long to wide (to have each snps as a column to do the lm())
  vg_g.filtered.wide <- reshape2::dcast(vg_g.filtered[,c('variant_id','Individual_ID','genotype')],
                                        Individual_ID ~ variant_id, value.var="genotype")
  
  # Subset res.md (residuals + metadata) to have the same individuals
  res.md.filtered <- droplevels(res.md[res.md$Individual_ID%in%vg_g.filtered$Individual_ID,])

  # Merge res.md.filtered (residuals + metadata) & variants (vg_g.filtered.wide) by Individual_ID
  df <- merge(res.md.filtered, vg_g.filtered.wide, by = 'Individual_ID')
  df <- df %>% mutate_if(is.character,as.factor)

  variants_def <- colnames(df)[!colnames(df)%in%c(covariates, individual_variables,'Individual_ID','residuals')]
  
  
  # 3. check variance of individual traits ---- Relevant for later on
  its_var <- sapply(c(covariates, individual_variables), function(it) var(as.numeric(df[[it]])), simplify=F)
  its_in <- names(its_var)[its_var>0]
  its_in.list <- its_var[names(its_var)%in%its_in]
  if(input=="residuals"){
    its_in  <- its_in[!its_in %in% covariates]
    its_in.list <- its_var[names(its_var)%in%its_in]
  }
  # 4. check variance of snps before being included in the 'basal' model ----
  snps_table <- sapply(variants_def, function(v) table(as.numeric(df[[v]])), simplify=F)
  snps_var <- sapply(variants_def, function(v) sum(table(as.numeric(df[[v]]))>2)==length(table(as.numeric(df[[v]]))) & 
                       length(table(as.numeric(df[[v]])))>1, simplify=F)  # at least 3 donors of each genotype
  snps_in <- names(snps_var)[which(snps_var==T)]
  snps_in.list <- snps_var[names(snps_var)%in%snps_in]
  snps.out <- NA
  
  # 5. Run models ----
  if(length(snps_in)>0){
    #Check if there are linearly dependent sVariants ----
    if(length(snps_in)>1){
      # cor matrix ----
      df.snps <- df[,colnames(df)%in%snps_in]
      df.snps %>% mutate_if(is.factor, as.numeric) -> df.snps
      correlationMatrix <- cor(df.snps)
      
      # findCorrelation from {caret} package (filter cutoff=correlation -> out) ----
      # highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.8, names = T)
      highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=correlation, names = T)
      snps.out <- ifelse(length(highlyCorrelated)==0, NA, paste(highlyCorrelated,collapse = ":"))
      
      if(!is.na(snps.out)){
        # Exclude correlated
        snps.f.out <- unlist(strsplit(snps.out, split = ":"))
        snps_in <- snps_in[!snps_in%in%snps.f.out]  
        print(paste0("Some SNPs are linearly dependent and thus excluded, we will keep: ", length(snps_in)))
        if(length(snps_in)==0){
          print(paste0(g, " has 0 eVariants with variance"))
          return(NA)
          break
        }
      }
    }else{
      snps.f.out <- c()
    } 
    # Build models and fit lm ----
    modelA <- paste(c(its_in[its_in!="EURv1"],snps_in),collapse = '+')
    modelB <- paste(c(its_in[its_in!="EURv1"],snps_in,"EURv1"),collapse = '+')
    lm_formula.y <- paste0('residuals') #Here residuals is the name of M values
    lm_formula.modelA <- paste(c(lm_formula.y, modelA),collapse='~')
    lm_formula.modelB <- paste(c(lm_formula.y, modelB),collapse='~')
    if(model=="DMP"){
      multi.fit.modelA <- lm(lm_formula.modelA, data=df)
      multi.fit.modelB <- lm(lm_formula.modelB, data=df)
    }else if(model=="DVP"){
      mod_A <- model.matrix( as.formula(lm_formula.modelA), data =  df)
      mod_B <- model.matrix( as.formula(lm_formula.modelB), data =  df)
      # df$residuals <-  M[g,] #Just to double check
      if(ncol(mod_B)>nrow(mod_B)){
        print("Number of variables is higher than the number of samples, which means that we do not have enough power to test")
        return(NA)
        next
      } 
      
      #Running an edited version of varFit from getAnywhere(varFit.default)
      data <- as.matrix(t(df$residuals))
      mod_A <- as.matrix(mod_A)
      coef <- c(1:ncol(mod_A))
      # z <- getLeveneResiduals(data, design = mod_A[, coef], #This functions outputs NAs in LogVarRatio because we are including more than one coefficient
      #                         type = NULL)
      #getLeveneResiduals performs matrix multiplications after qr()
      #It some column is linearly dependent to others, rank is smaller than the number of columns and we cannot do the matrix multiplication
      #We already filter variants that are correlated in a previous step, but the variants contain 1,2 and 3.
      #If a variant X has the same 3s as variant Y but different 1 and 2s, we will not have filtered out in the previous step, as correlation is not higher than 0.8 (or whatever threshold we use)
      #But when we create the desing matrix, we realize that two covariates, the ones that model 3s for X and 3s for Y are equal and then we cannot model
      #Instead of running the function getLeveneResiduals I will use an adapted code from getAnywhere(getLeveneResiduals) and remove the instances from the desing where this happen

      #Start of code adapted from getLeveneResiduals
      data <- as.matrix(data)
      y <- rowMeans(data)
      type = "AD"

      design <- mod_A[, coef]
      QR <- qr(design)
      if(ncol(QR$qr)>QR$rank){
        print(paste(ncol(QR$qr)-QR$rank, "design columns have been removed from the model due to linearly dependent variants")) #Correlated at the allele level, even though the two variants are different, one allele of each correlate 100%
        design <- design[,QR$pivot[seq_len(QR$rank)]] #Removing columns linearly depending to others
        QR <- qr(design)
      }
      V.inv <- chol2inv(QR$qr, size = QR$rank)
      des.y <- data %*% design
      beta.hat <- des.y %*% t(V.inv) #des.y <- 1 35, t(V.inv) 34 34. One variant is not popping up in V.inv (25 for the example I am checking)
        
      ybar <- beta.hat %*% t(design)
      n.inv <- diag(design %*% V.inv %*% t(design))
      n.inv[n.inv==1] <- 0.9999 #I added this step as well
      n <- 1/n.inv
      temporal <- n/(n - 1)
      temporal[temporal<0] <- 0 #I added this step as well
      lvg <- sqrt(temporal)
      lvg <- matrix(rep(lvg, nrow(data)), byrow = TRUE, nrow = nrow(data))
      z <- abs(data - ybar)
      z <- z * lvg
      #End of code adapted from getLeveneResiduals
    
      df_final <- df
      df_final$residuals <- as.vector(z) # we will model now z, but in very few cases a value for a donor is Inf. Find out why and fix it
      multi.fit.modelA <- lm(lm_formula.modelA, data= df_final) #Equivalent to limma::lmFit(z$data, mod_A, weights = NULL) from varFit(t(df$residuals), design = mod_A)
      
      mod_B <- as.matrix(mod_B)
      coef <- c(1:ncol(mod_B))

      #I copied the same chunck of code from before, in the future I will create it as my_levene_function
      design <- mod_B[, coef]
      QR <- qr(design)
      if(ncol(QR$qr)>QR$rank){
        print(paste(ncol(QR$qr)-QR$rank, "design columns have been removed from the model due to linearly depent variants")) #Correlated at the allele level, even though the two variants are different, one allele of each correlate 100%
        design <- design[,QR$pivot[seq_len(QR$rank)]] #Removing columns linearly depending to others
        QR <- qr(design)
      }
      V.inv <- chol2inv(QR$qr, size = QR$rank)
      des.y <- data %*% design
      beta.hat <- des.y %*% t(V.inv) #des.y <- 1 35, t(V.inv) 34 34. One variant is not popping up in V.inv (25 for the example I am checking)
      
      ybar <- beta.hat %*% t(design)
      n.inv <- diag(design %*% V.inv %*% t(design))
      n.inv[n.inv==1] <- 0.9999 #I added this step as well
      n <- 1/n.inv
      temporal <- n/(n - 1)
      temporal[temporal<0] <- 0 #I added this step as well
      lvg <- sqrt(temporal)
      lvg <- matrix(rep(lvg, nrow(data)), byrow = TRUE, nrow = nrow(data))
      z <- abs(data - ybar)
      z <- z * lvg
      #End of adapted code
      
      df_final <- df
      df_final$residuals <- as.vector(z)
      multi.fit.modelB <- lm(lm_formula.modelB, data= df_final) 
    }

    
    # 6. g report ----
    # report data frame
    g.df <- data.frame(Tissue = tissue,
                       ensembl.id = g,
                       indv_total = length(unique(metadata_2$Individual_ID)),
                       EUR_total = sum(res.md$EURv1>0.75),
                       AFR_total = sum(res.md$EURv1<0.75),
                       num_snps = length(variants_def),
                       num_snps_in = length(snps_in),
                       num_snps_lost =  length(variants_def)-length(snps_in),
                       num_snps_filtered_out =  ifelse(is.na(snps.out), NA, length(snps.f.out))
    )
    rownames(g.df) <- g
    
    # 8.output ----
    res <- list("report_summary" = g.df,
                "lm.modelA" = multi.fit.modelA,
                "lm.modelB" = multi.fit.modelB,
                "df" = df[,colnames(df) %in% c("Individual_ID", "residuals", its_in, snps_in)]
    )
    
    return(res)
    
  }else{
    print(paste0(g, " has 0 snps with variance"))
    return(NA)
    break
  }
}


# lm per event ----
g.lm_models <- sapply(CpGs, function(g) lm.cis_models(g), simplify = F)
names(g.lm_models) <-  CpGs
# CpGs that cannot be modelled ----
d[["Ancestry:DVP:NotModelled"]] <- sum(is.na(g.lm_models))
# CpGs modelled ----
g.lm_models <- g.lm_models[!is.na(g.lm_models)] # SNP with no variance or no SNP not correlated
d[["Ancestry:DVP:Modelled"]] <- length(g.lm_models)

# Compare modelA vs modelB ----
dvp <- names(g.lm_models)

#Testing:
# for(g in dvp){
#   print(anova(g.lm_models[[g]][["lm.modelA"]],g.lm_models[[g]][["lm.modelB"]], test="F")[2,])
# }
#I cannot do anova on some CpGs because of colinearity, but if I change threshold on correlations I can do it but for some others I used to do it I can no longer do it
# g <- "cg14674796"
# model_t <- lm.cis_models(g) #0.6 gets 22 variants and an NA on anova. With 0.8 30 variant I can test
# anova(model_t[["lm.modelA"]], model_t[["lm.modelB"]], test="F")[2,]

anova.P.Value <- sapply(dvp, function(g) anova(g.lm_models[[g]][["lm.modelA"]], g.lm_models[[g]][["lm.modelB"]], test="F")[2,6])
# Multiple testing correction
anova.adj.P.Val <- p.adjust(anova.P.Value, method = "BH")

# Parse results data ----
results.df <- cbind.data.frame(dvp, anova.adj.P.Val)
results.df$Class <- ifelse(results.df$anova.adj.P.Val < 0.05, "Not_cis-driven","Cis-driven")


# Save ancestry-DEG classified ----
saveRDS(results.df,
        # paste0(project_path,'Tissues/',tissue,'/Ancestry_DVP.Classified.rds'))
        paste0(project_path,'Tissues/',tissue,'/Testing_', model, '_', input, "_", correlation, '_Ancestry.Classified.qval005.rds'))
saveRDS(g.lm_models,
        # paste0(project_path,'Tissues/',tissue,'/Ancestry_DVP.g.lm_models.rds'))
        paste0(project_path,'Tissues/',tissue,'/Testing_', model, '_', input, "_", correlation, '_Ancestry.g.lm_models.qval005.rds'))


# Report summary ----
report_df <- do.call(rbind.data.frame,
                     lapply(dvp, function(gene)
                       g.lm_models[[gene]][["report_summary"]]))
saveRDS(report_df,
        # paste0(project_path,'Tissues/',tissue,'/DVP.imQTLs.Report_summary.rds'))
        paste0(project_path,'Tissues/',tissue,'/Testing_', model, '_', input, "_", correlation, '_imQTLs.Report_summary.qval005.rds'))

# d ----
# saveRDS(d, paste0(project_path,'Tissues/',tissue,'/Ancestry_DVP.Classification_summary.rds'))
saveRDS(d, paste0(project_path,'Tissues/',tissue,'/Testing_', model, '_', input, "_", correlation, '_Ancestry.Classification_summary.qval005.rds'))

# ---------------------- #
end_time <- Sys.time()
# ---------------------- #

# ---------------------- #
print("")
print("# ---- Elapsed time ---- #")
end_time - start_time
print("# ---------------------- #\n")
# ---------------------- #
