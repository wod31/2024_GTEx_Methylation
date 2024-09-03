# Tissues ---
first_dir <- "/gpfs/projects/bsc83/"

project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Age", "Ancestry", "BMI", "Sex")
#data <- matrix(nrow = length(tissues), ncol=4, dimnames = list(tissues, names))

tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Jose/00_Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

# Tissue metadata ----
metadata <- lapply(tissues, function(tissue) readRDS(paste0(project_path,"Tissues/", tissue, "/metadata.rds")))
names(metadata) <- tissues

# ---- Analysis ---- ####
# 1.1 mGene data ----
inpath_mqtls <- "/gpfs/scratch/bsc83/MN4/bsc83/bsc83535/GTEx/v9/mQTLs/"
mCpGs <- lapply(tissues, function(tis) {
  mgene_data <- as.data.frame(data.table::fread(paste0(inpath_mqtls,tis,".mQTLs.conditional.txt.gz")))
  mgene_data <- mgene_data[mgene_data$V7<0.05 & abs(mgene_data$V3)<250000,]
  mCpGs <- unique(gsub(':.*','',mgene_data$V1[mgene_data$V7<0.05 & abs(mgene_data$V3)<250000]))
  
  mgene_data_extra <- as.data.frame(data.table::fread(paste0(inpath_mqtls,tis,".regular.perm.fdr.txt")))
  mgene_data_extra <- mgene_data_extra[!mgene_data_extra$cpg_id %in% mCpGs,]
  mgene_data_extra <- mgene_data_extra[mgene_data_extra$qval<0.05,]
  mCpGs <- c(mCpGs,mgene_data_extra$cpg_id)
  mCpGs
  
})
names(mCpGs) <- tissues


# 1.2 Differential methylation results ----
results_DML <- lapply(tissues, function(tis) 
  readRDS(paste0(project_path,"Tissues/",tis,"/DML_results_5_PEERs_continous.rds")))
names(results_DML) <- tissues

# DMP ---
results_DML_all <- lapply(tissues, function(tis) 
  results_DML[[tis]][['EURv1']])
names(results_DML_all) <- tissues

## sig results 
DMP <- lapply(tissues, function(tis) 
  rownames(results_DML_all[[tis]][results_DML_all[[tis]]$adj.P.Val<0.05,]))
names(DMP) <- tissues

### DVPs ####
results_Dvp <- lapply(tissues, function(tis) 
  readRDS(paste0(project_path,"Tissues/",tis,"/DVP_Ancestry.rds")))
names(results_Dvp) <- tissues

res_hypergeometric_dmps <- list()
for(tissue in tissues) {
  sig <- DMP[[tissue]]
  all <- rownames(results_DML_all[[tissue]])

  mCpGs_t <- mCpGs[[tissue]]
  print(head(mCpGs_t))
  common <- sig[sig%in% mCpGs_t]

  # result <- phyper(q = length(common) - 1,
  #                  m = sum(all %in% mCpGs_t),
  #                  n = length(all) - sum(all %in% mCpGs_t),
  #                  k = length(sig),
  #                  lower.tail = FALSE)
  # res_hypergeometric_dmps[[tissue]] <- result

  ### test significance
  m <- matrix(c(length(common), length(sig[!sig%in% mCpGs_t]), length(mCpGs_t[!mCpGs_t %in% sig]), length(all[!all %in% c(sig,mCpGs_t)])), 2,2, byrow = T)
  print(m)

  m[is.na(m)] <- 0
  #m <- m[c(type,paste0('No ',type)),]
  print(m)
  f <- fisher.test(m)
  print(f)
  res_hypergeometric_dmps[[tissue]] <- (list("f" = f, "m" = length(common)))
}
saveRDS(res_hypergeometric_dmps, paste0(project_path,"Tissues/enrichment_mcpgs_dmps.fisher.qval005.rds"))


# res_hypergeometric_dvps <- list()
# for(tissue in tissues) {
#   sig <- results_Dvp[[tissue]]
#   all <- rownames(results_DML_all[[tissue]])
#   
#   mCpGs_t <- mCpGs[[tissue]]
#   common <- sig[rownames(sig) %in% mCpGs_t,]
#   # 
#   # result <- phyper(q = nrow(common) - 1,
#   #                  m = sum(all %in% mCpGs_t),
#   #                  n = length(all) - sum(all %in% mCpGs_t),
#   #                  k = nrow(sig),
#   #                  lower.tail = FALSE)
#   # res_hypergeometric_dvps[[tissue]] <- result
#   
#   m <- matrix(c(nrow(common), nrow(sig[!rownames(sig) %in% mCpGs_t,]), 
#                 length(mCpGs_t[!mCpGs_t %in% rownames(sig)]), length(all[!all %in% c(rownames(sig),mCpGs_t)])), 2,2, byrow = T)
#   print(m)  
#   
#   m[is.na(m)] <- 0
#   #m <- m[c(type,paste0('No ',type)),]
#   print(m)
#   f <- fisher.test(m)
#   print(f)
#   res_hypergeometric_dmps[[tissue]] <- (list("f" = f, "m" = length(common)))
# }
# saveRDS(res_hypergeometric_dvps, paste0(project_path,"Tissues/enrichment_mcpgs_dvps.fisher.rds"))
