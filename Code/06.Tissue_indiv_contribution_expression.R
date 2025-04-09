# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio
# setwd("..")
first_dir <- "/gpfs/"
setwd(paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Methylation/"))

# -------------- #
print(Sys.time())
#-------------- #

tissues <- list.dirs("Tissues/", full.names = F)[-1]
# tissues <- tissues[tissues!="KidneyCortex"]
# tissues <- tissues[33:length(tissues)]

#### reading metadata ####
metadata <- lapply(tissues, function(tissue) readRDS(paste0("Tissues/", tissue, "/metadata_expression.rds")))
names(metadata) <- tissues
metadata_m <- lapply(tissues, function(tissue) readRDS(paste0("Tissues/", tissue, "/metadata.rds")))
names(metadata_m) <- tissues

#### reading beta values ####
beta <- lapply(tissues, function(tissue) readRDS(paste0("Tissues/", tissue, "/counts.rds")))
names(beta) <- tissues

for (tissue in tissues) {
  metadata[[tissue]][,'Tissue'] <- tissue
}

### reading admixture and joining metadata #### 
metadata_df <- do.call("rbind",lapply(tissues, function(tissue) metadata[[tissue]][,c("Donor", "Tissue","RIN","ExonicRate","IschemicTime", "HardyScale")]))
head(metadata_df)
rownames(metadata_df) <- paste0(metadata_df$Donor,':',metadata_df$Tissue)
metadata_df$HardyScale <- droplevels(metadata_df$HardyScale)
metadata_df$IschemicTime <- as.numeric(metadata_df$IschemicTime)

metadata_df_m <- do.call("rbind",lapply(tissues, function(tissue) metadata_m[[tissue]][,c("SUBJID", "Tissue","PEER1", "PEER2", "PEER3", "PEER4", "PEER5")]))
head(metadata_df_m)
rownames(metadata_df_m) <- paste0(metadata_df_m$SUBJID,':',metadata_df_m$Tissue)

metadata_df <- metadata_df[rownames(metadata_df %in% rownames(metadata_df_m)),]

### create label with both subj-id pairs ####
genes <- list()
for (tissue in tissues) {
  genes[[tissue]] <- rownames(beta[[tissue]])
}

genes_all <- Reduce(intersect, genes)
head(genes_all)

for (tissue in tissues) {
  colnames(beta[[tissue]]) <- paste0(colnames(beta[[tissue]]),':',tissue)
  beta[[tissue]] <- beta[[tissue]][genes_all,]
}

### merging all data together ####
beta_df <- do.call("cbind",lapply(tissues, function(tissue) beta[[tissue]]))
beta_df <- beta_df[,rownames(metadata_df)]
head(beta_df)

### performing linear model #### 

print('data prepared')

### run model ####
library('variancePartition')
library('BiocParallel')
# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param = SnowParam(40, "SOCK", progressbar=TRUE)
register(param)

metadata_df$RIN <- scale(metadata_df$RIN)
metadata_df$ExonicRate <- scale(metadata_df$ExonicRate)
metadata_df$IschemicTime <- scale(metadata_df$IschemicTime)

print('create model')
form <- ~ RIN+ExonicRate+IschemicTime+(1|HardyScale) + (1 | Donor) + (1 | Tissue)
print(form)

### run in chuncks
# dfs <- split(as.data.frame(beta_df), (seq(nrow(beta_df))-1) %/% 50000) 
# for (i in c(1:length(dfs))) {
#   print('variance Partition') 
#   print(paste0('chunck ',i))
#   
#   vp = fitExtractVarPartModel(as.matrix(dfs[[i]]), form, metadata_df)
#   pl <- plotVarPart( sortCols(vp))
#   
#   saveRDS(vp, paste0(i,'_chunck_var_part.rds'))
#   
#   pdf(file = paste0("Plots/",i,"_chunckvar_part.pdf"), w = 6, h = 3.5)
#   print(pl)
#   dev.off()
# }
print('variance Partition')
vp = fitExtractVarPartModel(beta_df, form, metadata_df)
pl <- plotVarPart( sortCols(vp))

saveRDS(vp, 'var_part_expression.rds')

pdf(file = paste0("Plots/var_part_expr.pdf"), w = 6, h = 3.5)
print(pl)
dev.off()

