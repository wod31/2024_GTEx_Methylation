options(stringsAsFactors = F)
set.seed(1)

# Libraries ####
library(Rtsne)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(scales)

first_dir <- "/gpfs/"
setwd(paste0(first_dir, "/projects/bsc83/Projects/GTEx_v8/Methylation/"))# Paths ####
gtex_path <- "/gpfs/projects/bsc83/MN4/bsc83/Projects/GTEx_v8/Raquel/" 

# -------------- #
print(Sys.time())
#-------------- #

tissues <- list.dirs("Tissues/", full.names = F)[-1]
# tissues <- tissues[tissues!="KidneyCortex"]
# tissues <- tissues[33:length(tissues)]
# Filenames ####
TissueInfoFile <- "Tissue_info.rds"
TissueInfo <- readRDS(paste0(gtex_path,TissueInfoFile))
TissueInfo <- TissueInfo[TissueInfo$tissue_ID %in% tissues,]
tissues_cols <- TissueInfo$colcodes
names(tissues_cols) <- TissueInfo$tissue_ID

#### reading metadata ####
metadata <- lapply(tissues, function(tissue) readRDS(paste0("Tissues/", tissue, "/metadata.rds")))
names(metadata) <- tissues

#### reading beta values ####
beta <- lapply(tissues, function(tissue) readRDS(paste0("Tissues/", tissue, "/data.rds")))
names(beta) <- tissues

for (tissue in tissues) {
  metadata[[tissue]][,'Tissue'] <- tissue
}

### reading admixture and joining metadata #### 
metadata_df <- do.call("rbind",lapply(tissues, function(tissue) metadata[[tissue]][,c("SUBJID", "Tissue","PEER1", "PEER2", "PEER3", "PEER4", "PEER5")]))
head(metadata_df)
rownames(metadata_df) <- paste0(metadata_df$SUBJID,':',metadata_df$Tissue)

### create label with both subj-id pairs ####
for (tissue in tissues) {
  colnames(beta[[tissue]]) <- paste0(colnames(beta[[tissue]]),':',tissue)
}

### merging all data together ####
# beta_df <- do.call("cbind",lapply(tissues, function(tissue) beta[[tissue]]))
# 
# probes <- rownames(beta_df)
# beta <- sapply(beta_df, as.numeric)
# rownames(beta) <- probes
# 
# ## select most variable ####
# cv <- apply(beta, 1, var)
# summary(cv)
# 
# ## Select 10% of most var
# cutoff <- 1
# summary(cv[rank(cv) / length(cv) > 1 - cutoff])
# 
# meth_variable <- beta[rank(cv) / length(cv) > 1 - cutoff, ]
# dim(meth_variable)
# 
# #M_vals <- log2(meth_variable/(1-meth_variable))
# 
# #par(mfrow=c(1,1))
# 
# # tsne <- Rtsne(t(meth_variable),#
# #               dims = 2, perplexity=30, verbose=TRUE, max_iter = 1000,partial_pca = T)
# # pdf( paste0('Plots/', "tSNE.Methylationprofile.coloured_by_tissue.50var.pdf"),
# #      width = 7.5,
# #      height = 7.5)
# # par(mfrow=c(1,1))
# # plot(tsne$Y, main="Expression profiles",pch=16,
# #      col=alpha( tissues_cols[sapply(colnames(meth_variable), function(sample) unlist(strsplit(sample,split = ":"))[[2]] )],0.5),cex=0.75,
# #      xlab = "t-SNE1",ylab="t-SNE2")
# # dev.off()
# library(ggpubr)
# library(RColorBrewer)
# library(viridis)
# library(gridExtra)
# #library(FactoMineR)
# #library(factoextra)
# 
# # Merge with sample metadata --
# meth_variable <- as.data.frame(t(meth_variable))
# # n <- dim(meth_variable)[2]
# # meth_variable$SUBJID <- rownames(meth_variable)
# # meth_variable <- cbind.data.frame(meth_variable, metadata_df[meth_variable$SUBJID,1:3])
# # Run PCA --
# #res.pca <- PCA(exprs.data, ncp = 100, graph = F, quali.sup = c((n+1):ncol(meth_variable)) )
# 
# pca_res <- prcomp((meth_variable), scale. = TRUE)
# # Results for individuals --
# res.ind <- pca_res$x
# pca <- as.data.frame(res.ind[,1:2])
# #pca$SUBJID <- rownames(pca)
# #pca <- merge(pca, metadata_df, by = "SUBJID")
# 
# svg( paste0('Plots/', "PCA.Methylationprofile.coloured_by_tissue.all.svg"),
#      width = 7.5,
#      height = 7.5)
# par(mfrow=c(1,1))
# plot(pca, main="methylation profiles",pch=16,
#      col=alpha( tissues_cols[sapply(rownames(meth_variable), function(sample) unlist(strsplit(sample,split = ":"))[[2]] )],0.5),cex=0.75,
#      xlab = "PCA1",ylab="PCA2")
# dev.off()
# 
# # saveRDS(pca,'PCA_first2.rds')
# 
# ### do hierarchical clustering ####
# M_vals <- log2(meth_variable/(1-meth_variable))
# M_vals_scaled <- scale((meth_variable))
# dist_mat <- dist(M_vals_scaled, method = 'euclidean')
# 
# hclust_avg <- hclust(dist_mat, method = 'average')
# clusterCut <- cutree(hclust_avg, 9)
# #table(clusterCut, metadata_df$tissue)
# library(rafalib)
# pdf( paste0('Plots/', "hierarchicalClustrer.all.pdf"),
#      width = 15,
#      height = 7.5)
# par(mfrow=c(1,1))
# myplclust(hclust_avg, lab = sapply(rownames(meth_variable), function(sample) unlist(strsplit(sample,split = ":"))[[1]] ), 
#           lab.col = tissues_cols[sapply(rownames(meth_variable), function(sample) unlist(strsplit(sample,split = ":"))[[2]] )])
# dev.off()

### plot methylation levels ####
beta_df <- list()
for (tissue in tissues) {
  #colnames(beta[[tissue]]) <- paste0(colnames(beta[[tissue]]),':',tissue)
  probes <- rownames(beta[[tissue]])
  print(tissue)
  print(head(beta[[tissue]]))
  beta[[tissue]] <- sapply(beta[[tissue]], as.numeric)
  rownames(beta[[tissue]]) <- probes
  means <- rowMeans(beta[[tissue]])
  beta_df[[tissue]] <- as.data.frame(means)
  beta_df[[tissue]]$cpg <- rownames(beta_df[[tissue]])
  beta_df[[tissue]]$tis <- tissue
  colnames(beta_df[[tissue]]) <- c('beta','cpg','tissue')
  
}

### merging all data together ####
to_plot <- do.call("rbind.data.frame",beta_df)
#to_plot %>%  group_by(tissue) %>% summarise(n=n(), median=mean(beta)) ->Summary.data
library(ggplot2)
 g <- ggplot(to_plot) + geom_violin(aes(tissue, beta), scale="width", fill="lightgray") + xlab("") + ylab("Beta values") +
   #geom_text(data=Summary.data ,aes(x = tissue, y = 1.1, label=n), fontface =2, size = 3) +
   theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                      panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),)
 ggsave("Plots/methylation_all_tissues.pdf", g, device = "pdf", width = 8, height = 4)
