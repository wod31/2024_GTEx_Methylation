#### enrichments #####

library(missMethyl)
library(ggplot2)

### load ranges 
first_dir <- "/gpfs/projects/bsc83/"
#first_dir <- "~/marenostrum/"

project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Age", "Ancestry", "BMI", "Sex")

results_DML <- lapply(tissues, function(tis) 
  readRDS(paste0(project_path,"Tissues/",tis,"/DML_results_5_PEERs_continous.rds")))
names(results_DML) <- tissues

tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Methylation/Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

#### differentiate between up and down 

hypo <- lapply(tissues, function(tissue) sapply(names(results_DML[[tissue]])[!is.na(results_DML[[tissue]])], function(trait)
  results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$logFC<0,]))
names(hypo) <- tissues
hyper <- lapply(tissues, function(tissue) sapply(names(results_DML[[tissue]])[!is.na(results_DML[[tissue]])], function(trait)
  results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$logFC>0,]))
names(hyper) <- tissues

### perform enrichments 
# GO
traits_to_use <- c('EURv1','SEX2','AGE','BMI')
# all DMRs

GOenrichments <- list()
#GOenrichments_prom <- list()

#annotation <- read.csv("/gpfs/scratch/bsc83/bsc83535/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")

for (tissue in tissues) {
  print(tissue)
  #betas <- rownames(readRDS(paste0(project_path,'Tissues/', tissue, "/data.rds")))
  GOenrichments[[tissue]] <- list()
  for (trait in names(results_DML[[tissue]])[!is.na(results_DML[[tissue]])]) {
    print(trait)
    if (!trait %in% traits_to_use) {
      next}
    if (nrow(results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$adj.P.Val<0.05 & results_DML[[tissue]][[trait]]$logFC>0,]) < 1 ) {
      next
    }
    tryCatch(
      {res <- gometh(rownames(results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$adj.P.Val<0.05 & results_DML[[tissue]][[trait]]$logFC>0,]), all.cpg=rownames(results_DML[[tissue]][[trait]]),
                           collection="GO", array.type="EPIC")
      res <- res[res$ONTOLOGY=="BP",]
      print(table(res$FDR<0.05))
      if (sum(res$FDR<0.05) > 0) {
            GOenrichments[[tissue]][[trait]] <- res[res$FDR<0.05,]
          }
      },  error=function(cond) {
        message("Error")
        message("Here's the original error message:")
        message(cond)
        # Choose a return value in case of error
        return(NA)
      })
    # if (!is.na(nrow(gst.region))) {
    #   gst.region <- gst.region[gst.region$ONTOLOGY=="BP",]
    #   print(table(gst.region$FDR<0.05))
    #
    #   if (sum(gst.region$FDR<0.05) > 0) {
    #     GOenrichments[[tissue]][[trait]] <- gst.region[gst.region$FDR<0.05,]
    #   }
    # }

    #topgo <- topGSA(gst.region, n=20)

    # pdf(paste0(project_path, "Plots/",tissue,"_",trait,"_GO_DMP.pdf"), width = 8, height = 5)
    # print(ggplot(data = topgo, aes(x=DE/N, y = factor(TERM, levels=rev(TERM)),
    #                                color = -log10(FDR), size = DE)) +
    #         geom_point() + scale_color_gradient(low = "red", high = "blue") +
    #         theme_bw() +  ylab("") +  xlab("Gene Ratio") +
    #         ggtitle(paste0(tissue,' ',trait," GO top all")))
    # dev.off()
    #
    # gst.region_prom <- gometh(rownames(results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$adj.P.Val<0.05,]), all.cpg=(betas),
    #                      collection="GO", array.type="EPIC")#, genomic.features = c("TSS200",
    #                                                                               #"TSS1500",
    #                                                                               #"1stExon"))
    # gst.region_prom <- gst.region_prom[gst.region_prom$ONTOLOGY=="BP",]
    # print(table(gst.region_prom$FDR<0.05))
    #
    # if (sum(gst.region_prom$FDR<0.05) > 0) {
    #   GOenrichments_prom[[tissue]][[trait]] <- gst.region_prom[gst.region_prom$FDR<0.05,]
    # }
    #
    # topgo <- topGSA(gst.region_prom, n=20)
    #
    # pdf(paste0(project_path, "Plots/",tissue,"_",trait,"_GO_DMR.prom.pdf"), width = 8, height = 5)
    # print(ggplot(data = topgo, aes(x=DE/N, y = factor(TERM, levels=rev(TERM)),
    #                                color = -log10(FDR), size = DE)) +
    #         geom_point() + scale_color_gradient(low = "red", high = "blue") +
    #         theme_bw() +  ylab("") +  xlab("Gene Ratio") +
    #         ggtitle(paste0(tissue,' ',trait," GO top all")))
    # dev.off()

  }
}

saveRDS(GOenrichments, 'GO_enrichments_DMP_peers5.hyper.continous.rds')

# for (tissue in tissues) {
#   print(tissue)
#   betas <- rownames(readRDS(paste0(project_path,'Tissues/', tissue, "/data.rds")))
#   GOenrichments[[tissue]] <- list()
#   for (trait in names(results_DML[[tissue]])[!is.na(results_DML[[tissue]])]) {
#     print(trait)
#     if (!trait %in% traits_to_use) {
#       next}
#     if (nrow(results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$adj.P.Val<0.05 & results_DML[[tissue]][[trait]]$logFC<0,]) < 1 ) {
#       next
#     }
#     tryCatch(
#       {res <- gometh(rownames(results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$adj.P.Val<0.05 & results_DML[[tissue]][[trait]]$logFC<0,]), all.cpg=(betas),
#                      collection="GO", array.type="EPIC")
#       res <- res[res$ONTOLOGY=="BP",]
#       print(table(res$FDR<0.05))
#       if (sum(res$FDR<0.05) > 0) {
#         GOenrichments[[tissue]][[trait]] <- res[res$FDR<0.05,]
#       }
#       },  error=function(cond) {
#         message("Error")
#         message("Here's the original error message:")
#         message(cond)
#         # Choose a return value in case of error
#         return(NA)
#       })
#     # if (!is.na(nrow(gst.region))) {
#     #   gst.region <- gst.region[gst.region$ONTOLOGY=="BP",]
#     #   print(table(gst.region$FDR<0.05))
#     #
#     #   if (sum(gst.region$FDR<0.05) > 0) {
#     #     GOenrichments[[tissue]][[trait]] <- gst.region[gst.region$FDR<0.05,]
#     #   }
#     # }
# 
#     #topgo <- topGSA(gst.region, n=20)
# 
#     # pdf(paste0(project_path, "Plots/",tissue,"_",trait,"_GO_DMP.pdf"), width = 8, height = 5)
#     # print(ggplot(data = topgo, aes(x=DE/N, y = factor(TERM, levels=rev(TERM)),
#     #                                color = -log10(FDR), size = DE)) +
#     #         geom_point() + scale_color_gradient(low = "red", high = "blue") +
#     #         theme_bw() +  ylab("") +  xlab("Gene Ratio") +
#     #         ggtitle(paste0(tissue,' ',trait," GO top all")))
#     # dev.off()
#     #
#     # gst.region_prom <- gometh(rownames(results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$adj.P.Val<0.05,]), all.cpg=(betas),
#     #                      collection="GO", array.type="EPIC")#, genomic.features = c("TSS200",
#     #                                                                               #"TSS1500",
#     #                                                                               #"1stExon"))
#     # gst.region_prom <- gst.region_prom[gst.region_prom$ONTOLOGY=="BP",]
#     # print(table(gst.region_prom$FDR<0.05))
#     #
#     # if (sum(gst.region_prom$FDR<0.05) > 0) {
#     #   GOenrichments_prom[[tissue]][[trait]] <- gst.region_prom[gst.region_prom$FDR<0.05,]
#     # }
#     #
#     # topgo <- topGSA(gst.region_prom, n=20)
#     #
#     # pdf(paste0(project_path, "Plots/",tissue,"_",trait,"_GO_DMR.prom.pdf"), width = 8, height = 5)
#     # print(ggplot(data = topgo, aes(x=DE/N, y = factor(TERM, levels=rev(TERM)),
#     #                                color = -log10(FDR), size = DE)) +
#     #         geom_point() + scale_color_gradient(low = "red", high = "blue") +
#     #         theme_bw() +  ylab("") +  xlab("Gene Ratio") +
#     #         ggtitle(paste0(tissue,' ',trait," GO top all")))
#     # dev.off()
# 
#   }
# }
# 
# saveRDS(GOenrichments, 'GO_enrichments_DMP_peers5.hypo.continous.rds')

for (tissue in tissues) {
  print(tissue)
  betas <- rownames(readRDS(paste0(project_path,'Tissues/', tissue, "/data.rds")))
  GOenrichments[[tissue]] <- list()
  for (trait in names(results_DML[[tissue]])[!is.na(results_DML[[tissue]])]) {
    print(trait)
    if (!trait %in% traits_to_use) {
      next}
    if (nrow(results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$adj.P.Val<0.05,]) < 1 ) {
      next
    }
    tryCatch(
      {res <- gometh(rownames(results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$adj.P.Val<0.05,]), all.cpg=(betas),
                     collection="GO", array.type="EPIC")
      res <- res[res$ONTOLOGY=="BP",]
      print(table(res$FDR<0.05))
      if (sum(res$FDR<0.05) > 0) {
        GOenrichments[[tissue]][[trait]] <- res[res$FDR<0.05,]
      }
      },  error=function(cond) {
        message("Error")
        message("Here's the original error message:")
        message(cond)
        # Choose a return value in case of error
        return(NA)
      })
    # if (!is.na(nrow(gst.region))) {
    #   gst.region <- gst.region[gst.region$ONTOLOGY=="BP",]
    #   print(table(gst.region$FDR<0.05))
    #
    #   if (sum(gst.region$FDR<0.05) > 0) {
    #     GOenrichments[[tissue]][[trait]] <- gst.region[gst.region$FDR<0.05,]
    #   }
    # }

    #topgo <- topGSA(gst.region, n=20)

    # pdf(paste0(project_path, "Plots/",tissue,"_",trait,"_GO_DMP.pdf"), width = 8, height = 5)
    # print(ggplot(data = topgo, aes(x=DE/N, y = factor(TERM, levels=rev(TERM)),
    #                                color = -log10(FDR), size = DE)) +
    #         geom_point() + scale_color_gradient(low = "red", high = "blue") +
    #         theme_bw() +  ylab("") +  xlab("Gene Ratio") +
    #         ggtitle(paste0(tissue,' ',trait," GO top all")))
    # dev.off()
    #
    # gst.region_prom <- gometh(rownames(results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$adj.P.Val<0.05,]), all.cpg=(betas),
    #                      collection="GO", array.type="EPIC")#, genomic.features = c("TSS200",
    #                                                                               #"TSS1500",
    #                                                                               #"1stExon"))
    # gst.region_prom <- gst.region_prom[gst.region_prom$ONTOLOGY=="BP",]
    # print(table(gst.region_prom$FDR<0.05))
    #
    # if (sum(gst.region_prom$FDR<0.05) > 0) {
    #   GOenrichments_prom[[tissue]][[trait]] <- gst.region_prom[gst.region_prom$FDR<0.05,]
    # }
    #
    # topgo <- topGSA(gst.region_prom, n=20)
    #
    # pdf(paste0(project_path, "Plots/",tissue,"_",trait,"_GO_DMR.prom.pdf"), width = 8, height = 5)
    # print(ggplot(data = topgo, aes(x=DE/N, y = factor(TERM, levels=rev(TERM)),
    #                                color = -log10(FDR), size = DE)) +
    #         geom_point() + scale_color_gradient(low = "red", high = "blue") +
    #         theme_bw() +  ylab("") +  xlab("Gene Ratio") +
    #         ggtitle(paste0(tissue,' ',trait," GO top all")))
    # dev.off()

  }
}

saveRDS(GOenrichments, 'GO_enrichments_DMP_peers5.continous.rds')

# #KEGG 
# library(limma)
# library(DMRcate)
# library(ExperimentHub)
# # 
# KEGG_enrichments <- list()
# ## I Download KEGG because missmethyl function goregion with "KEGG" option fails
# kegg <- readRDS(paste0(project_path,"Data/Hs.c2.cp.kegg.v7.1.entrez.rds"))
# 
# for (tissue in tissues) {
#   print(tissue)
#   betas <- rownames(readRDS(paste0(project_path,'Tissues/', tissue, "/data.rds")))
#   KEGG_enrichments[[tissue]] <- list()
#   for (trait in names(results_DML[[tissue]])[!is.na(results_DML[[tissue]])]) {
#     print(trait)
#     if (!trait %in% traits_to_use) {
#       next}
#     if (nrow(results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$adj.P.Val<0.05 & results_DML[[tissue]][[trait]]$logFC<0,]) < 1 ) {
#       next
#     }
#     tryCatch(
#       {res <- gsameth(rownames(results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$adj.P.Val<0.05 & results_DML[[tissue]][[trait]]$logFC<0,]), all.cpg=(betas),
#                collection=kegg,array.type="EPIC")
#       print(table(res$FDR<0.05))
#       if (sum(res$FDR<0.05) > 0) {
#             KEGG_enrichments[[tissue]][[trait]] <- res[res$FDR<0.05,]
#           }
#       },  error=function(cond) {
#         message("Error")
#         message("Here's the original error message:")
#         message(cond)
#         # Choose a return value in case of error
#         return(NA)
#       })
#     # if (!is.na(nrow(gst.region))) {
#     #   print(table(gsa.region$FDR<0.05))
#     #   
#     #   if (sum(gst.region$FDR<0.05) > 0) {
#     #     KEGG_enrichments[[tissue]][[trait]] <- gsa.region[gsa.region$FDR<0.05,]
#     #   }
#     # }
# 
#     # topgo <- topGSA(gsa.region, n=20)
#     # topgo$Description <- rownames(topgo)
# 
#     # pdf(paste0(project_path, "Plots/",tissue,"_",trait,"_KEGG_DMP.pdf"), width = 8, height = 5)
#     # print(ggplot(data = topgo, aes(x=DE/N, y = factor(Description, levels=rev(Description)),
#     #                                color = -log10(FDR), size = DE)) +
#     #         geom_point() + scale_color_gradient(low = "red", high = "blue") +
#     #         theme_bw() +  ylab("") +  xlab("Gene Ratio") +
#     #         ggtitle(paste0(tissue,' ',trait," KEGG top all")))
#     # dev.off()
# 
#   }
# }
# 
# saveRDS(KEGG_enrichments, 'KEGG_enrichments_DMP.hypo.rds')
# 
# #hallmark
# 
# hallmark_enrichments <- list()
# hallmark <- readRDS(paste0(project_path,"Data/Hs.h.all.v7.1.entrez.rds"))
# 
# for (tissue in tissues[5:9]) {
#   print(tissue)
#   betas <- rownames(readRDS(paste0(project_path,'Tissues/', tissue, "/data.rds")))
#   hallmark_enrichments[[tissue]] <- list()
#   for (trait in names(results_DML[[tissue]])[!is.na(results_DML[[tissue]])]) {
#     print(trait)
#     if (!trait %in% traits_to_use) {
#       next}
#     if (nrow(results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$adj.P.Val<0.05 & results_DML[[tissue]][[trait]]$logFC<0,]) < 1 ) {
#       next
#     }
#     tryCatch(
#       {res <- gsameth(rownames(results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$adj.P.Val<0.05 & results_DML[[tissue]][[trait]]$logFC<0,]), all.cpg=(betas),
#               collection=hallmark,array.type="EPIC")
#       print(table(res$FDR<0.05))
#       if (sum(res$FDR<0.05) > 0) {
#         hallmark_enrichments[[tissue]][[trait]] <- res[res$FDR<0.05,]}
#       },  error=function(cond) {
#         message("Error")
#         message("Here's the original error message:")
#         message(cond)
#         # Choose a return value in case of error
#         #return(NA)
#       })
#     # if (!is.na(nrow(gst.region))) {
#     #   print(table(gsa.region$FDR<0.05))
#     #   
#     #   if (sum(gst.region$FDR<0.05) > 0) {
#     #     hallmark_enrichments[[tissue]][[trait]] <- gsa.region[gsa.region$FDR<0.05,]
#     #   }
#     # }
#     # topgo <- topGSA(gsa.region, n=20)
#     # topgo$Description <- rownames(topgo)
# 
#     # pdf(paste0(project_path, "Plots/",tissue,"_",trait,"_Hallmark_DMP.pdf"), width = 8, height = 5)
#     # print(ggplot(data = topgo, aes(x=DE/N, y = factor(Description, levels=rev(Description)),
#     #                                color = -log10(FDR), size = DE)) +
#     #         geom_point() + scale_color_gradient(low = "red", high = "blue") +
#     #         theme_bw() +  ylab("") +  xlab("Gene Ratio") +
#     #         ggtitle(paste0(tissue,' ',trait," Hallmark top all")))
#     # dev.off()
# 
#   }
# }
# 
# saveRDS(hallmark_enrichments, 'Hallmarks_enrichments_DMP.hypo.rds')
# 
