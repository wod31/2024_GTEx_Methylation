#### enrichments #####
#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Get funtional enrichments of DMPs for GO pathways
# @software version: R=4.2.2

library(missMethyl)
library(ggplot2)

### load ranges 
first_dir <- "/gpfs/projects/bsc83/"

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

### perform enrichments ##########
# GO
traits_to_use <- c('EURv1','SEX2','AGE','BMI')
# all DMRs

GOenrichments <- list()

for (tissue in tissues) {
  print(tissue)
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

  }
}

saveRDS(GOenrichments, 'GO_enrichments_DMP_peers5.hyper.continous.rds')

for (tissue in tissues) {
  print(tissue)
  GOenrichments[[tissue]] <- list()
  for (trait in names(results_DML[[tissue]])[!is.na(results_DML[[tissue]])]) {
    print(trait)
    if (!trait %in% traits_to_use) {
      next}
    if (nrow(results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$adj.P.Val<0.05 & results_DML[[tissue]][[trait]]$logFC<0,]) < 1 ) {
      next
    }
    tryCatch(
      {res <- gometh(rownames(results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$adj.P.Val<0.05 & results_DML[[tissue]][[trait]]$logFC<0,]), all.cpg=rownames(results_DML[[tissue]][[trait]]),
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
  
  }
}

saveRDS(GOenrichments, 'GO_enrichments_DMP_peers5.hypo.continous.rds')

for (tissue in tissues) {
  print(tissue)
  GOenrichments[[tissue]] <- list()
  for (trait in names(results_DML[[tissue]])[!is.na(results_DML[[tissue]])]) {
    print(trait)
    if (!trait %in% traits_to_use) {
      next}
    if (nrow(results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$adj.P.Val<0.05,]) < 1 ) {
      next
    }
    tryCatch(
      {res <- gometh(rownames(results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$adj.P.Val<0.05,]), all.cpg=rownames(results_DML[[tissue]][[trait]]),
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
  
  }
}

saveRDS(GOenrichments, 'GO_enrichments_DMP_peers5.continous.rds')

