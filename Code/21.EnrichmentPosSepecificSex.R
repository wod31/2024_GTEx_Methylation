###### Enrichment per region ######
#### enrichments #####

library(missMethyl)
library(ggplot2)
library(dplyr)
library(tidyr)

### load ranges 
first_dir <- "/gpfs/projects/bsc83/"
#first_dir <- "~/marenostrum/"

project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("Lung", "ColonTransverse", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "WholeBlood")
names <- c( "Sex")

results_DML <- lapply(tissues, function(tis) 
  readRDS(paste0(project_path,"Tissues/",tis,"/DML_results_5_PEERs_continous.rds")))
names(results_DML) <- tissues

tissue_info <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Jose/00_Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]

files <- list.files('/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Data/EpiMap/', pattern='.bed.gz',full.names=T)

names_chrom <- c("Lung", "ColonTransverse", "Breast", "MuscleSkeletal", "KidneyCortex", "PBMC")
chromhmm <- lapply(names_chrom, function(tis)
  read.delim(files[grep(tis, files)], sep='\t', header=F))
names(chromhmm) <- tissues

annotation <- read.csv("/gpfs/scratch/bsc83/bsc83535/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")

ann_bed <- annotation[!is.na(annotation$MAPINFO) & !is.na(annotation$CHR),] %>%
  dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>% distinct()
ann_bed$chrom <- paste0('chr',ann_bed$chrom)
head(ann_bed)
ann_bed$start <- ann_bed$start-1

library(valr)
chromhmm_cpgs <- lapply(tissues, function(tis) {
  chrom_df <- chromhmm[[tis]][,c(1:4)]
  colnames(chrom_df) <- c('chrom','start','end','region')
  bed_intersect(ann_bed, chrom_df, suffix = c("_ann", "_chromhmm"))})
names(chromhmm_cpgs) <- tissues

#### differentiate between up and down 

# hypo <- lapply(tissues, function(tissue) sapply(names(results_DML[[tissue]])[!is.na(results_DML[[tissue]])], function(trait)
#   results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$logFC<0,]))
# names(hypo) <- tissues
# hyper <- lapply(tissues, function(tissue) sapply(names(results_DML[[tissue]])[!is.na(results_DML[[tissue]])], function(trait)
#   results_DML[[tissue]][[trait]][results_DML[[tissue]][[trait]]$logFC>0,]))
# names(hyper) <- tissues

### perform enrichments 
# GO
traits_to_use <- c('SEX2')
# all DMRs

GOenrichments <- list()

 for (tissue in tissues) {
   print(tissue)
   GOenrichments[[tissue]] <- list()
   #betas <- rownames(readRDS(paste0(project_path,'Tissues/', tissue, "/data.rds")))
   
   chrom_tissue <- as.data.frame(chromhmm_cpgs[[tissue]])
   chrom_tissue$region_chromhmm_new <- chrom_tissue$region_chromhmm
   # 
   chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TssFlnkD", "TssFlnk", "TssFlnkU","TssA")] <- "TSS"
   chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("EnhA2", "EnhA1","EnhWk","EnhG1", "EnhG2")] <- "Enh"
   
   chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("ReprPCWk","ReprPC")] <- "ReprPC"
   
   chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TxWk","Tx")] <- "Tx"
   
   for (type in c('Enh','EnhBiv','Het','Quies','ReprPC','TSS','TssBiv','Tx','ZNF/Rpts')) {
     #subset betas on the region
   #betas <- betas[betas %in% chrom_tissue$name_ann[chrom_tissue$region_chromhmm_new == type]]
   
    for (trait in names(results_DML[[tissue]])[!is.na(results_DML[[tissue]])]) {
     print(trait)
      if (!trait %in% traits_to_use) {
        next}
     res <- results_DML[[tissue]][[trait]]
     res <- res[rownames(res) %in% chrom_tissue$name_ann[chrom_tissue$region_chromhmm_new == type],]
     #betas <- rownames(res)
     if (nrow(res[res$adj.P.Val<0.05 & res$logFC>0,]) < 1 ) {
       next
     }
     tryCatch(
       {res_enr <- gometh(rownames(res[res$adj.P.Val<0.05 & res$logFC>0,]), all.cpg=rownames(res),
                            collection="GO", array.type="EPIC")
       res_enr <- res_enr[res_enr$ONTOLOGY=="BP",]
       print(table(res_enr$FDR<0.05))
       if (sum(res_enr$FDR<0.05) > 0) {
             GOenrichments[[tissue]][[type]] <- res_enr[res_enr$FDR<0.05,]
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
 }
#
 saveRDS(GOenrichments, 'GO_enrichments_DMP_peers5.hyper.SexChromhmm.rds')

 for (tissue in tissues) {
   print(tissue)
   GOenrichments[[tissue]] <- list()
   #betas <- rownames(readRDS(paste0(project_path,'Tissues/', tissue, "/data.rds")))
   
   chrom_tissue <- as.data.frame(chromhmm_cpgs[[tissue]])
   chrom_tissue <- chrom_tissue[chrom_tissue$.overlap ==1,]
   chrom_tissue$region_chromhmm_new <- chrom_tissue$region_chromhmm
   # 
   chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TssFlnkD", "TssFlnk", "TssFlnkU","TssA")] <- "TSS"
   chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("EnhA2", "EnhA1","EnhWk","EnhG1", "EnhG2")] <- "Enh"
   
   chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("ReprPCWk","ReprPC")] <- "ReprPC"
   
   chrom_tissue$region_chromhmm_new[chrom_tissue$region_chromhmm %in% c("TxWk","Tx")] <- "Tx"
   
   for (type in c('Enh','EnhBiv','Het','Quies','ReprPC','TSS','TssBiv','Tx','ZNF/Rpts')) {
     #subset betas on the region
     #betas <- betas[betas %in% chrom_tissue$name_ann[chrom_tissue$region_chromhmm_new == type]]
     
     for (trait in names(results_DML[[tissue]])[!is.na(results_DML[[tissue]])]) {
       print(trait)
       if (!trait %in% traits_to_use) {
         next}
       res <- results_DML[[tissue]][[trait]]
       res <- res[rownames(res) %in% chrom_tissue$name_ann[chrom_tissue$region_chromhmm_new == type],]
       #betas <- rownames(res)
       if (nrow(res[res$adj.P.Val<0.05,]) < 1 ) {
         next
       }
       tryCatch(
         {res_enr <- gometh(rownames(res[res$adj.P.Val<0.05,]), all.cpg=rownames(res),
                            collection="GO", array.type="EPIC")
         res_enr <- res_enr[res_enr$ONTOLOGY=="BP",]
         print(table(res_enr$FDR<0.05))
         if (sum(res_enr$FDR<0.05) > 0) {
           GOenrichments[[tissue]][[type]] <- res_enr[res_enr$FDR<0.05,]
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
 }
     
# 
 saveRDS(GOenrichments, 'GO_enrichments_DMP_peers5.hypo.SexChromhmm.rds')