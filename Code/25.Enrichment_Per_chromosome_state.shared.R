
# library(minfi)
# library(lumi)
library(dplyr)
library(tidyr)
library(tidyverse)
## Input data #####
#library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Age", "Ancestry", "BMI", "Sex")
traits_to_use <- c('EURv1','SEX2','AGE','BMI')

#### enrichment -- chi-squared first ####
### read sharing ####
sharing <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Sharing_DMP.rds')
shared_cpgs <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/chromHMM_shared_9tissues.rds')

betas <- read.csv("~/marenostrum_scratch/MN4/bsc83/bsc83535/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")$Name


### enrichment shared positions
my_fisher <- function(type, trait){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  res <- sharing[sharing$trait == trait,]
  chrom_tissue <- shared_cpgs

  type_df <- chrom_tissue[chrom_tissue$region_chromhmm == type,]
  type_diff <- (type_df[type_df$name_ann %in% res$CG[res$number>=3 & res$dir==-1],])

  
  if (nrow(type_diff) < 1 ) {
    return(NA)
  }
  tryCatch(
    {res_2 <- gometh(type_diff$name_ann, all.cpg=chrom_tissue$name_ann,
                   collection="GO", array.type="EPIC")
    res <- res#[res$ONTOLOGY=="BP",]
    print(table(res$FDR<0.05))
    return(res)
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
families <- c('Enh','EnhBiv','Het','Quies','ReprPC','TSS','TssBiv','Tx','ZNF/Rpts')
fisher_results <- lapply(c("EURv1" ,"SEX2" , "AGE"  , "BMI"), function(trait) lapply(families, function(region) my_fisher(region,trait)))
names(fisher_results) <- c("Ancestry" ,"Sex" , "Age"  , "BMI")

for (name in c("Ancestry" ,"Sex" , "Age"  , "BMI")) {
  names(fisher_results[[name]]) <- families
}

saveRDS(fisher_results, '/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/Tissues/enrichment_chromhmm_hypo_shared_CI.continous.rds')


library(missMethyl)

betas <- read.csv("~/marenostrum_scratch/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")$Name
to_plot$number <- as.numeric(to_plot$number)
res <- gometh(to_plot$CG[to_plot$trait=='AGE' & to_plot$number>=3], all.cpg=(betas),
              collection="GO", array.type="EPIC")
res <- res[res$ONTOLOGY=="BP",]
print(table(res$FDR<0.05))

topgo <- topGSA(res, n=20)
(ggplot(data = topgo[topgo$FDR<0.05,], aes(x=DE/N, y = factor(TERM, levels=rev(TERM)),
                                           color = -log10(FDR), size = DE)) +
    geom_point() + scale_color_gradient(low = "red", high = "blue") +
    theme_bw() +  ylab("") +  xlab("Gene Ratio") +
    ggtitle(paste0('AGE Shared'," GO BP sig")))

## enrichment DMPs individually variable 
sharing <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Sharing_DMP.rds')