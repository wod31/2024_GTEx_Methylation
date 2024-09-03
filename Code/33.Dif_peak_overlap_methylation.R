#### merge dna methylation with differential peak calling

#first_dir <- "/gpfs/projects/bsc83/"
first_dir <- "~/marenostrum/"
#annotation <- read.csv("/gpfs/scratch/bsc83/bsc83535/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")
annotation <- read.csv("~/marenostrum_scratch/MN4/bsc83/bsc83535/GTEx/v9/Oliva/GPL21145_MethylationEPIC_15073387_v-1-0_processed.csv")

project_path <- paste0(first_dir, "Projects/GTEx_v8/Methylation/")

tissues <- c("MuscleSkeletal",'KidneyCortex',"BreastMammaryTissue","Lung", "ColonTransverse", "WholeBlood")
names <- c("Sex")
traits_to_use <- c('SEX2')

results_DML <- lapply(tissues, function(tis) 
  readRDS(paste0("~/marenostrum/Projects/GTEx_v8/Methylation/Tissues/",tis,"/DML_results_5_PEERs_continous.rds")))
names(results_DML) <- tissues

### overlap chromhmm and cpgs tested #### 
# ann_granges <- annotation[!is.na(annotation$MAPINFO),] %>%
#   dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>% distinct() %>% 
#   makeGRangesFromDataFrame(keep.extra.columns=T)
ann_bed <- annotation[!is.na(annotation$MAPINFO) & !is.na(annotation$CHR),] %>%
  dplyr::select(chrom=CHR, start=MAPINFO, end=MAPINFO, name=IlmnID) %>% distinct()
ann_bed$chrom <- paste0('chr',ann_bed$chrom)
head(ann_bed)
ann_bed$start <- ann_bed$start-1

## read differential binding
dif_bind <- read.delim('~/marenostrum/Projects/GTEx_v8/Methylation/Data/male_vs_female_chip_heart.bed', header = F)

library(valr)
colnames(dif_bind) <- c('chrom','start','end','Normalized_count_male','Normalized_count_female','Log2FC','logPval','logQval')
intersect <- bed_intersect(ann_bed, dif_bind, suffix = c("_ann", "_chip"))
head(intersect)
intersect <- intersect[intersect$.overlap ==1,]

### read sharing <- 
sharing <- readRDS('~/marenostrum/Projects/GTEx_v8/Methylation/Data/Sharing_DMP.rds')

### overlap with DMA 
shared_chip <- merge(sharing[sharing$trait=='SEX2',], intersect, by.x='CG', by.y='name_ann')
