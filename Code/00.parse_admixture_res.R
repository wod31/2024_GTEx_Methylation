#!/usr/bin/env Rscript
# @Author: Winona Oliveros Diez
# @E-mail: winn95@gmail.com
# @Description: Parse Admixture annotation
# @software version: R=4.2.2


install.packages(c("fields","RColorBrewer","mapplots"))
BiocManager::install("LEA")
library(LEA)

tbl=read.table("~/marenostrum_scratch/GTEx/v8/genotype_data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phasedLDPrunechr21.2.Q")
tbl2=read.table("~/marenostrum_scratch/GTEx/v8/genotype_data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phasedLDPrunechr21.3.Q")
tbl3=read.table("~/marenostrum_scratch/GTEx/v8/genotype_data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phasedLDPrunechr21.4.Q")

barplot(t(as.matrix(tbl2)), col=rainbow(3),
          xlab="Individual #", ylab="Ancestry", border=NA)


#match sample order
samps<-read.table("~/marenostrum_scratch/GTEx/v8/genotype_data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phasedLDPrunechr21.fam")[,1]
head(samps)

#### ancestries #### 
sample_metadata <- read.table('~/marenostrum_scratch/GTEx/v8/metadata/inferred_ancestry_838donors.txt', header = T)
head(sample_metadata)

rownames(sample_metadata) <- sample_metadata$ID
sample_metadata <- sample_metadata[samps,]
rownames(tbl) <- samps
rownames(tbl2) <- samps
rownames(tbl3) <- samps
sample_metadata <- sample_metadata[order(sample_metadata$inferred_ancestry),]

barplot(t(as.matrix(tbl2[sample_metadata$ID,])), col=rainbow(3),
        xlab="Individual #", ylab="Ancestry", border=NA, )

#### compare with and without asians ####
tbl2$ID <- rownames(tbl2)
common <- merge(tbl, tbl2, by='ID')

ggplot(common, aes(x = V1.x, y = V1.y, color=inferred_ancestry)) + 
  geom_point(alpha=0.6, size=3)

library("ggpubr")
ggscatter(common, x = "V2.x", y = "V2.y", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "EUR v1", ylab = "EUR v2")

#### explore ancestries compared to gtex #####
tbl_f <- as.data.frame(tbl2)
head(tbl_f)
tbl_f$ID <- rownames(tbl_f)
tbl_f <- merge(tbl_f, sample_metadata)
head(tbl_f)

# ggplot(tbl_f, aes(x = V1, y = V2, color=inferred_ancestry)) + 
#   geom_point(alpha=0.6, size=3)

ggplot(tbl_f,aes(x = V3, y = inferred_ancestry, fill=inferred_ancestry)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()+
  theme_classic()

library(reshape2)
tbl_m <- melt(tbl_f, id.vars = c("ID", "inferred_ancestry"))
head(tbl_m)

pdf('marenostrum/Projects/GTEx_v8/Methylation/Plots/admixture_k4.pdf')
ggplot(tbl_m,aes(x = value, y = inferred_ancestry, fill=variable, color=variable)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, position = position_jitterdodge(), alpha=0.5)+
  theme_classic()
dev.off()

library(plotly)

plot_ly(x=tbl_f$V1, y=tbl_f$V2, z=tbl_f$V3, type="scatter3d", mode="markers", color=tbl_f$inferred_ancestry)


#### remove asians 
samps<-read.table("~/marenostrum_scratch/GTEx/v8/genotype_data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phasedLDPrunechr21.fam")
head(samps)
samps <- samps[samps$V1 %in% tbl$ID[tbl$inferred_ancestry!='ASN'],]
write.table(samps,'~/marenostrum_scratch/GTEx/v8/genotype_data/samples_to_keep.txt', sep = '\t', col.names = F, row.names = F, quote = F)

#### save results #### 
head(tbl)
write.table(common,'~/marenostrum_scratch/GTEx/v8/genotype_data/admixture_inferred_ancestry.txt', sep = '\t', col.names = F, row.names = F, quote = F)
