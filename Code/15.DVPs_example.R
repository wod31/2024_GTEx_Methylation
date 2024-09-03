rm(list=ls())
Sys.time()
library(ggplot2)
project_path <- "/gpfs/projects/bsc83/Projects/GTEx_v8/Methylation/"
# project_path <- "~/Documents/mn4/projects/bsc83/Projects/GTEx_v8/Methylation/"

tissues <- c("Lung", "ColonTransverse", "Ovary", "Prostate", "BreastMammaryTissue", "MuscleSkeletal", "KidneyCortex", "Testis", "WholeBlood")
names <- c("Age", "Ancestry", "BMI", "Sex")
data <- matrix(nrow = length(tissues), ncol=4, dimnames = list(tissues, names))

first_path <- "/gpfs"
# first_path <- "~/Documents/mn4/"

tissue_info <- readRDS(paste0(first_path, "/projects/bsc83/Projects/GTEx_v8/Jose/00_Data/Tissue_info_whole.rds"))

tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissue_info <- tissue_info[tissue_info$tissue_ID %in% tissues,]


# Ancestry --
# cpgname <- "cg00918796" #9
# cpgname <- 'cg01453458' #9
# cpgname <- 'cg04193820' ### multiple studies etnhicty, autoimmunity, hypomethylated AFR, RNF135 gene, we find it in 4 tissues

#cis-driven shared:
cpgname <- "cg12052018" #3
# cpgname <- "cg00349687" #2


print("About to read betas")
betas <- lapply(tissues, function(tissue) readRDS(paste0(project_path,"Tissues/", tissue, "/data.rds"))[cpgname,])
names(betas) <- tissues

Sys.time()
print("About to read metadata")
metadata <- lapply(tissues, function(tissue) readRDS(paste0(project_path, "Tissues/",tissue, "/metadata.rds")))
names(metadata) <- tissues

print("About to read admixture data")
admixture_ancestry <- read.table(paste0(first_path, '/scratch/bsc83/bsc83535/GTEx/v8/genotype_data/admixture_inferred_ancestry.txt'))
colnames(admixture_ancestry) <- c('SUBJID','AFRv1','EURv1','inferred_ancestry','AFRv2','EURv2')

EA_value <- sapply(tissues, function(tissue)
  median(as.numeric(betas[[tissue]][,colnames(betas[[tissue]]) %in% admixture_ancestry$SUBJID[admixture_ancestry$EURv1>=0.5]]))
)
AA_value <- sapply(tissues, function(tissue)
  median(as.numeric(betas[[tissue]][,colnames(betas[[tissue]]) %in% admixture_ancestry$SUBJID[admixture_ancestry$EURv1<0.5]]))
)

save.image(file = paste0(project_path, "/Plots/DVP_", cpgname, ".RData"))

data <- cbind.data.frame("Ancestry" = c(rep("EA",length(EA_value)),
                                        rep("AA",length(AA_value))),
                         "value" = c(EA_value, AA_value),
                         "Tissue" = rep(tissues, 2))
data$Ancestry <- factor(data$Ancestry, levels = c("EA", "AA"), order = T)
data$Tissue <- factor(data$Tissue, levels = tissues,
                      order = T)
tissue_cols <- tissue_info[tissues, "colcodes"]
names(tissue_cols) <- as.character(levels(data$Tissue))
data$cpg <- rep(cpgname, nrow(data))
data$cpg <- factor(data$cpg)

traits_cols <- c('#C49122','#4B8C61','#70A0DF')
names(traits_cols) <- c('Ancestry','Sex','Age')

p3 <- ggplot(data = data,
             aes(x = Ancestry,
                 y = (value),
                 col = Tissue),
) +
  geom_violin(col = "black") +
  geom_boxplot(col = "black",
               fill = "white",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  geom_point(aes(col = Tissue),
             size = 0.8) +
  geom_line(aes(group=Tissue)) +
  xlab("") +
  ylab("Beta (median)") +
  scale_color_manual(values = tissue_cols) +
  labs(title="") +
  facet_grid(~cpg) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill=traits_cols["Ancestry"]),
        legend.position = "none")
pdf(file = paste0(project_path, "/Plots/DVP_", cpgname, ".pdf"), w = 3, h = 3)
p3
dev.off()


# #### example continous ancestry_tissue
# tissue_data <- as.data.frame(t(betas$Lung)) ## need residuals
# tissue_data$SUBJID <- rownames(tissue_data)
# tissue_data <- merge(tissue_data, admixture_ancestry)
# tissue_data$cg04193820 <- as.numeric(tissue_data$cg04193820)
# library(ggpubr)
# 
# ggscatter(tissue_data, x = "EURv1", y = "cg04193820",
#           add = "reg.line", conf.int = TRUE,
#           cor.coef = TRUE, cor.method = "pearson",
#           xlab = "Prop EUR Ancestry", ylab = "Beta")



#Only one tissue:
tissue <- "ColonTransverse"
to_plot <- cbind(t(betas$ColonTransverse), "EA")
colnames(to_plot) <- c("Betas", "Ancestry")
to_plot[rownames(to_plot) %in% admixture_ancestry$SUBJID[admixture_ancestry$EURv1<=0.25], 2] <- "AA"
to_plot <- as.data.frame(to_plot)
to_plot$Betas <- as.numeric(to_plot$Betas)

to_plot$cpg <- cpgname

ggplot(as.data.frame(to_plot), aes(Ancestry, Betas)) +     
  geom_violin(aes(fill = Ancestry), col = "black", scale = "width") +
  geom_jitter(col = "black",
              alpha = 0.5,
              size = 0.8) +
  theme_bw() +
  facet_grid(~cpg) +
  theme(        
    strip.background = element_rect(fill=traits_cols["Ancestry"]),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", linewidth=1)
)

#Plot 9 facets, one per tissue

beta_values <- c()
tissues <- c()
donors <- c()

for (i in seq_along(betas)) {
  vector <- betas[[i]]
  list_name <- names(betas)[i]
  
  beta_values <- c(beta_values, vector)
  tissues <- c(tissues, rep(list_name, length(vector)))
  donors <- c(donors, names(vector))
}

# Create a data frame
to_plot <- as.data.frame(cbind(
  Betas = beta_values,
  Tissue = tissues,
  SUBJID = donors,
  Ancestry = "EA"
))

to_plot$Ancestry[to_plot$SUBJID %in% admixture_ancestry$SUBJID[admixture_ancestry$EURv1<=0.25]] <- "AA"
to_plot <- as.data.frame(to_plot)
to_plot$Betas <- as.numeric(unlist(to_plot$Betas))
to_plot$Ancestry <- unlist(to_plot$Ancestry)
to_plot$Tissue <- unlist(to_plot$Tissue)
ggplot(to_plot, aes(Ancestry, Betas)) +     
  geom_violin(aes(fill = Ancestry), col = "black", scale = "width") +
  geom_jitter(col = "black",
              alpha = 0.5,
              size = 0.8) +
  theme_bw() +
  facet_grid(~Tissue) +
  theme(        
    strip.background = element_rect(fill=traits_cols["Ancestry"]),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", linewidth=1)
  )

#What if I keep the same donors across tissues
to_plot$SUBJID <- unlist(to_plot$SUBJID)
donors <- to_plot$SUBJID[to_plot$Tissue %in% c("Lung", "ColonTransverse", "Ovary")]
# donors <- to_plot$SUBJID[to_plot$Tissue %in% c("Lung", "ColonTransverse")]
shared_donors <- names(table(donors))[table(donors)>2]
to_plot_s <- to_plot[to_plot$SUBJID %in% shared_donors,]
ggplot(to_plot_s, aes(Ancestry, Betas)) +     
  geom_violin(aes(fill = Ancestry), col = "black", scale = "width") +
  geom_jitter(col = "black",
              alpha = 0.5,
              size = 0.8) +
  theme_bw() +
  facet_grid(~Tissue) +
  theme(        
    strip.background = element_rect(fill=traits_cols["Ancestry"]),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", linewidth=1)
  )

#Plot with SNP. In ovary we only had 10:

# [1] "chr14_74413354_G_A_b38"  "chr14_74510149_G_A_b38"  "chr14_74565195_T_G_b38"  "chr14_74588063_G_A_b38" 
# [5] "chr14_74617925_G_GC_b38" "chr14_74626565_C_A_b38"  "chr14_74643734_G_A_b38"  "chr14_74646793_A_G_b38" 
# [9] "chr14_74730275_C_G_b38"  "chr14_74772526_G_T_b38"
# 
# but the significant ones are:
# chr14_74730275_C_G_b381 #this is by far the one
# chr14_74588063_G_A_b381
# chr14_74617925_G_GC_b381
# chr14_74643734_G_A_b381 #0.1
# chr14_74646793_A_G_b381 #0.1

#read donor and variant data
genotypes <- read.table(paste0(first_path, "/projects/bsc83/Projects/GTEx_v8/Methylation/Data/Fst/imQTL_GT/Ovary.imQTLs.Donor_genotypes.txt.gz"))
genotypes <- read.table(paste0(first_path, "/scratch/bsc83/bsc83535/GTEx/v8/genotype_data/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf.gz"))
