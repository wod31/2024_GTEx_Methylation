#### main plot overview varPart ####
### Plot varPart #####
library(variancePartition)
first_dir <- "~/"
setwd(paste0(first_dir, "marenostrum/Projects/GTEx_v8/Methylation/"))

# -------------- #
print(Sys.time())
#-------------- #

tissues <- list.dirs("Tissues/", full.names = F)[-1]

#### reading varPart values ####
beta <- readRDS(paste0("var_part_expression.rds"))
beta$gene <- gsub('\\..*','',rownames(beta))

library(reshape2)
betas_locations_m <- melt(data = beta[,c("gene", "Donor", "Tissue")], id.vars = c('gene'),
                          variable.name = 'Variable', value.name = 'Proportion')
head(betas_locations_m)

library("dplyr")
mu <- betas_locations_m %>% 
  group_by(Variable) %>%
  summarise(grp.mean = mean(Proportion))
mu

pdf(file = paste0('Plots/',"/VarPart_all_main2_expression.pdf"), w = 4, h = 3)
ggplot(betas_locations_m, aes(x = Proportion, fill=Variable)) + 
  geom_density() +
  scale_fill_brewer(palette="Dark2") +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Variable),
             linetype="dashed") +
  geom_text(
    size    = 5,
    data    = mu,
    mapping = aes(x = Inf, y = Inf, label = round(mu$grp.mean, digits = 3)),
    hjust   = 1.05,
    vjust   = 1.5
  )+
  scale_color_brewer(palette="Dark2") +
  theme_bw()
dev.off()
