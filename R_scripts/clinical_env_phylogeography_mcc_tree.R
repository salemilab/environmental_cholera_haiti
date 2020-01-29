setwd('/Users/tpaisie/Dropbox (UFL)/cholera/may_2019/env_paper/ouest_only/mcc_location_tree/')

library(ggtree)
library(ape)
library(ggplot2)
library(colorspace)
library(Biostrings)
library(phytools)
library(treeio)
library(dplyr)
library(readr)
library(ggnewscale)

#molecular clock clinical-environmental mcc tree
mcc <- read.beast(file = "OUEST-ONLY-snps_RC_BSP_ASYM_BSVSS_PSSS_MCC.tree")
loc <- read.table(file = "ouest_only_zone_location.csv", sep = ',', header = TRUE, row.names = 1)
tree <- read.tree(file = 'OUEST-ONLY-snps_RC_BSP_ASYM_BSVSS_PSSS_MCC.nwk')

get.fields(mcc)


p <- ggtree(mcc, aes(color = location), mrsd = "2015-12-05", size = 0.75) + scale_color_manual(values=c("mediumpurple4", "springgreen3")) +
  geom_point2(aes(label = posterior, subset = posterior >= 0.9), size = 2.5) +
  scale_x_continuous(breaks = seq(2010, 2017)) + theme_tree2() + geom_tiplab(size=2.5, color="black")

p <- ggtree(mcc, mrsd = "2015-12-05", size = 1, ladderize = TRUE) + geom_tiplab(size=2.5, color="black")
plot(p)

pp <- (p + scale_y_continuous(expand=c(0, 0.3))) %>% gheatmap(loc, offset=0, width = 0.05, colnames = FALSE) + 
  scale_color_manual(breaks = c("Ouest", "D", "Sud", "Artibonite", "B", "C", "D", "A", "E", "F")) + theme_tree2(legend.position='right')
plot(pp)