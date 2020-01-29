setwd('/Users/tpaisie/Dropbox (UFL)/cholera/may_2019/env_paper/ouest_only/ml_tree/ouest_all/')

library(ggtree)
library(ape)
library(ggplot2)
library(colorspace)
library(Biostrings)
library(phytools)
library(treeio)
library(dplyr)
library(readr)
library(phangorn)
library(geiger)

source <- read.csv(file = "ouest_only_source_053019.csv", header = TRUE, sep = ",")
ml_tree <- read.tree('OUEST-ONLY-snps.fa.tree')

q <- ggtree(midpoint(ml_tree), ladderize = TRUE, layout = "circular")
d <- q$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 90,]
p <- q

dd <- data.frame(source, check.rows = TRUE, check.names = TRUE)
row.names(dd) <- NULL

cb_palette <- c('C' = 'mediumpurple4', 'E' = 'springgreen3')

p <- q %<+% dd + geom_tippoint(aes(color=source), size = 3.5) +
  geom_nodepoint(data=d, aes(label=label), size = 1.5, shape = 19) + 
  scale_color_manual(values = cb_palette) +
  geom_strip(54, 66, label = '2014-2015', align = T, color = 'pink', barsize = 3) +
  geom_strip(83, 53, align = T, label = '2013', color = '#3FA9F5', barsize = 3) +
  geom_strip(111, 91, label = '2010-2012', align = T, color = 'darkgrey', barsize = 3) + 
  geom_strip(34, 87, label = '2012-2013', align = T, color = 'gold', barsize = 3) + 
  geom_hilight(node = 123, fill = "pink", alpha = 0.2) +
  geom_hilight(node = 169, fill = 'darkgray', alpha = 0.2) + 
  geom_hilight(node = 150, fill = '#3FA9F5', alpha = 0.2) + 
  geom_hilight(node = 204, fill = 'gold', alpha = 0.2) + geom_treescale()

plot(p)  