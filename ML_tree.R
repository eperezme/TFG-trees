install.packages(c("ape", "ggplot2","BiocManager", 'svglite', 'dplyr'))
library(BiocManager)
BiocManager::install(c("ggtree", 'Biostrings' ))
library(ggplot2)
library(treeio)
library(ggtree)
library(dplyr)
library(svglite)
library(cowplot)


file_ml <- "C:/Users/eduar/Google Drive/Uni/TFG/Seqs/Named/Aaaaaaaaa/Trees/ML_Tree.trees" #Windows
file_ml <- "/Users/eduard/Google Drive/Uni/TFG/Seqs/Named/Aaaaaaaaa/Trees/ML_Tree.trees" #Mac

ml_trees <- read.beast(file_ml)
as.phylo(ml_trees)
get.fields(ml_trees)
ml_trees@data$Lineage[is.na(ml_trees@data$Lineage)] <- 'NO'
ml_trees@phylo <-  root(ml_trees@phylo, node = 337)
# "SAL_AF053591", "SAL_EU492280", "SAL_JX960834"
drop_sal <- c("SAL_EU492280", "SAL_JX960834")
drop_DD <- c('DD_LT617541', 'DD_LT617544', 'DD_LT617542', 'DD_LT617543')
ml_trees <- drop.tip(ml_trees, drop_sal)
ml_trees <- drop.tip(ml_trees, drop_DD)
tree_ml <- ggtree(ml_trees) + geom_tiplab()
tree_ml

ti_ml_lenght <- 0.00120 + 0.00120 + 0.000002 + 0.00120 + 0.001203

ml_tree_collapsed <- ml_trees
ml_tree_collapsed@phylo$edge.length[223] <- ti_ml_lenght #Change OTU name
ml_tree_collapsed@data$Lineage[ml_tree_collapsed@data$node == 278] <- 'TI' #Change OTU name

tree_ml <- ggtree(ml_tree_collapsed, aes(color = Lineage))
tree_ml <- tree_ml %>%
  collapse(node =du.node.ml) %>%
  collapse(node =227) %>%
  collapse(node =ad.node.ml) %>%
  collapse(node =ma.node.ml) %>%
  collapse(node =na.node.ml) %>%
  # collapse(node =279) %>%
  collapse(node =326) %>%
  collapse(node =263) %>%
  collapse(node =278) %>%
  collapse(node = 172)



tree_ml$data$Lineage[is.na(tree_ml$data$Lineage)] <- 'NO'

tree3 <- tree_ml + theme_tree() + 
  scale_colour_manual(values = c('#DC4437', '#109D59', '#F5B400', '#4385F5', '#FF8483', '#58508D',
                                 '#C67ECA', '#0154A4', '#3B3B3B' ,'#ff297f'), 
                      breaks = c("AD", 'AT', 'DA', 'DS','DU', 'MA', 'ME', 'NA',
                                 'TI') ) +
  theme(legend.title = element_text(colour = 'black', size = 10, face='bold'), legend.position = 'bottom', 
        legend.box = 'horizontal')


tree3 <- tree3 +
  geom_cladelabel(node = 322, label = 'Duero Lineage', color = '#FF6362', align = T, offset.text = .001, offset = 0.004, barsize = 1.5) +
  geom_cladelabel(node = 172, label = 'Danubian Lineage', color = '#F5B400', align = T, offset.text = .001, offset = 0.004, barsize = 1.5) +
  geom_cladelabel(node = 234, label = 'Adriatic Lineage', color = '#DC4437', align = T, offset.text = .001, offset = 0.004, barsize = 1.5) +
  geom_cladelabel(node = 222, label = 'Marmoratus Lineage', color = '#58508D', align = T, offset.text = .001, offset = 0.004, barsize = 1.5) +
  geom_cladelabel(node = 278, label = 'Tigris Lineage', color = '#86D112', align = T, offset.text = .001, offset = 0.004, barsize = 1.5) +
  geom_cladelabel(node = 227, label = 'Mediterranean Lineage', color = '#C67ECA', align = T, offset.text = .001, offset = 0.004, barsize = 1.5) +
  geom_cladelabel(node = 57, label = 'Dades Lineage', color = '#4385F5', align = T, offset.text = .001, offset = 0.004, barsize = 1.5) +
  geom_cladelabel(node = 263, label = 'Atlantic Lineage', color = '#109D59', align = T, offset.text = .001, offset = 0.004, barsize = 1.5) +
  geom_cladelabel(node = 245, label = 'North-African Lineage', color = '#0154A4', align = T, offset.text = .001, offset = 0.004, barsize = 1.5) +
  geom_cladelabel(node = 326, label = expression(italic('S.ohridanus')), color = '#3B3B3B', align = T, offset.text = 0, offset = -.001, barsize = 1.5) +
  geom_cladelabel(node = 164, label = expression(italic('S.salar')), color = '#3B3B3B', align = T, offset.text = 0, offset = -.002, barsize = 1.5) +
  xlim(-.001,.1)

tree3

ggsave(filename = "ML_collapsed.svg")

ml_tree_ant <- tree3 +
 # geom_text2(aes(label = node))
  geom_text2(aes(label = bootstrap, subset = node %in% c(234, 263)), color = 'black', hjust = -.3) +
  geom_text2(aes(label = bootstrap, subset = node %in% c(176, 177, 178, 179, 180)), color = 'black', hjust = 1.1, vjust = -.2)



### SAVE FILES ###
ml_tree_ant
ggsave(filename = "ML_collapsed_bootsrap.svg")

gridplot <- plot_grid(ml_tree_ant, bi_tree_ant, labels = "AUTO")
ggsave(filename = "ML&BI_Trees.svg", plot = gridplot, width = 18, height = 8, dpi = 600)

ggsave(filename = "ML&BI_Trees.jpeg", plot = gridplot, width = 18, height = 8, dpi = 600)

ggsave(filename = "ML&BI_Trees.svg", plot = gridplot, width = 18, height = 8, dpi = 600)









geom_text2(data=d, aes(label=label, subset = node %in% c(183)), hjust = -0.3, color = 'black', size = 3) + #DA+AD 0.09
  geom_text2(data=d, aes(label=label, subset = node %in% c(177, 178, 179, 180, 181, 182, 184, 195, 196, 197, 204, 255, 267, 268, 292, 298)), hjust = 1.15, vjust = -0.5, color = 'black', size = 3) +
  geom_text2(data=d, aes(label=label, subset = node %in% c(248)), hjust = -0.15, vjust = 0.2, color = 'black', size = 3) +
  


ggdensitree(ml_trees, alpha=.3, colour='steelblue') + 
  geom_tiplab(size=3) + xlim(0, 45)


