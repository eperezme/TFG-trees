install.packages(c("ape", "ggplot2","BiocManager", 'svglite', 'dplyr', 'ggpubr'))
library(BiocManager)
BiocManager::install(c("ggtree", 'Biostrings' ))
library(ggplot2)
library(treeio)
library(ggtree)
library(svglite)
library(dplyr)
library(ggpubr)

# file <- "/Users/eduard/Google Drive/Uni/TFG/Seqs/Named/Trees/Relaxed_clock.tree"
# file <- "C:/Users/eduar/Google Drive/Uni/TFG/Seqs/Named/Trees/Relaxed_clock.tree"
file <- "C:/Users/eduar/Google Drive/Uni/TFG/Seqs/Named/Trees/annotated.tree" #Windows
file <- "/Users/eduard/Google Drive/Uni/TFG/Seqs/Named/Trees/annotated.tree" #Mac

beast_tree <- read.beast(file)
as.phylo(beast_tree)
get.fields(beast_tree)
get.data(beast_tree)

beast_tree@data$Lineage[is.na(beast_tree@data$Lineage)] <- 'NO'

tree0 <- ggtree(beast_tree)

d <- tree0$data
d <- d[!d$isTip,]
d$label <- round(as.numeric(d$posterior),2)
d <- d[d$label > 0.50,]

tree0 + geom_text(data=d, aes(label=label), hjust=-.2)
#### LINEAGE NODES  ####

du.node <- MRCA(beast_tree, c("DU_LT617553", "DU_LT617552")) #DUERO
da.node <- MRCA(beast_tree, c("DA_KF985732", "DA_FJ608998")) #DANUBIAN
ad.node <- MRCA(beast_tree, c("AD_LT617529", "AD_LT617522")) #ADRIATIC
ma.node <- MRCA(beast_tree, c("MA_LT617576", "MA_LT617574")) #MARMORATUS
ti.node <- MRCA(beast_tree, c("TU_LT617587", "TU_LT617585")) #TIGRIS
me.node <- MRCA(beast_tree, c("ME_LT617583", "ME_LT617579")) #MEDITERRANEAN
dd.node <- MRCA(beast_tree, c("DD_LT617545", "DD_LT617543")) #DADES
at.node <- MRCA(beast_tree, c("AT_KF985738", "AT_MN528759")) #ATLANTIC
na.node <- MRCA(beast_tree, c("NA_LT617562", "NA_LT617563")) #NORTH-AFRICA
ohr.node <- MRCA(beast_tree, c("Ohr_AF053590", "Ohr_JX960764")) #Ohridanus
sal.node <- MRCA(beast_tree, c("SAL_JX960834", "SAL_AF053591"))

tiplab <- beast_tree@phylo$tip.label
da.lab <- tiplab[grep("^DA", tiplab)]
ad.lab <- tiplab[grep("^AD", tiplab)]
ma.lab <- tiplab[grep("^MA", tiplab)]
ti.lab <- tiplab[grep("^TU", tiplab)]
me.lab <- tiplab[grep("^ME", tiplab)]
dd.lab <- tiplab[grep("^DD", tiplab)]
at.lab <- tiplab[grep("^AT", tiplab)]
na.lab <- tiplab[grep("^NA", tiplab)]
ohr.lab <- tiplab[grep("^Ohr", tiplab)]
sal.lab <- tiplab[grep("^SAL", tiplab)]




#### TREEE ####

tree <- ggtree(beast_tree, aes(color = Lineage)) + theme_tree() + 
  geom_tiplab(size=1.5, color="black") +
  geom_rootedge(.001) +
  geom_cladelabel(node = 344, label = 'Duero Lineage', color = '#FF6362', align = T, offset.text = .001, offset = 0.004, barsize = 1.5) +
  geom_cladelabel(node = da.node, label = 'Danubian Lineage', color = '#F5B400', align = T, offset.text = .001, offset = 0.004, barsize = 1.5) +
  geom_cladelabel(node = ad.node, label = 'Adriatic Lineage', color = '#DC4437', align = T, offset.text = .001, offset = 0.004, barsize = 1.5) +
  geom_cladelabel(node = ma.node, label = 'Marmoratus Lineage', color = '#58508D', align = T, offset.text = .001, offset = 0.004, barsize = 1.5) +
  geom_cladelabel(node = ti.node, label = 'Tigris Lineage', color = '#86D112', align = T, offset.text = .001, offset = 0.004, barsize = 1.5) +
  geom_cladelabel(node = me.node, label = 'Mediterranean Lineage', color = '#C67ECA', align = T, offset.text = .001, offset = 0.004, barsize = 1.5) +
  geom_cladelabel(node = dd.node, label = 'Dades Lineage', color = '#4385F5', align = T, offset.text = .001, offset = 0.004, barsize = 1.5) +
  geom_cladelabel(node = at.node, label = 'Atlantic Lineage', color = '#109D59', align = T, offset.text = .001, offset = 0.004, barsize = 1.5) +
  geom_cladelabel(node = na.node, label = 'North-African Lineage', color = '#0154A4', align = T, offset.text = .001, offset = 0.004, barsize = 1.5) +
  geom_text2(data=d, aes(label=label, subset = node %in% c(183)), hjust = -0.3, color = 'black', size = 3) + #DA+AD 0.09
  geom_text2(data=d, aes(label=label, subset = node %in% c(177, 178, 179, 180, 181, 182, 184, 195, 196, 197, 204, 255, 267, 268, 292, 298)), 
             hjust = 1.15, vjust = -0.5, color = 'black', size = 3) +
  geom_text2(data=d, aes(label=label, subset = node %in% c(248)), hjust = -0.15, vjust = 0.2, color = 'black', size = 3) +
  scale_colour_manual(values = c('#DC4437', '#109D59', '#F5B400', '#4385F5', '#FF8483', '#58508D',
                                 '#C67ECA', '#0154A4', '#3B3B3B' ,'#86D112'), 
                      breaks = c("AD", 'AT', 'DA', 'DD','DU', 'MA', 'ME', 'NA',
                                 'TI') ) +
  theme(legend.title = element_text(colour = 'black', size = 10, face='bold'), legend.position = 'bottom', 
      legend.box = 'horizontal') +
  xlim_tree(0.04)
tree
#geom_highlight(node = 344, fill = '#FF6362') +
  #geom_highlight(node = da.node, fill = '#F5B400') +
  #geom_highlight(node = ad.node, fill = '#DC4437') +
   #geom_nodelab(aes(x=branch, label=posterior, subset = !is.na(as.numeric(posterior)) & as.numeric(posterior) > 80), vjust=-.5, size = 3)
#geom_text2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 80))

 # geom_highlight(node=da.node, fill = 'gold')


tree3 <- ggtree(beast_tree, aes(color = Lineage)) + theme_tree() + 
  geom_tiplab(size=1.5, color="black") +
  scale_colour_manual(values = c('#DC4437', '#109D59', '#F5B400', '#4385F5', '#FF8483', '#58508D',
                                 '#C67ECA', '#0154A4', '#3B3B3B' ,'#86D112'), 
                      breaks = c("AD", 'AT', 'DA', 'DD','DU', 'MA', 'ME', 'NA',
                                 'TI') ) +
  theme(legend.title = element_text(colour = 'black', size = 10, face='bold'), legend.position = 'bottom', 
        legend.box = 'horizontal') +
  xlim_tree(0.035)


## COLLAPSE TREE ##
tree3 <- tree3 %>%
  collapse(node = da.node) %>%
  collapse(node = ad.node) %>%
  collapse(node = ma.node) %>%
  collapse(node = ti.node) %>%
  collapse(node = me.node) %>%
  collapse(node = dd.node) %>%
  collapse(node = at.node) %>%
  collapse(node = na.node) %>%
  collapse(node = du.node) %>%
  collapse(node = ohr.node) %>%
  collapse(node = sal.node)

## ADD CLADES TO TREE ##
tree3 <- tree3 + 
  geom_cladelabel(node = da.node, label = 'Danubian Lineage', color = '#F5B400', align = T, offset.text = 0, offset = -.0025, barsize = 1.5) +
  geom_cladelabel(node = ad.node, label = 'Adriatic Lineage', color = '#DC4437', align = T, offset.text = 0, offset = -.002, barsize = 1.5) +
  geom_cladelabel(node = ma.node, label = 'Marmoratus Lineage', color = '#58508D', align = T, offset.text = 0, offset = -.001, barsize = 1.5) +
  geom_cladelabel(node = ti.node, label = 'Tigris Lineage', color = '#86D112', align = T, offset.text = 0, offset = -.002, barsize = 1.5) +
  geom_cladelabel(node = me.node, label = 'Mediterranean Lineage', color = '#C67ECA', align = T, offset.text = 0, offset = -.001, barsize = 1.5) +
  geom_cladelabel(node = dd.node, label = 'Dades Lineage', color = '#4385F5', align = T, offset.text = 0, offset = -.0008, barsize = 1.5) +
  geom_cladelabel(node = at.node, label = 'Atlantic Lineage', color = '#109D59', align = T, offset.text = 0, offset = -.002, barsize = 1.5) +
  geom_cladelabel(node = na.node, label = 'North-African Lineage', color = '#0154A4', align = T, offset.text = 0, offset = -.0015, barsize = 1.5) +
  geom_cladelabel(node = du.node, label = 'Duero Lineage', color = '#FF6362', align = T, offset.text = 0, offset = -.0005, barsize = 1.5) +
  geom_cladelabel(node = ohr.node, label = 'S.ohridanus', color = '#3B3B3B', align = T, offset.text = 0, offset = -.001, barsize = 1.5) +
  geom_cladelabel(node = sal.node, label = 'S.salar', color = '#3B3B3B', align = T, offset.text = 0, offset = -.003, barsize = 1.5) 
  
tree3$data$posterior <- round(as.numeric(tree3$data$posterior), 2) # ROUND POSTERIOR TO 2 DIGITS

 ## ADD POSTERIOR VALUES ##
bi_tree_ant <- tree3 +
  #geom_text(aes(label = node))
  geom_text2(aes(label=posterior, subset = node %in% c(177, 178, 179, 180, 181, 196, 197, 267)), hjust = 1.15, vjust = -0.5, color = 'black', size = 3) +
  geom_text2(aes(label=posterior, subset = node %in% c(182)), hjust = 1.15, vjust = 1.5, color = 'black', size = 3) +
  geom_text2(aes(label=posterior, subset = node %in% c(195)), hjust = -.15, color = 'black', size = 3)


bi_tree_ant

ggsave(filename = "BI_collapsed_posterior.svg") # SAVE TREE #

####  CIRCULAR TREE #####



to_drop <- c("SAL_AF053591", "SAL_EU492280", "SAL_JX960834")

beast_tree2 <- drop.tip(beast_tree, to_drop)

beast_tree2 



## GENERATE TREE ##
tree2 <- ggtree(beast_tree2, aes(color = Lineage), layout = 'circular') +
  geom_rootedge(0.002) +
  geom_tiplab(size=2, color="black") +
  geom_cladelabel(node = 340, label = 'Duero Lineage', color = '#FF6362',
                  align = T, hjust = 'centre', offset.text = 0.001, offset = 0.0025, barsize = 1.5,
                  angle = -69, family = 'Times New Roman') +
  geom_cladelabel(node = 211, label = 'Danubian Lineage', color = '#F5B400', 
                  align = T, hjust = 'centre', offset.text = 0.001, offset = 0.0025, barsize = 1.5,
                  angle = 30, family = 'Times New Roman') +
  geom_cladelabel(node = 179, label = 'Adriatic Lineage', color = '#DC4437', 
                  align = T, hjust = 'centre', offset.text = 0.001, offset = 0.0025, barsize = 1.5,
                  angle = -75, family = 'Times New Roman') +
  geom_cladelabel(node = 205, label = 'Marmoratus Lineage', color = '#58508D', 
                  align = T, hjust = 'centre', offset.text = 0.002, offset = 0.0025, barsize = 1.5,
                  angle = -48, family = 'Times New Roman') +
  geom_cladelabel(node = 209, label = 'Tigris Lineage', color = '#86D112', 
                  align = T, hjust = 'centre', offset.text = 0.001, offset = 0.0025, barsize = 1.5,
                  angle = -57, family = 'Times New Roman') +
  geom_cladelabel(node = 198, label = 'Mediterranean Lineage', color = '#C67ECA', 
                  align = T, hjust = 'centre', offset.text = 0.002, offset = 0.0025, barsize = 1.5,
                  angle = -27, family = 'Times New Roman') +
  geom_cladelabel(node = 194, label = 'Dades Lineage', color = '#4385F5', 
                  align = T, hjust = 'centre', offset.text = 0.001, offset = 0.0025, barsize = 1.5,
                  angle = -40, family = 'Times New Roman') +
  geom_cladelabel(node = 264, label = 'Atlantic Lineage', color = '#109D59', 
                  align = T, hjust = 'centre', offset.text = 0.001, offset = 0.0025, barsize = 1.5,
                  angle = 30, family = 'Times New Roman') +
  geom_cladelabel(node = 323, label = 'North-African Lineage', color = '#0154A4', 
                  align = T, hjust = 'centre', offset.text = 0.001, offset = 0.0025, barsize = 1.5,
                  angle = -47, family =  'Times New Roman') 
  geom_text2(aes(label = posterior), color = "black", size = 3.5)
tree2

tree2$data$posterior <- round(as.numeric(tree2$data$posterior), 2) # ROUND POSTERIOR #
tree2$data$posterior <- tree2$data$posterior*100 # MAKE POSTERIOR AS PERCENTAGE #
# 
# tree2_nodes <- tree2 + geom_text2(aes(label = node)) ## Tree with node numbers ##


## ADD POSTERIOR VALUES ##
tree5 <- tree2 +
  # #geom_text2(aes(label = node))
  # geom_text2(aes(label = posterior, subset = node %in% 
  #                          c(177)), hjust = 1.1,vjust= 1.5, color = 'black', size = 2.5, angle = 25) +
  # geom_text2(aes(label = posterior, subset = node %in% 
  #                          c(178)), hjust = 1.2, vjust = -.2, color = 'black', size = 2.5, angle = 46) +
  # geom_text2(aes(label = posterior, subset = node %in% 
  #                          c(179)), hjust = 1.2, vjust = -.1, color = 'black', size = 2.5, angle = 80) +
  # geom_text2(aes(label = posterior, subset = node %in% 
  #                          c(180)), hjust = -.1, vjust = 1.2, color = 'black', size = 2.5, angle = -25) +
  # geom_text2(aes(label = posterior, subset = node %in% 
  #                          c(181)), hjust = -.3, vjust = 1.3, color = 'black', size = 2.5, angle = 58) +
  # geom_text2(aes(label = posterior, subset = node %in% 
  #                          c(182)), vjust = 1.3, color = 'black', size = 2.5, angle = -65) +
  # #geom_text2(aes(label = posterior, subset = node %in% 
  # #                         c(195)), vjust = 1.2, color = 'black', size = 2.5, angle = -40) +
  # geom_text2(aes(label = posterior, subset = node %in% 
  #                          c(196)), vjust = 1.3, color = 'black', size = 2.5, angle = -40) +
  # geom_text2(aes(label = posterior, subset = node %in% 
  #                          c(197)), vjust = 1.3, hjust = -.3, color = 'black', size = 2.5, angle = 58) +
  # geom_text2(aes(label = posterior, subset = node %in% 
  #                          c(267)), vjust = 1.3, hjust = 1.1, color = 'black', size = 2.5, angle = 66) +
  
  scale_colour_manual(values = c('#DC4437', '#109D59', '#F5B400', '#4385F5', '#FF8483', '#58508D',
                                 '#C67ECA', '#0154A4', '#3B3B3B' ,'#86D112'), 
                      breaks = c("AD", 'AT', 'DA', 'DS','DU', 'MA', 'ME', 'NA',
                                 'TI') ) +
  theme(legend.title = element_text(colour = 'black', size = 10, face='bold'), legend.position = 'bottom', 
        legend.box = 'horizontal')

tree5


## SAVE SVG ##
ggsave(filename = "Circular_to_be_edited.png", width = 10, height = 10, dpi = 320)


ggarrange(tree2_nodes,tree5, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

geom_text2(aes(label = posterior, subset = node %in% 
                 c(268)), vjust = -.4, color = 'black', size = 3, angle = 10) +
  
  geom_text2(aes(label = posterior, subset = node %in% c(248)), vjust = 1.3, 
             color = 'black', size = 3, angle = -59) +
  geom_text2(aes(label = posterior, subset = node %in% 
                   c(183)), vjust = 1.3, color = 'black', size = 3, angle = -17) + #DA+AD 0.09
  
  
