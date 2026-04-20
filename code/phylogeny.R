#Phylogenetic tree
library(ggtree)
library(treeio)
library(ggplot2)
library(ggbreak)
library(grid)
library(stringr)

tree = read.tree("sequence_data/all_outgroups/all_outgroups_aligned.fasta.treefile")

# df = read.csv("sequence_data/metadata.csv")
# df = rbind(df, c("Gannoruwa_bat_lyssavirus_NC_031988.1", "Other lyssavirus", NA,  "bat"))
# df = rbind(df, c("NC_003243.1_ABLV", "Other lyssavirus", NA,  "bat"))
# df = rbind(df, c("NC_009527.1_eblv1", "Other lyssavirus", NA,  "bat"))
# 
# x = full_join(tree, df, by = c("label" = "Accession"))
# 
# ggtree(x)+
#   geom_nodelab()+
#   geom_tippoint(aes(color=Clade))

for(i in 1:length(tree$node.label)){
  if(i %in% c(1,2,3,4,73,197,198,295,294,293,284)){
    tree$node.label[i] = strsplit(tree$node.label[i], "/")[[1]][2]
  }
  else{tree$node.label[i] = NA}
}

p2 = ggtree(tree) +
  geom_tippoint(size = 0.25, colour = "black")+
  geom_cladelabel(node=(416+82), label="Cosmo\nAF1b", 
                  color='#332288', fontsize=3, 
                  offset = 0.005, offset.text = 0.003)+
  geom_cladelabel(node=(416+74), label="Cosmo AM2a", 
                  color='#88CCEE', fontsize=3, 
                  offset = 0.005, offset.text = 0.003)+
  geom_cladelabel(node=(416+5), label="Arctic A", 
                  color='#CCDDAA', fontsize=3, 
                  offset = 0.005, offset.text = 0.003)+
  geom_cladelabel(node=(416+199), label="Asian\nSEA2a", 
                  color='#44AA99', fontsize=3, 
                  offset = 0.005, offset.text = 0.003)+
  geom_cladelabel(node=(416+274), label="Asian SEA2a", 
                  color='#44AA99', fontsize=3, 
                  offset = 0.005, offset.text = 0.003)+
  geom_cladelabel(node=(416+246), label="Asian\nSEA2b", 
                  color='#117733', fontsize=3, 
                  offset = 0.005, offset.text = 0.003)+
  geom_cladelabel(node=(416+296), label="Bat DR", 
                  color='#DDCC77', fontsize=3, 
                  offset = 0.005, offset.text = 0.003)+
  geom_cladelabel(node=(416+335), label="Bat TB1", 
                  color='#999933', fontsize=3, 
                  offset = 0.005, offset.text = 0.003)+
  geom_cladelabel(node=(416+365), label="Bat LC", 
                  color='#AA4499', fontsize=3, 
                  offset = 0.005, offset.text = 0.003)+
  geom_cladelabel(node=(416+380), label="Bat EF-E2", 
                  color='#CC6677', fontsize=3, 
                  offset = 0.005, offset.text = 0.003)+
  geom_cladelabel(node=(416+285), label="RAC-SK\nSCSK", 
                  color='#882255', fontsize=3, 
                  offset = 0.005, offset.text = 0.003)+ 
  geom_nodelab(hjust = 0, size = 2)+
  xlim(0, 0.52)+ 
  geom_treescale(x = 0.1, y = 50, fontsize = 4, 
                                offset = 7)
p2

source("code/ENC.R")

library(ggmap)
p3 = p2 +
  inset(ggplotGrob(p), xmin = -0.02, xmax = 0.25, ymin = 120, ymax = 410)
p3

png("plots/Figure 1.png", width = 10, height = 7, units = 'in', res = 600)
p3
dev.off()
