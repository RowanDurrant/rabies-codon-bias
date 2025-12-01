#Phylogenetic tree
library(ggtree)
library(treeio)
library(ggplot2)
library(ggbreak)
library(grid)
library(stringr)

tree = read.tree("sequence_data/trees/all_seqs.fasta.treefile")

for(i in 1:length(tree$node.label)){
  if(i %in% c(2,3,4,113,114,299,300,301,302)){
    tree$node.label[i] = strsplit(tree$node.label[i], "/")[[1]][2]
  }
  else{tree$node.label[i] = NA}
}

tipcolours = c()

df = read.csv("sequence_data/metadata.csv")
df = rbind(df, c("gannoruwa_outgroup", NA, "Gannoruwa bat lyssavirus", "bat", "bat"))

for(i in 1:length(tree$tip.label)){
  tipcolours[i] = df$Clade[df$Accession == tree$tip.label[i]]
}
d <- data.frame(node=c(1:(Nnode(tree)+length(tree$tip.label))), 
                clade = c(tipcolours, rep("black", Nnode(tree))))
d$clade = factor(d$clade, c("Cosmo AF1b", "Cosmo AM2a", "Arctic A", "Asian SEA2a", 
                            "Asian SEA2b", 
                            "Bat DR","Bat TB1", "Bat LC",
                            "Bat EF-E2","RAC-SK SCSK",
                          "Gannoruwa bat lyssavirus"
))

p2 = ggtree(tree) +
  ylim(0, 440) +
  geom_tippoint(size = 0.25, colour = "black")+
  geom_cladelabel(node=431, label="NC_031988", 
                  color='black', fontsize=3)+
  geom_cladelabel(node=706, label="Cosmopolitan AM2a", 
                  color='#88CCEE', fontsize=4, 
                  offset = 0.005, offset.text = 0.005)+
  geom_cladelabel(node=714, label="Arctic A", 
                  color='#CCDDAA', fontsize=4, 
                  offset = 0.005, offset.text = 0.005)+
  geom_cladelabel(node=436, label="Asian SEA2a", 
                  color='#44AA99', fontsize=4, 
                  offset = 0.005, offset.text = 0.005)+
  geom_cladelabel(node=517, label="Asian SEA2a", 
                  color='#44AA99', fontsize=4, 
                  offset = 0.005, offset.text = 0.005)+
  geom_cladelabel(node=483, label="Asian SEA2b", 
                  color='#117733', fontsize=4, 
                  offset = 0.005, offset.text = 0.005)+
  geom_cladelabel(node=734, label="Bat DR", 
                  color='#DDCC77', fontsize=4, 
                  offset = 0.005, offset.text = 0.005)+
  geom_cladelabel(node=772, label="Bat TB1", 
                  color='#999933', fontsize=4, 
                  offset = 0.005, offset.text = 0.005)+
  geom_cladelabel(node=803, label="Bat LC", 
                  color='#AA4499', fontsize=4, 
                  offset = 0.005, offset.text = 0.005)+
  geom_cladelabel(node=818, label="Bat EF-E2", 
                  color='#CC6677', fontsize=4, 
                  offset = 0.005, offset.text = 0.005)+
  geom_cladelabel(node=853, label="RAC-SK SCSK", 
                  color='#882255', fontsize=4, 
                  offset = 0.005, offset.text = 0.005)+
  geom_cladelabel(node=546, label="Cosmopolitan\nAF1b", 
                  color='#332288', fontsize=4, 
                  offset = 0.005, offset.text = 0.005)+
  xlim(0, 0.465)+ geom_treescale(x = 0, y = 320, fontsize = 4, 
                                offset = 7)+ 
  geom_nodelab(hjust = 0, fontsize = 2)
p2

source("code/ENC.R")

library(ggmap)
p3 = p2 +
  inset(ggplotGrob(p), xmin = 0.24, xmax = 0.49, ymin = 20, ymax = 430)
p3


png("plots/Figure 1.png", width = 10, height = 7, units = 'in', res = 600)
p3
dev.off()
