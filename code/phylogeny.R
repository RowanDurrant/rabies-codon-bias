#Phylogenetic tree
library(ggtree)
library(treeio)
library(ggplot2)
library(ggbreak)
library(grid)
library(stringr)

tree = read.tree("sequence_data/tree/all_seqs.fasta.treefile")

for(i in 1:length(tree$node.label)){
  if(i %in% c(1,2,3,4,73, 197,198,284,285,286,287)){
    tree$node.label[i] = strsplit(tree$node.label[i], "/")[[1]][2]
  }
  else{tree$node.label[i] = NA}
}

tipcolours = c()

df = read.csv("sequence_data/metadata.csv")
df = rbind(df, c("Gannoruwa_bat_lyssavirus_NC_031988.1", "Gannoruwa bat lyssavirus", NA,  "bat"))

for(i in 1:length(tree$tip.label)){
  tipcolours[i] = df$Clade[df$Accession == tree$tip.label[i]]
}
d <- data.frame(node=c(1:(Nnode(tree)+length(tree$tip.label))), 
                clade = c(tipcolours, rep("black", Nnode(tree))))
d$clade = factor(d$clade, c("Cosmopolitan AF1b", "Cosmopolitan AM2a", 
                            "Arctic A", "Asian SEA2a", 
                            "Asian SEA2b", 
                            "Bats DR","Bats TB1", "Bats LC",
                            "Bats EF-E2","RAC-SK SCSK",
                          "Gannoruwa bat lyssavirus"
))

p2 = ggtree(tree) +
  ylim(0, 415) +
  geom_tippoint(size = 0.25, colour = "black")+
  geom_cladelabel(node=414, label="NC_031988", 
                  color='black', fontsize=3)+
  geom_cladelabel(node=496, label="Cosmo\nAF1b", 
                  color='#332288', fontsize=3, 
                  offset = 0.005, offset.text = 0.003)+
  geom_cladelabel(node=488, label="Cosmo AM2a", 
                  color='#88CCEE', fontsize=3, 
                  offset = 0.005, offset.text = 0.003)+
  geom_cladelabel(node=419, label="Arctic A", 
                  color='#CCDDAA', fontsize=3, 
                  offset = 0.005, offset.text = 0.003)+
  geom_cladelabel(node=613, label="Asian\nSEA2a", 
                  color='#44AA99', fontsize=3, 
                  offset = 0.005, offset.text = 0.003)+
  geom_cladelabel(node=272, label="Asian SEA2a", 
                  color='#44AA99', fontsize=3, 
                  offset = 0.005, offset.text = 0.003)+
  geom_cladelabel(node=660, label="Asian\nSEA2b", 
                  color='#117733', fontsize=3, 
                  offset = 0.005, offset.text = 0.003)+
  geom_cladelabel(node=702, label="Bat DR", 
                  color='#DDCC77', fontsize=3, 
                  offset = 0.005, offset.text = 0.003)+
  geom_cladelabel(node=740, label="Bat TB1", 
                  color='#999933', fontsize=3, 
                  offset = 0.005, offset.text = 0.003)+
  geom_cladelabel(node=770, label="Bat LC", 
                  color='#AA4499', fontsize=3, 
                  offset = 0.005, offset.text = 0.003)+
  geom_cladelabel(node=785, label="Bat EF-E2", 
                  color='#CC6677', fontsize=3, 
                  offset = 0.005, offset.text = 0.003)+
  geom_cladelabel(node=819, label="RAC-SK SCSK", 
                  color='#882255', fontsize=3, 
                  offset = 0.005, offset.text = 0.003)+
  xlim(0, 0.465)+ geom_treescale(x = 0, y = 320, fontsize = 4, 
                                offset = 7)+ 
  geom_nodelab(hjust = 0)
p2

source("code/ENC.R")

library(ggmap)
p3 = p2 +
  inset(ggplotGrob(p), xmin = 0.255, xmax = 0.49, ymin = 20, ymax = 430)
p3


png("plots/Figure 1.png", width = 10, height = 7, units = 'in', res = 600)
p3
dev.off()
