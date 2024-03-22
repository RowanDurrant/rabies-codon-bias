#Phylogenetic tree
library(ggtree)
library(treeio)
library(ggplot2)
tree = read.tree("sequence_data/trees/rooted_N/all_seqs.fasta.treefile")
tipcolours = c()

df = read.csv("sequence_data/metadata.csv")
df = rbind(df, c("gannoruwa_outgroup", NA, "Gannoruwa bat lyssavirus", "bat", "bat"))

for(i in 1:length(tree$tip.label)){
  tipcolours[i] = df$Clade[df$Accession == tree$tip.label[i]]
}
d <- data.frame(node=c(1:(Nnode(tree)+length(tree$tip.label))), 
                clade = c(tipcolours, rep("black", Nnode(tree))))
d$clade = factor(d$clade, c("Cosmo AF1b", "Cosmo AM2a", "Asian SEA2a", 
                            "Asian SEA2b", 
                          "Bat TB1",
                          "Bat DR", "Bat EF-E2","RAC-SK SCSK", "Bat LC",
                          "Gannoruwa bat lyssavirus"
))
p = ggtree(tree) +
  theme_tree2()
p = p %<+% d  +  
   theme(legend.position="bottom") +
  ylim(0, 420) +
  geom_tippoint(aes(color=clade), size = 1) +
  scale_color_manual(values = c("#332288","#88CCEE","#44AA99","#117733","#999933",
                                "#DDCC77","#CC6677","#882255","#AA4499", "black"), 
                                name = "Clade", guide = guide_legend(),
                     labels = c("Cosmo AF1b\n(dog)",
                                "Cosmo AM2a\n(mongoose)",
                                "Asian SEA2a\n(dog)",
                                "Asian SEA2b\n(CFB)",
                                "Bat TB1\n(Mexican free-tailed bat)",
                                "Bat DR\n(vampire bat)",
                                "Bat EF-E2\n(big brown bat)",
                                "RAC-SK SCSK\n(skunk)",
                                "Bat LC\n(hoary bat)",
                                "Gannoruwa bat lyssavirus"
                     ))
p

png("plots/Figure 1.png", width = 7.5, height = 7.5, units = 'in', res = 600)
p
dev.off()
