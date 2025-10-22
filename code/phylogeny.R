#Phylogenetic tree
library(ggtree)
library(treeio)
library(ggplot2)
library(ggbreak)
library(grid)
tree = read.tree("sequence_data/trees/all_seqs.fasta.treefile")
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
p = ggtree(tree) +
  theme_tree2()
p = p %<+% d  +  
   theme(legend.position="bottom") +
  ylim(0, 440) +
  geom_tippoint(aes(color=clade), size = 1) +
  scale_color_manual(values = c("#332288","#88CCEE","#CCDDAA","#44AA99","#117733",  
                                "#DDCC77","#999933", "#AA4499","#CC6677","#882255", "black"), 
                                name = "Clade", guide = guide_legend(),
                     labels = c("Cosmo AF1b\n(dog)",
                                "Cosmo AM2a\n(mongoose)",
                                "Arctic A\n(arctic fox)",
                                "Asian SEA2a\n(dog)",
                                "Asian SEA2b\n(CFB)",
                                "Bat DR\n(vampire bat)",
                                "Bat TB1\n(Mexican free\n-tailed bat)",
                                "Bat LC\n(hoary bat)",
                                "Bat EF-E2\n(big brown bat)",
                                "RAC-SK SCSK\n(skunk)",
                                "Gannoruwa bat\nlyssavirus"
                     ))
p

# png("plots/Figure 1.png", width = 7.5, height = 7.5, units = 'in', res = 600)
# p
# dev.off()


p2 = ggtree(tree) +
  ylim(0, 440) +
  #geom_text(aes(label=node), cex = 2, alpha = 0.5)
  geom_tippoint(size = 0.25, colour = "black")+
  # geom_hilight(node=546, fill = NA, color="#332288", alpha=0.4) +
  # geom_hilight(node=706, fill = NA, color="#88CCEE", alpha=0.4) +
  # geom_hilight(node=714, fill = NA, color="#CCDDAA", alpha=0.4) +
  # geom_hilight(node=436, fill = NA, color="#44AA99", alpha=0.4) +
  # geom_hilight(node=517, fill = NA, color="#44AA99", alpha=0.4) +
  # geom_hilight(node=483, fill = NA, color="#117733", alpha=0.4) +
  # geom_hilight(node=734, fill = NA, color="#DDCC77", alpha=0.4) +
  # geom_hilight(node=772, fill = NA, color="#999933", alpha=0.4) +
  # geom_hilight(node=803, fill = NA, color="#AA4499", alpha=0.4) +
  # geom_hilight(node=818, fill = NA, color="#CC6677", alpha=0.4) +
  # geom_hilight(node=853, fill = NA, color="#882255", alpha=0.4) +
  geom_cladelabel(node=431, label="NC_031988", 
                  color='black', fontsize=2)+
  geom_cladelabel(node=706, label="Cosmopolitan AM2a", 
                  color='#88CCEE', fontsize=3, 
                  offset = 0.005, offset.text = 0.005)+
  geom_cladelabel(node=714, label="Arctic A", 
                  color='#CCDDAA', fontsize=3, 
                  offset = 0.005, offset.text = 0.005)+
  geom_cladelabel(node=436, label="Asian SEA2a", 
                  color='#44AA99', fontsize=3, 
                  offset = 0.005, offset.text = 0.005)+
  geom_cladelabel(node=517, label="Asian SEA2a", 
                  color='#44AA99', fontsize=3, 
                  offset = 0.005, offset.text = 0.005)+
  geom_cladelabel(node=483, label="Asian SEA2b", 
                  color='#117733', fontsize=3, 
                  offset = 0.005, offset.text = 0.005)+
  geom_cladelabel(node=734, label="Bat DR", 
                  color='#DDCC77', fontsize=3, 
                  offset = 0.005, offset.text = 0.005)+
  geom_cladelabel(node=772, label="Bat TB1", 
                  color='#999933', fontsize=3, 
                  offset = 0.005, offset.text = 0.005)+
  geom_cladelabel(node=803, label="Bat LC", 
                  color='#AA4499', fontsize=3, 
                  offset = 0.005, offset.text = 0.005)+
  geom_cladelabel(node=818, label="Bat EF-E2", 
                  color='#CC6677', fontsize=3, 
                  offset = 0.005, offset.text = 0.005)+
  geom_cladelabel(node=853, label="RAC-SK SCSK", 
                  color='#882255', fontsize=3, 
                  offset = 0.005, offset.text = 0.005)+
  geom_cladelabel(node=546, label="Cosmopolitan AF1b", 
                  color='#332288', fontsize=3, 
                  offset = 0.005, offset.text = 0.005)+
  xlim(0, 0.465)+ geom_treescale(x = 0, y = 320, fontsize = 4, 
                                offset = 7)
p2

source("code/ENC.R")

library(ggmap)
p2 +
  inset(ggplotGrob(p), xmin = 0.24, xmax = 0.48, ymin = 20, ymax = 430)
