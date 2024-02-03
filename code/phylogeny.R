#Phylogenetic tree
library(ggtree)
library(treeio)
tree = read.tree("trees/rooted_N/all_seqs.fasta.treefile")
tipcolours = c()

library("readxl")
df = as.data.frame(read_excel("data/Codon_usage_N.xlsx"))
for(i in 3:ncol(df)){
  colnames(df)[i] = paste0(as.character(df[1,i]),"_",colnames(df)[i])
}
colnames(df)[2] = "Accession no."
df = df[2:nrow(df),]
df = rbind(df, c("Gannoruwa bat lyssavirus", "gannoruwa_outgroup"))
for(i in 1:length(tree$tip.label)){
  tipcolours[i] = df$Host[df$`Accession no.` == tree$tip.label[i]]
}
d <- data.frame(node=c(1:(Nnode(tree)+length(tree$tip.label))), 
                host = c(tipcolours, rep("black", Nnode(tree))))
d$host = factor(d$host, c("Dog (AF1b)", "Mongoose","Dog (SEA2a)", "Chinese ferret badger",
                          "Free-tailed bat",
                          "Vampire bat", "Big brown bat","Skunk", "Hoary bat",
                          "Gannoruwa bat lyssavirus"
))
p = ggtree(tree) + 
  theme_tree2()
p = p %<+% d  +  
   theme(legend.position="bottom") +
  ylim (0,420)+
  geom_tippoint(aes(color=host), size = 2) +
  scale_color_manual(values = c("#332288","#88CCEE","#44AA99","#117733","#999933",
                                "#DDCC77","#CC6677","#882255","#AA4499", "black"), name = "Clade", guide = guide_legend(),
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

png("plots/Figure 1.png", width = 7.5, height = 5, units = 'in', res = 600)
p
dev.off()
