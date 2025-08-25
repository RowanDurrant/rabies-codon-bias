library(ggplot2)
library(ggpubr)
library(ggtree)
library(treeio)

df = read.csv("output_data/myseqs_RSDUc_cpg_tpa.tsv", sep='\t', row.names = 1, header= TRUE)
df$Accession = rownames(df)
meta = read.csv("sequence_data/metadata.csv")
meta = rbind(meta, c("gannoruwa_outgroup", NA, "Gannoruwa bat lyssavirus", "bat", "bat"))

df = merge(meta, df, by = "Accession")

tree = read.tree("sequence_data/trees/all_seqs.fasta.treefile")

tipcolours = c()
for(i in 1:length(tree$tip.label)){
  if(tree$tip.label[i] == "gannoruwa_outgroup"){tipcolours[i] = "Gannoruwa bat lyssavirus"}
  else{tipcolours[i] = df$Clade[df$Accession == tree$tip.label[i]]}
}
d <- data.frame(node=c(1:(Nnode(tree)+length(tree$tip.label))), 
                clade = c(tipcolours, rep("black", Nnode(tree))))
d$clade = factor(d$clade, c("Cosmo AF1b", "Cosmo AM2a", "Asian SEA2a", 
                            "Asian SEA2b", "Arctic A",
                            "Bat TB1",
                            "Bat DR", "Bat EF-E2","RAC-SK SCSK", "Bat LC",
                            "Gannoruwa bat lyssavirus"
))

p = ggtree(tree, ladderize = F)+
  theme_tree2() 
p = p %<+% d  +  
  theme(legend.position="bottom") +
  ylim(0, 431) +
  geom_tippoint(aes(color=clade), size = 1) +
  scale_color_manual(values = c("#332288","#88CCEE","#44AA99","#117733", "#CCDDAA", 
                                "#999933", "#DDCC77","#CC6677","#882255","#AA4499", "black"), 
                     name = "Clade", guide = guide_legend(),
                     labels = c("Cosmo AF1b\n(dog)",
                                "Cosmo AM2a\n(mongoose)",
                                "Asian SEA2a\n(dog)",
                                "Asian SEA2b\n(CFB)",
                                "Arctic A\n(arctic fox)",
                                "Bat TB1\n(Mexican free-tailed bat)",
                                "Bat DR\n(vampire bat)",
                                "Bat EF-E2\n(big brown bat)",
                                "RAC-SK SCSK\n(skunk)",
                                "Bat LC\n(hoary bat)",
                                "Gannoruwa bat\nlyssavirus"
                     ))

d3 <- data.frame(id=df$Accession, value = df$TpApos2)
p3 <- facet_plot(p, panel='UpApos2', data=d3, geom=geom_segment, 
                 aes(x=0, xend=value, y=y, yend=y, colour = tipcolours))

d4 <- data.frame(id=df$Accession, value = df$TpAbridge)
p4 <- facet_plot(p3, panel='UpAbridge', data=d4, geom=geom_segment, 
                 aes(x=0, xend=value, y=y, yend=y, colour = tipcolours))

#CpG
d5 <- data.frame(id=df$Accession, value = df$CpGpos1)
p5 <- facet_plot(p4, panel='CpGpos1', data=d5, geom=geom_segment, 
                 aes(x=0, xend=value, y=y, yend=y, colour = tipcolours))

d6 <- data.frame(id=df$Accession, value = df$CpGpos2)
p6 <- facet_plot(p5, panel='CpGpos2', data=d6, geom=geom_segment, 
                 aes(x=0, xend=value, y=y, yend=y, colour = tipcolours))

d7 <- data.frame(id=df$Accession, value = df$CpGbridge)
p7 <- facet_plot(p6, panel='CpGbridge', data=d7, geom=geom_segment, 
                 aes(x=0, xend=value, y=y, yend=y, colour = tipcolours))
p7+
  theme(panel.spacing.x = unit(3, "mm"))

