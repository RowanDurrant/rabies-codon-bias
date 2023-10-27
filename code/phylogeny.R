#Phylogenetic tree
library(ggtree)
library(treeio)
tree = read.tree("trees/rooted_N/all_N_unique.fasta.treefile")
tipcolours = c()

library("readxl")
df = as.data.frame(read_excel("Codon_usage_N.xlsx"))
for(i in 3:ncol(df)){
  colnames(df)[i] = paste0(as.character(df[1,i]),"_",colnames(df)[i])
}
colnames(df)[2] = "Accession no."
df = df[2:nrow(df),]
df = rbind(df, c("Gannoruwa bat lyssavirus", "NC_031988_gannoruwa_outgroup"))
for(i in 1:length(tree$tip.label)){
  tipcolours[i] = df$Host[df$`Accession no.` == tree$tip.label[i]]
}
d <- data.frame(node=c(1:(Nnode(tree)+length(tree$tip.label))), 
                host = c(tipcolours, rep("black", Nnode(tree))))
d$host = factor(d$host, c("Big brown bat", "Hoary bat", "Free-tailed bat",
                              "Vampire bat", "Dog (SEA2a)", "Chinese ferret badger",
                              "Dog (AF1b)", "Mongoose", "Skunk", "Gannoruwa bat lyssavirus"))
p = ggtree(tree) + 
  theme_tree2()
p = p %<+% d  +  
  geom_tippoint(aes(color=host), size = 2) +
  scale_color_manual(values = c("#004949","#009292","#ff6db6","#ffb6db",
                                "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
                                "#920000","#924900","#db6d00","#24ff24","#ffff6d",
                                "#555555", "#999999"), name = "Host", guide = guide_legend()) +
  theme(legend.position="bottom") 
p


tree = read.tree("trees/rooted_G/all_G_unique.fasta.treefile")
tipcolours = c()

library("readxl")
df = as.data.frame(read_excel("Codon_usage_G.xlsx"))
for(i in 3:ncol(df)){
  colnames(df)[i] = paste0(as.character(df[1,i]),"_",colnames(df)[i])
}
colnames(df)[2] = "Accession no."
df = df[2:nrow(df),]
df = rbind(df, c("Gannoruwa bat lyssavirus", "NC_031988_gannoruwa_outgroup"))
for(i in 1:length(tree$tip.label)){
  tipcolours[i] = df$Host[df$`Accession no.` == tree$tip.label[i]]
}
d <- data.frame(node=c(1:(Nnode(tree)+length(tree$tip.label))), 
                host = c(tipcolours, rep("black", Nnode(tree))))
d$host = factor(d$host, c("Big brown bat", "Hoary bat", "Free-tailed bat",
                          "Vampire bat", "Dog (SEA2a)", "Chinese ferret badger",
                          "Dog (AF1b)", "Skunk", "Gannoruwa bat lyssavirus", "Mongoose"))
p = ggtree(tree) + 
  theme_tree2()
p = p %<+% d  +  
  geom_tippoint(aes(color=host), size = 2) +
  scale_color_manual(values = c("#004949","#009292","#ff6db6","#ffb6db",
                                "#490092","#006ddb","#b66dff","#b6dbff",
                                "#920000","#6db6ff","#924900","#db6d00","#24ff24","#ffff6d",
                                "#555555", "#999999"), name = "Host", guide = guide_legend()) +
  theme(legend.position="bottom") 
p
