tab=read.table('sequence_data/trees/all_seqs.fasta.state',header=TRUE)

node = c()
seqs = c()

for(i in unique(tab$Node)){
  node = append(node, i)
  seqs = append(seqs, paste(tab$State[tab$Node == i], collapse = ""))
}

#cpg calcs
library("stringr")

CpG = function(x){
  x = unname(as.character(x))
  nG = str_count(x, "G")
  nC = str_count(x, "C")
  nCpG = str_count(x, "CG")
  N = nchar(x)
  ObsExpCpG = (nCpG/(nC * nG))*N
  return(ObsExpCpG)
}

UpA = function(x){
  x = unname(as.character(x))
  nT = str_count(x, "T")
  nA = str_count(x, "A")
  nUpA = str_count(x, "TA")
  N = nchar(x)
  ObsExpUpA = (nUpA/(nT * nA))*N
  return(ObsExpUpA)
}

GCcontent = function(x){
  x = unname(as.character(x))
  nG = str_count(x, "G")
  nC = str_count(x, "C")
  N = nchar(x)
  gccont = (nC + nG)/N
  return(gccont)
  
}

CpG_actual = function(x){
  x = unname(as.character(x))
  nCpG = str_count(x, "CG")
  return(nCpG)
}

cpg = c()
upa = c()
gc = c()
cpg_actual = c()
accessions = node

for(j in 1:length(seqs)){
  cpg = append(cpg, CpG(seqs[j]))
  upa = append(upa, UpA(seqs[j]))
  gc = append(gc, GCcontent(seqs[j]))
  cpg_actual = append(cpg_actual, CpG_actual(seqs[j]))
  
}
internal_nodes = data.frame(accessions, cpg, upa, gc, cpg_actual)

library(Biostrings)
seqs = readDNAStringSet("sequence_data/all_seqs.fasta")
for(j in 1:length(seqs)){
  cpg = append(cpg, CpG(seqs[j]))
  upa = append(upa, UpA(seqs[j]))
  gc = append(gc, GCcontent(seqs[j]))
  cpg_actual = append(cpg_actual, CpG_actual(seqs[j]))
  accessions = append(accessions, names(seqs[j]))
  
}
tips = data.frame(accessions, cpg, upa, gc, cpg_actual)
all = tips
#tree

library(ggtree)
library(treeio)
library(viridis)
library(ggplot2)
library(ggbreak)

tree = read.tree("sequence_data/trees/all_seqs.fasta.treefile")
tipcolours = c()
internal_nodes = internal_nodes[order(match(internal_nodes$accessions,
                                            tree$node.label)),]
for(i in 1:length(tree$tip.label)){
  tipcolours[i] = all$cpg[all$accessions == tree$tip.label[i]]
}
d <- data.frame(node=c(1:(Nnode(tree)+length(tree$tip.label))), 
                cpg = c(tipcolours, internal_nodes$cpg))
p = ggtree(tree) + 
  theme_tree2()
p = p %<+% d  +  
  ylim (0,431)+
  geom_point(aes(color=cpg), size = 2) +
  scale_color_viridis(name = "Obs/Exp CpG", guide = "colorbar",
                      breaks = c(0.45, 0.55, 0.65)) +
  xlim(0,0.46)+ scale_x_break(c(0.21,0.395), ticklabels = c(0.4,0.45))+
  theme(legend.position="bottom",
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank())
p

tipcolours = c()
internal_nodes = internal_nodes[order(match(internal_nodes$accessions,
                                            tree$node.label)),]
for(i in 1:length(tree$tip.label)){
  tipcolours[i] = all$upa[all$accessions == tree$tip.label[i]]
}
d2 <- data.frame(node=c(1:(Nnode(tree)+length(tree$tip.label))), 
                upa = c(tipcolours, internal_nodes$upa))
p2 = ggtree(tree) + 
  theme_tree2()
p2 = p2 %<+% d2  +  
  ylim (0,431)+
  geom_point(aes(color=upa), size = 2) +
  scale_color_viridis(name = "Obs/Exp UpA", guide = "colorbar",
                      breaks = c(0.6, 0.65,0.7)) +
  scale_x_break(c(0.205,0.395), ticklabels = c(0.4,0.45))+
  xlim(0.46,0)+ 
  theme(legend.position="bottom",
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank())
p2

p+p2

png("plots/Figure 8.png", width = 7, height = 7, units = 'in', res = 600)
p+p2
dev.off()
