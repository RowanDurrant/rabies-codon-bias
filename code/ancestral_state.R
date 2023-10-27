tab=read.table('trees/rooted/all_seqs.fasta.state',header=TRUE)

node = c()
seqs = c()

for(i in unique(tab$Node)){
  node = append(node, i)
  seqs = append(seqs, paste(tab$State[tab$Node == i], collapse = ""))
}

#cpg calcs

CpG = function(x){
  x = unname(as.character(x))
  nG = str_count(x, "G")
  nC = str_count(x, "C")
  nCpG = str_count(x, "CG")
  N = nchar(x)
  ObsExpCpG = (nCpG/(nC * nG))*N
  return(ObsExpCpG)
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
gc = c()
cpg_actual = c()
accessions = node

for(j in 1:length(seqs)){
  cpg = append(cpg, CpG(seqs[j]))
  gc = append(gc, GCcontent(seqs[j]))
  cpg_actual = append(cpg_actual, CpG_actual(seqs[j]))
  
}
internal_nodes = data.frame(accessions, cpg, gc, cpg_actual)


seqs = readDNAStringSet("trees/rooted/all_seqs.fasta")
for(j in 1:length(seqs)){
  cpg = append(cpg, CpG(seqs[j]))
  gc = append(gc, GCcontent(seqs[j]))
  cpg_actual = append(cpg_actual, CpG_actual(seqs[j]))
  accessions = append(accessions, names(seqs[j]))
  
}
tips = data.frame(accessions, cpg, gc, cpg_actual)

all = rbind(tips, internal_nodes)

#tree

library(ggtree)
library(treeio)
library(viridis)

tree = read.tree("trees/rooted/all_seqs.fasta.treefile")
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
  geom_point(aes(color=cpg), size = 2) +
  scale_color_viridis(name = "Obs/Exp CpG", guide = "colorbar", 
                      breaks =c(0.425, 0.525, 0.625)) +
  theme(legend.position="bottom")
p

