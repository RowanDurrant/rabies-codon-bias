#cpg calcs
library("stringr")
library(seqinr)
library(ggtree)
library(treeio)
library(viridis)
library(ggplot2)
library(ggbreak)

tab=read.table('sequence_data/gannoruwa_outgroup/all_seqs.fasta.state',header=TRUE)

reps = 100
node = rep(unique(tab$Node), each = reps)
seqs = c()

for(i in unique(tab$Node)){
  print(i)
  for(j in 1:reps){
    newseq = tab$State[tab$Node == i]
    variable_sites = tab$Site[tab$Node == i & 
                                tab$p_A != 1 & 
                                tab$p_C != 1 & 
                                tab$p_G != 1 & 
                                tab$p_T != 1]
    
    for(k in variable_sites){
      probs = unname(tab[tab$Node == i & 
                          tab$Site == k, 4:7])

      newseq[k] = sample(c("A", "C", "T", "G"), 1,
                                          prob = probs)
    }
    seqs = append(seqs, paste(newseq, collapse = ""))
  }
}

CpG = function(x){
  x = toupper(unname(as.character(x)))
  nG = str_count(x, "G")
  nC = str_count(x, "C")
  nCpG = str_count(x, "CG")
  N = nchar(x)
  ObsExpCpG = (nCpG/(nC * nG))*N
  return(ObsExpCpG)
}

UpA = function(x){
  x = toupper(unname(as.character(x)))
  nT = str_count(x, "T")
  nA = str_count(x, "A")
  nUpA = str_count(x, "TA")
  N = nchar(x)
  ObsExpUpA = (nUpA/(nT * nA))*N
  return(ObsExpUpA)
}

cpg = c()
upa = c()
accessions = node

for(j in 1:length(seqs)){
  cpg = append(cpg, CpG(seqs[j]))
  upa = append(upa, UpA(seqs[j]))
}
internal_nodes_long = data.frame(accessions, cpg, upa)

cpg = c()
cpg_upper = c()
cpg_lower = c()
upa = c()
upa_upper = c()
upa_lower = c()
accessions = c()
for(i in unique(internal_nodes_long$accessions)){
  accessions = append(accessions, i)
  
  all_cpg = internal_nodes_long$cpg[internal_nodes_long$accessions == i]
  cpg = append(cpg, mean(all_cpg))
  cpg_lower = append(cpg_lower, unname(quantile(all_cpg, 0.025)))
  cpg_upper = append(cpg_upper, unname(quantile(all_cpg, 0.975)))
  
  all_upa = internal_nodes_long$upa[internal_nodes_long$accessions == i]
  upa = append(upa, mean(all_upa))
  upa_lower = append(upa_lower, unname(quantile(all_upa, 0.025)))
  upa_upper = append(upa_upper, unname(quantile(all_upa, 0.975)))
}
internal_nodes = data.frame(accessions, cpg, cpg_lower, cpg_upper,
                               upa, upa_lower, upa_upper)

seqs = seqinr::read.fasta("sequence_data/all_seqs.fasta", as.string = T)
for(j in 1:length(seqs)){
  cpg = append(cpg, CpG(seqs[[j]][1]))
  upa = append(upa, UpA(seqs[[j]][1]))
  gc = append(gc, GCcontent(seqs[[j]][1]))
  cpg_actual = append(cpg_actual, CpG_actual(seqs[[j]][1]))
  accessions = append(accessions, names(seqs[j]))
  
}
tips = data.frame(accessions, cpg, upa, gc, cpg_actual)
all = tips

#tree
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

d$cpg[d$node == 432] #common ancestor
d$cpg[d$node == 433] #carnivore ancestor
d$cpg[d$node == 730] #bat ancestor

d2$upa[d2$node == 432] #common ancestor
d2$upa[d2$node == 433] #carnivore ancestor
d2$upa[d2$node == 730] #bat ancestor
