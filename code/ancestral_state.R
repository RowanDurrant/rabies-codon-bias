library("stringr")
library(seqinr)
library(ggtree)
library(treeio)
library(viridis)
library(ggplot2)
library(ggbreak)

tab=read.table('sequence_data/all_outgroups/all_outgroups_aligned.fasta.state',header=TRUE)

reps = 100
node = rep(unique(tab$Node), each = reps)
seqs = c()

for(i in unique(tab$Node)){
  print(i)
  variable_sites = tab$Site[tab$Node == i & 
                              tab$p_A != 1 & 
                              tab$p_C != 1 & 
                              tab$p_G != 1 & 
                              tab$p_T != 1]
  for(j in 1:reps){
    newseq = tab$State[tab$Node == i]
    for(k in variable_sites){
      probs = unlist(unname(tab[tab$Node == i &
                                  tab$Site == k, 4:7]))
      newseq[k] = sample(c("A", "C", "G", "T"), 1,
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

cpg = internal_nodes$cpg
upa = internal_nodes$upa
accessions = internal_nodes$accessions

seqs = seqinr::read.fasta("sequence_data/all_outgroups/all_outgroups.fasta", as.string = T)
for(j in 1:length(seqs)){
  cpg = append(cpg, CpG(seqs[[j]][1]))
  upa = append(upa, UpA(seqs[[j]][1]))
  accessions = append(accessions, names(seqs[j]))
  
}
tips = data.frame(accessions, cpg, upa)
all = tips

write.csv(internal_nodes_long, "output_data/asr_values.csv")

#tree
tree = read.tree("sequence_data/all_outgroups/all_outgroups_aligned.fasta.treefile")
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
                      breaks = c(0.4,0.5,0.6)) +
  geom_cladelabel(node=(416+82), label="", 
                  color='#332288', fontsize=3, 
                  align = T, barsize = 1.5, offset = 0.01)+
  geom_cladelabel(node=(416+74), label="", 
                  color='#88CCEE', fontsize=3, 
                  align = T, barsize = 1.5, offset = 0.01)+
  geom_cladelabel(node=(416+5), label="", 
                  color='#CCDDAA', fontsize=3, 
                  align = T, barsize = 1.5, offset = 0.01)+
  geom_cladelabel(node=(416+199), label="", 
                  color='#44AA99', fontsize=3, 
                  align = T, barsize = 1.5, offset = 0.01)+
  geom_cladelabel(node=(416+274), label="", 
                  color='#44AA99', fontsize=3, 
                  align = T, barsize = 1.5, offset = 0.01)+
  geom_cladelabel(node=(416+246), label="", 
                  color='#117733', fontsize=3, 
                  align = T, barsize = 1.5, offset = 0.01)+
  geom_cladelabel(node=(416+296), label="", 
                  color='#DDCC77', fontsize=3, 
                  align = T, barsize = 1.5, offset = 0.01)+
  geom_cladelabel(node=(416+335), label="", 
                  color='#999933', fontsize=3, 
                  align = T, barsize = 1.5, offset = 0.01)+
  geom_cladelabel(node=(416+365), label="", 
                  color='#AA4499', fontsize=3, 
                  align = T, barsize = 1.5, offset = 0.01)+
  geom_cladelabel(node=(416+380), label="", 
                  color='#CC6677', fontsize=3, 
                  align = T, barsize = 1.5, offset = 0.01)+
  geom_cladelabel(node=(416+285), label="", 
                  color='#882255', fontsize=3, 
                  align = T, barsize = 1.5, offset = 0.01)+ 
  scale_x_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5))+
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
  scale_color_viridis(name = "Obs/Exp UpA", guide = "colorbar") +
  geom_cladelabel(node=(416+82), label="Cosmo\nAF1b", 
                  color='#332288', fontsize=3, 
                  align = T, barsize = 1.5, offset = 0.01)+
  geom_cladelabel(node=(416+74), label="Cosmo AM2a", 
                  color='#88CCEE', fontsize=3, 
                  align = T, barsize = 1.5, offset = 0.01)+
  geom_cladelabel(node=(416+5), label="Arctic A", 
                  color='#CCDDAA', fontsize=3, 
                  align = T, barsize = 1.5, offset = 0.01)+
  geom_cladelabel(node=(416+199), label="Asian\nSEA2a", 
                  color='#44AA99', fontsize=3, 
                  align = T, barsize = 1.5, offset = 0.01)+
  geom_cladelabel(node=(416+274), label="Asian SEA2a", 
                  color='#44AA99', fontsize=3, 
                  align = T, barsize = 1.5, offset = 0.01)+
  geom_cladelabel(node=(416+246), label="Asian\nSEA2b", 
                  color='#117733', fontsize=3, 
                  align = T, barsize = 1.5, offset = 0.01)+
  geom_cladelabel(node=(416+296), label="Bat DR", 
                  color='#DDCC77', fontsize=3, 
                  align = T, barsize = 1.5, offset = 0.01)+
  geom_cladelabel(node=(416+335), label="Bat TB1", 
                  color='#999933', fontsize=3, 
                  align = T, barsize = 1.5, offset = 0.01)+
  geom_cladelabel(node=(416+365), label="Bat LC", 
                  color='#AA4499', fontsize=3, 
                  align = T, barsize = 1.5, offset = 0.01)+
  geom_cladelabel(node=(416+380), label="Bat EF-E2", 
                  color='#CC6677', fontsize=3, 
                  align = T, barsize = 1.5, offset = 0.01)+
  geom_cladelabel(node=(416+285), label="RAC-SK\nSCSK", 
                  color='#882255', fontsize=3, 
                  align = T, barsize = 1.5, offset = 0.01)+ 
  scale_x_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5),
                     limits = c(NA,0.58))+
  theme(legend.position="bottom",
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank())
p2

ggpubr::ggarrange(p,p2)
png("plots/Figure 6.png", width = 11, height = 7, units = 'in', res = 600)
ggpubr::ggarrange(p,p2)
dev.off()

internal_nodes[internal_nodes$accessions == "Node3",] #RABV MRCA
internal_nodes[internal_nodes$accessions == "Node4",] #Carnivore MRCA
internal_nodes[internal_nodes$accessions == "Node285",] #Bat MRCA
internal_nodes[internal_nodes$accessions == "Node294",] #Bat MRCA without RAC-SK
