library(stringr)
library(ggpubr)
library(seqinr)
library(ggplot2)

seqs = read.fasta("sequence_data/all_seqs.fasta",as.string = T)
seqs = seqs[1:length(seqs)-1]
metadata = read.csv("sequence_data/metadata.csv")

accession = c()
clade = c()
motif = c()
start_loc = c()

for(i in 1:length(seqs)){
  x = toupper(seqs[[i]][1])
  
  cpg_locs = str_locate_all(x, "CG")
  accession = append(accession, rep(names(seqs)[i], length(cpg_locs[[1]][,1])))
  motif = append(motif, rep("CpG", length(cpg_locs[[1]][,1])))
  clade = append(clade, rep(metadata$Clade[metadata$Accession==names(seqs)[i]], 
                            length(cpg_locs[[1]][,1])))
  start_loc = append(start_loc, cpg_locs[[1]][,1])
  
  ZAP_opt_locs = str_locate_all(x, c("C.......G.CG", "C....G.CG", "C.....G.CG",
                                 "C......G.CG", "C........G.CG"))
  accession = append(accession, rep(names(seqs)[i], length(ZAP_opt_locs[[1]][,1])))
  motif = append(motif, rep("ZAP optimal motifs", length(ZAP_opt_locs[[1]][,1])))
  clade = append(clade, rep(metadata$Clade[metadata$Accession==names(seqs)[i]], 
                            length(ZAP_opt_locs[[1]][,1])))
  start_loc = append(start_loc, ZAP_opt_locs[[1]][,1])
  
  ZAP_subopt_locs = str_locate_all(x, c("C.......C.CG", "C....C.CG", "C.....C.CG",
                                     "C......C.CG", "C........C.CG"))
  accession = append(accession, rep(names(seqs)[i], length(ZAP_subopt_locs[[1]][,1])))
  motif = append(motif, rep("ZAP suboptimal motifs", length(ZAP_subopt_locs[[1]][,1])))
  clade = append(clade, rep(metadata$Clade[metadata$Accession==names(seqs)[i]], 
                            length(ZAP_subopt_locs[[1]][,1])))
  start_loc = append(start_loc, ZAP_subopt_locs[[1]][,1])
  
}


df = as.data.frame(cbind(accession, clade, motif, start_loc))

clade = c()
motif = c()
prop = c()
start_loc = c()
for(i in unique(df$clade)){
  for(j in unique(df$motif)){
    for(k in unique(df$start_loc[df$clade == i & df$motif == j])){
      prop = append(prop, nrow(df[df$clade == i & df$motif == j & df$start_loc == k,])/
        length(unique(df$accession[df$clade == i & df$motif == j])))
      clade = append(clade, i)
      motif = append(motif, j)
      start_loc = append(start_loc, k)
    }
  }
}

df2 = as.data.frame(cbind(prop, clade, motif, start_loc))
df2$clade = factor(df2$clade, c("Cosmo AF1b", "Cosmo AM2a", "Arctic A", "Asian SEA2a", 
                              "Asian SEA2b", 
                              "Bat DR","Bat TB1", "Bat LC",
                              "Bat EF-E2","RAC-SK SCSK"))
mypal = c(rgb(0.9,0.9,0.9, alpha = 0.5),
          rgb(216/256, 27/256, 96/256, alpha = 0.8),
          rgb(30/256, 136/256, 229/256, alpha = 0.8))

g = ggplot(data = df2, aes(x = clade, y = as.numeric(start_loc), 
                       colour = motif, size = as.numeric(prop)))+
  geom_point()+
  scale_size_area()+
  scale_x_discrete(limits=rev)+
  scale_y_continuous(limits = c(0, 1353),
                     breaks = c(0,250, 500, 750, 1000, 1250, 1353))+
  coord_flip()+
  theme_pubclean()+
  scale_colour_manual(values = mypal)+
  labs(y = "Motif start location", x = "Clade",
    colour = "Motif", size = "Proportion of\nsequences")+
  theme(legend.position="bottom")
