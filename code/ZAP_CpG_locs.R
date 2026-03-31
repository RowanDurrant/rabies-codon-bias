library(stringr)
library(ggpubr)
library(seqinr)
library(ggplot2)

seqs = seqinr::read.fasta("sequence_data/all_seqs.fasta",as.string = T)
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
  
   tab_list = str_locate_all(x, c("C.......G.CG", "C....G.CG", "C.....G.CG",
                                 "C......G.CG", "C........G.CG"))
   ZAP_opt_locs = unname(c(tab_list[[1]][,1],
                           tab_list[[2]][,1],
                           tab_list[[3]][,1],
                           tab_list[[4]][,1],
                           tab_list[[5]][,1]))
  
  accession = append(accession, rep(names(seqs)[i], length(ZAP_opt_locs)))
  motif = append(motif, rep("ZAP optimal motifs", length(ZAP_opt_locs)))
  clade = append(clade, rep(metadata$Clade[metadata$Accession==names(seqs)[i]], 
                            length(ZAP_opt_locs)))
  start_loc = append(start_loc, ZAP_opt_locs)
}


df3 = as.data.frame(cbind(accession, clade, motif, start_loc))

clade = c()
motif = c()
prop = c()
start_loc = c()
for(i in unique(df3$clade)){
  for(j in unique(df3$motif)){
    for(k in unique(df3$start_loc[df3$clade == i & df3$motif == j])){
      prop = append(prop, length(unique(df3$accession[df3$clade == i & df3$motif == j & df3$start_loc == k]))/
        length(unique(df3$accession[df3$clade == i])))
      clade = append(clade, i)
      motif = append(motif, j)
      start_loc = append(start_loc, k)
    }
  }
}

df2 = as.data.frame(cbind(prop, clade, motif, start_loc))
df2$clade = factor(df2$clade, c("Cosmopolitan AF1b", "Cosmopolitan AM2a", "Arctic A", "Asian SEA2a", 
                              "Asian SEA2b", 
                              "Bats DR","Bats TB1", "Bats LC",
                              "Bats EF-E2","RAC-SK SCSK"))
mypal = c(rgb(0.9,0.9,0.9, alpha = 0.6),
          rgb(216/256, 27/256, 96/256, alpha = 0.6))

mylabels = c("Cosmo AF1b",
             "Cosmo AM2a",
             "Arctic A",
             "Asian SEA2a",
             "Asian SEA2b",
             "Bat DR",
             "Bat TB1",
             "Bat LC",
             "Bat EF-E2",
             "RAC-SK SCSK")

g = ggplot(data = df2, aes(x = clade, y = as.numeric(start_loc), 
                       colour = motif, size = as.numeric(prop)))+
  geom_point(stroke = 0)+
  scale_size_area()+
  scale_x_discrete(limits=rev,
                   labels = rev(mylabels))+
  scale_y_continuous(limits = c(0, 1353),
                     breaks = c(0,250, 500, 750, 1000, 1250, 1353))+
  scale_radius()+
  coord_flip()+
  theme_pubclean()+
  scale_colour_manual(values = mypal)+
  labs(y = "Motif start location", x = "Clade",
    colour = "Motif", size = "Proportion of\nsequences")+
  theme(legend.position="bottom")
g
