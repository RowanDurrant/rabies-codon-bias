##CpG 
# Observed to Expected CpG is calculated as below : 
# Obs/Exp CpG = Number of CpG * N / (Number of C * Number of G) 
# where N = length of sequence.

library(stringr)

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

library("Biostrings")

seqs = readDNAStringSet("sequence_data/all_seqs.fasta")
seqs = seqs[1:length(seqs)-1]
metadata = read.csv("sequence_data/metadata.csv")

cpg = c()
gc = c()
cpg_actual = c()
accessions = c()
clade = c()
host_group = c()

  for(j in 1:length(seqs)){
       cpg = append(cpg, CpG(seqs[j]))
       gc = append(gc, GCcontent(seqs[j]))
       cpg_actual = append(cpg_actual, CpG_actual(seqs[j]))
       accessions = append(accessions, names(seqs[j]))
       clade = append(clade, metadata$Clade[metadata$Accession==names(seqs[j])])
       host_group = append(host_group, metadata$Group[metadata$Accession==names(seqs[j])])
    }
  
df = data.frame(accessions, clade, host_group, cpg, gc, cpg_actual)
#write.csv(df, "N_CpG.csv")
df$clade = factor(df$clade, c("Cosmo AF1b", "Cosmo AM2a","Asian SEA2a", 
                            "Asian SEA2b", 
                          "Bat TB1",
                          "Bat DR", "Bat EF-E2","RAC-SK SCSK", "Bat LC"))

library(ggplot2)
p1= ggplot(data = df, aes(x = clade, y = cpg))+
  geom_boxplot()+ 
  geom_jitter(aes(color = clade), size = 0.5, 
              width = 0.4, height = 0) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")+
  ylab("Obs/Exp CpG") + xlab("Clade") +
  scale_color_manual(values = c("#332288","#88CCEE","#44AA99","#117733","#999933",
                                "#DDCC77","#CC6677","#882255","#AA4499"), 
                     name = "Clade", guide = guide_legend(),
                     labels = c("Cosmo AF1b\n(dog)",
                                "Cosmo AM2a\n(mongoose)",
                                "Asian SEA2a\n(dog)",
                                "Asian SEA2b\n(CFB)",
                                "Bat TB1\n(Mexican free-tailed bat)",
                                "Bat DR\n(vampire bat)",
                                "Bat EF-E2\n(big brown bat)",
                                "RAC-SK SCSK\n(skunk)",
                                "Bat LC\n(hoary bat)"
                                )) +
  scale_x_discrete(labels = c("Cosmo AF1b\n(dog)",
                              "Cosmo AM2a\n(mongoose)",
                              "Asian SEA2a\n(dog)",
                              "Asian SEA2b\n(CFB)",
                              "Bat TB1\n(Mexican free-tailed bat)",
                              "Bat DR\n(vampire bat)",
                              "Bat EF-E2\n(big brown bat)",
                              "RAC-SK SCSK\n(skunk)",
                              "Bat LC\n(hoary bat)"
  ))

p2= ggplot(data = df, aes(x = clade, y = gc))+
  geom_boxplot()+ 
  geom_jitter(aes(color = clade), size  = 0.5, 
              width = 0.4, height = 0) + theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")+
  ylab("GC content") + xlab("Clade") +
  scale_color_manual(values = c("#332288","#88CCEE","#44AA99","#117733","#999933",
                                "#DDCC77","#CC6677","#882255","#AA4499"), 
                     name = "Clade", guide = guide_legend(),
                     labels = c("Cosmo AF1b\n(dog)",
                                "Cosmo AM2a\n(mongoose)",
                                "Asian SEA2a\n(dog)",
                                "Asian SEA2b\n(CFB)",
                                "Bat TB1\n(Mexican free-tailed bat)",
                                "Bat DR\n(vampire bat)",
                                "Bat EF-E2\n(big brown bat)",
                                "RAC-SK SCSK\n(skunk)",
                                "Bat LC\n(hoary bat)"
                     )) +
  scale_x_discrete(labels = c("Cosmo AF1b\n(dog)",
                              "Cosmo AM2a\n(mongoose)",
                              "Asian SEA2a\n(dog)",
                              "Asian SEA2b\n(CFB)",
                              "Bat TB1\n(Mexican free\n-tailed bat)",
                              "Bat DR\n(vampire bat)",
                              "Bat EF-E2\n(big brown bat)",
                              "RAC-SK SCSK\n(skunk)",
                              "Bat LC\n(hoary bat)"
  )) 

p3 = ggplot(data = df, aes(x = clade, y = cpg_actual))+
  geom_boxplot()+ 
  geom_jitter(aes(color = clade), size  = 0.5, 
              width = 0.4, height = 0) + theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")+
  ylab("No. CpGs") + xlab("Clade") +
  scale_color_manual(values = c("#332288","#88CCEE","#44AA99","#117733","#999933",
                                "#DDCC77","#CC6677","#882255","#AA4499"), 
                     name = "Clade", guide = guide_legend(),
                     labels = c("Cosmo AF1b\n(dog)",
                                "Cosmo AM2a\n(mongoose)",
                                "Asian SEA2a\n(dog)",
                                "Asian SEA2b\n(CFB)",
                                "Bat TB1\n(Mexican free-tailed bat)",
                                "Bat DR\n(vampire bat)",
                                "Bat EF-E2\n(big brown bat)",
                                "RAC-SK SCSK\n(skunk)",
                                "Bat LC\n(hoary bat)"
                     )) +
  scale_x_discrete(labels = c("Cosmo AF1b\n(dog)",
                              "Cosmo AM2a\n(mongoose)",
                              "Asian SEA2a\n(dog)",
                              "Asian SEA2b\n(CFB)",
                              "Bat TB1\n(Mexican free\n-tailed bat)",
                              "Bat DR\n(vampire bat)",
                              "Bat EF-E2\n(big brown bat)",
                              "RAC-SK SCSK\n(skunk)",
                              "Bat LC\n(hoary bat)"
  ))
library(ggpubr)
ggarrange(p1, ggarrange(p2, p3, labels = c("B", "C"), ncol = 1, nrow = 2), labels = c("A"))

png("plots/Figure 7.png", width = 8, height = 5, units = 'in', res = 600)
ggarrange(p1, ggarrange(p2, p3, labels = c("B", "C"), ncol = 1, nrow = 2), labels = c("A"))

dev.off()