##CpG 
# Observed to Expected CpG is calculated as below : 
# Obs/Exp CpG = Number of CpG * N / (Number of C * Number of G) 
# where N = length of sequence.

library(stringr)
library(ggpubr)
library("Biostrings")
library(ggplot2)
library(egg)

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

ZAP_optimal = function(x){
  x = unname(as.character(x))
  n_optimal = str_count(x, pattern = "C.......G.CG")
  n_optimal = n_optimal + str_count(x, pattern = "C....G.CG")
  n_optimal = n_optimal + str_count(x, pattern = "C.....G.CG")
  n_optimal = n_optimal + str_count(x, pattern = "C......G.CG")
  n_optimal = n_optimal + str_count(x, pattern = "C........G.CG")
  return(n_optimal)
}

ZAP_suboptimal = function(x){
  x = unname(as.character(x))
  n_optimal = str_count(x, pattern = "C.......C.CG")
  n_optimal = n_optimal + str_count(x, pattern = "C....C.CG")
  n_optimal = n_optimal + str_count(x, pattern = "C.....C.CG")
  n_optimal = n_optimal + str_count(x, pattern = "C......C.CG")
  n_optimal = n_optimal + str_count(x, pattern = "C........C.CG")
  return(n_optimal)
}

TpA = function(x){
  x = unname(as.character(x))
  nT = str_count(x, "T")
  nA = str_count(x, "A")
  nTpA = str_count(x, "TA")
  N = nchar(x)
  ObsExpTpA = (nTpA/(nT * nA))*N
  return(ObsExpTpA)
}

TAcontent = function(x){
  x = unname(as.character(x))
  nT = str_count(x, "T")
  nA = str_count(x, "A")
  N = nchar(x)
  tacont = (nT + nA)/N
  return(tacont)
  
}

TpA_actual = function(x){
  x = unname(as.character(x))
  nTpA = str_count(x, "TA")
  return(nTpA)
}

seqs = readDNAStringSet("sequence_data/all_seqs.fasta")
seqs = seqs[1:length(seqs)-1]
metadata = read.csv("sequence_data/metadata.csv")

cpg = c()
gc = c()
cpg_actual = c()
accessions = c()
clade = c()
host_group = c()
ZAP_optimal_motifs = c()
ZAP_suboptimal_motifs = c()
tpa = c()
ta = c()
tpa_actual = c()

  for(j in 1:length(seqs)){
       cpg = append(cpg, CpG(seqs[j]))
       tpa = append(tpa, TpA(seqs[j]))
       gc = append(gc, GCcontent(seqs[j]))
       ta = append(ta, TAcontent(seqs[j]))
       cpg_actual = append(cpg_actual, CpG_actual(seqs[j]))
       tpa_actual= append(tpa_actual, TpA_actual(seqs[j]))
       accessions = append(accessions, names(seqs[j]))
       clade = append(clade, metadata$Clade[metadata$Accession==names(seqs[j])])
       host_group = append(host_group, metadata$Group[metadata$Accession==names(seqs[j])])
       ZAP_optimal_motifs = append(ZAP_optimal_motifs, ZAP_optimal(seqs[j]))
       ZAP_suboptimal_motifs = append(ZAP_suboptimal_motifs, ZAP_suboptimal(seqs[j]))
    }
  
df = data.frame(accessions, clade, host_group, cpg, gc, cpg_actual, 
                ZAP_optimal_motifs, ZAP_suboptimal_motifs, tpa, ta, tpa_actual)
write.csv(df, "output_data/N_CpG.csv")
df$clade = factor(df$clade, c("Cosmo AF1b", "Cosmo AM2a", "Arctic A", "Asian SEA2a", 
                            "Asian SEA2b", 
                            "Bat DR","Bat TB1", "Bat LC",
                           "Bat EF-E2","RAC-SK SCSK"))

mypal = c("#332288","#88CCEE","#CCDDAA","#44AA99","#117733",  
          "#DDCC77","#999933", "#AA4499","#CC6677","#882255")
mylabels = c("Cosmo AF1b\n(dog)",
             "Cosmo AM2a\n(mongoose)",
             "Arctic A\n(arctic fox)",
             "Asian SEA2a\n(dog)",
             "Asian SEA2b\n(CFB)",
             "Bat DR\n(vampire bat)",
             "Bat TB1\n(Mexican free\n-tailed bat)",
             "Bat LC\n(hoary bat)",
             "Bat EF-E2\n(big brown bat)",
             "RAC-SK SCSK\n(skunk)")

p1= ggplot(data = df, aes(x = clade, y = cpg, fill = clade))+
  geom_boxplot(outlier.size = 0.1, width = 0.5)+
  geom_jitter(colour = "black", size = 0.5,
              width = 0.25, height = 0, alpha= 0.6) + theme_bw()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  ylab("Obs/Exp CpG") + xlab("Clade") +
  scale_fill_manual(values = mypal,
                    name = "Clade", guide = guide_legend(),
                    labels = mylabels)

# p2= ggplot(data = df, aes(x = clade, y = gc, fill = clade))+
#   geom_boxplot()+ 
#   geom_jitter(aes(color = clade), size = 1,
#               width = 0.4, height = 0, alpha= 0.6) + theme_bw()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         legend.position = "none")+
#   ylab("GC content") + xlab("Clade") +
#   scale_fill_manual(values = mypal, 
#                      name = "Clade", guide = guide_legend(),
#                      labels = mylabels) +
#   scale_x_discrete(labels = mylabels) 

p3 = ggplot(data = df, aes(x = clade, y = cpg_actual, fill = clade))+
  geom_boxplot(outlier.size = 0.1, width = 0.5)+
  geom_jitter(colour = "black", size = 0.5,
              width = 0.25, height = 0, alpha= 0.6) + theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")+
  ylab("No. CpGs") + xlab("Clade") +
  scale_fill_manual(values = mypal,
                    name = "Clade", guide = guide_legend(),
                    labels = mylabels) +
  scale_x_discrete(labels = mylabels)

# ggarrange(p1, ggarrange(p2, p3, labels = c("B", "C"), ncol = 1, nrow = 2), labels = c("A"))
# 
# png("plots/Figure 7.png", width = 8, height = 5, units = 'in', res = 600)
# ggarrange(p1, ggarrange(p2, p3, labels = c("B", "C"), ncol = 1, nrow = 2), labels = c("A"))
# 
# dev.off()

p6 = ggplot(data = df, aes(x = clade, y = tpa, fill = clade))+
  geom_boxplot(outlier.size = 0.1, width = 0.5)+
  geom_jitter(colour = "black", size = 0.5,
              width = 0.25, height = 0, alpha= 0.6) + theme_bw()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  ylab("Obs/Exp UpA") + xlab("Clade") +
  scale_fill_manual(values = mypal,
                    name = "Clade", guide = guide_legend(),
                    labels = mylabels)


# p7= ggplot(data = df, aes(x = clade, y = ta))+
#   geom_boxplot()+
#   geom_jitter(aes(color = clade), size = 1,
#               width = 0.4, height = 0, alpha= 0.6) + theme_bw()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         legend.position = "none")+
#   ylab("UA content") + xlab("Clade") +
#   scale_color_manual(values = mypal,
#                      name = "Clade", guide = guide_legend(),
#                      labels = mylabels) +
#   scale_x_discrete(labels = mylabels)

p8 = ggplot(data = df, aes(x = clade, y = tpa_actual, fill = clade))+
  geom_boxplot(outlier.size = 0.1, width = 0.5)+
  geom_jitter(colour = "black", size = 0.5,
              width = 0.25, height = 0, alpha= 0.6) + theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")+
  ylab("No. UpAs") + xlab("Clade") +
  scale_fill_manual(values = mypal,
                     name = "Clade", guide = guide_legend(),
                     labels = mylabels) +
  scale_x_discrete(labels = mylabels)
#ggarrange(p6, ggarrange(p7, p8, labels = c("B", "C"), ncol = 1, nrow = 2), labels = c("A"))

# png("plots/tpa_fig_7.png", width = 8, height = 5, units = 'in', res = 600)
# ggarrange(p6, ggarrange(p7, p8, labels = c("B", "C"), ncol = 1, nrow = 2), labels = c("A"))
# 
# dev.off()

egg::ggarrange(p1,p6,p3,p8, labels = c("A", "B", "C", "D"))

png("plots/Figure 6.png", width = 8, height = 8, units = 'in', res = 600)
egg::ggarrange(p1,p6,p3,p8, labels = c("A", "B", "C", "D"))
dev.off()


library(tidyr)
df$ratio = df$ZAP_optimal_motifs/df$ZAP_suboptimal_motifs

p4 = ggplot(data = df, aes(x = clade, y = ratio, fill = clade))+
  geom_hline(yintercept = 1, colour = "black", linetype = "dashed")+
  geom_boxplot(outlier.size = 0.1, width = 0.5)+
  geom_jitter(width = 0.3, height = 0, 
              alpha = 0.6, colour = "black",
              size = 1) + 
  theme_bw()+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Ratio optimal:suboptimal motifs") + xlab("Clade") +
  ylim(-0.02,8)+
  scale_fill_manual(values = mypal, 
                     name = "Clade", guide = guide_legend(),
                     labels = mylabels) +
  scale_x_discrete(labels = mylabels) +
  scale_y_continuous(trans = "log", 
                     breaks = c(1/4,1/2,1,2,4))+
  guides(fill = "none")
p4

png("plots/Figure 8.png", width = 5, height = 5.5, units = 'in', res = 600)
p4
dev.off()

mean(df$cpg[df$host_group == "Carnivore"])
mean(df$cpg[df$host_group == "Bat"])

mean(df$tpa[df$host_group == "Carnivore"])
mean(df$tpa[df$host_group == "Bat"])

t.test(data = df, cpg ~ host_group)
t.test(data = df, tpa ~ host_group)

summary(lm(data = df, ZAP_optimal_motifs ~ cpg_actual))
mean(ZAP_optimal_motifs/cpg_actual)

summary(lm(data= df, ZAP_optimal_motifs ~ cpg_actual))
summary(lm(data= df, ZAP_suboptimal_motifs ~ cpg_actual))

p5 = ggplot(data = df, aes(y = ZAP_optimal_motifs, 
                      x = cpg_actual, colour = clade))+
  geom_jitter(alpha = 0.6,
              size = 2, height = 0.3, width = 0.3) + 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("No. ZAP optimal motifs") + xlab("No. CpG dinucleotides") +
  scale_colour_manual(values = mypal, 
                    name = "Clade", guide = guide_legend(),
                    labels = mylabels) +
  geom_smooth(method = "lm", se = F, col = "black")

p6 = ggplot(data = df, aes(y = ZAP_suboptimal_motifs, 
                      x = cpg_actual, colour = clade))+
  geom_jitter(alpha = 0.6,
              size = 2, height = 0.3, width = 0.3) + 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("No. ZAP suboptimal motifs") + xlab("No. CpG dinucleotides") +
  scale_colour_manual(values = mypal, 
                      name = "Clade", guide = guide_legend(),
                      labels = mylabels) +
  geom_smooth(method = "lm", se = F, col = "black")+
  theme(legend.position = "bottom")

ggpubr::ggarrange(p5,p6, labels = c("A", "B"), 
                  common.legend = T, legend = "bottom")

p7 = ggplot(data = df, aes(y = ZAP_suboptimal_motifs, 
                      x = ZAP_optimal_motifs, colour = clade))+
  geom_jitter(alpha = 0.6,
              size = 2, height = 0.3, width = 0.3) + 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("No. ZAP suboptimal motifs") + 
  xlab("No. ZAP optimal motifs") +
  scale_colour_manual(values = mypal, 
                      name = "Clade", guide = guide_legend(),
                      labels = mylabels) +
  geom_smooth(method = "lm", se = F, col = "black")+
  theme(legend.position = "bottom")

source("code/ZAP CpG locs.R")

ggpubr::ggarrange(p4, g, nrow = 2, labels = "AUTO")
