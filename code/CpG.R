##CpG 
# Observed to Expected CpG is calculated as below : 
# Obs/Exp CpG = Number of CpG * N / (Number of C * Number of G) 
# where N = length of sequence.

library(stringr)
library(ggpubr)
library(seqinr)
library(ggplot2)
library(egg)
library(ggimage)
library(tibble)

CpG = function(x){
  x = unname(as.character(x))
  nG = str_count(x, "g")
  nC = str_count(x, "c")
  nCpG = str_count(x, "cg")
  N = nchar(x)
  ObsExpCpG = (nCpG/(nC * nG))*N
  return(ObsExpCpG)
}

GCcontent = function(x){
  x = unname(as.character(x))
  nG = str_count(x, "g")
  nC = str_count(x, "c")
  N = nchar(x)
  gccont = (nC + nG)/N
  return(gccont)
  
}

CpG_actual = function(x){
  x = unname(as.character(x))
  nCpG = str_count(x, "cg")
  return(nCpG)
}

ZAP_optimal = function(x){
  x = unname(as.character(x))
  n_optimal = str_count(x, pattern = "c.......g.cg")
  n_optimal = n_optimal + str_count(x, pattern = "c....g.cg")
  n_optimal = n_optimal + str_count(x, pattern = "c.....g.cg")
  n_optimal = n_optimal + str_count(x, pattern = "c......g.cg")
  n_optimal = n_optimal + str_count(x, pattern = "c........g.cg")
  return(n_optimal)
}

ZAP_suboptimal = function(x){
  x = unname(as.character(x))
  n_optimal = str_count(x, pattern = "c.......c.cg")
  n_optimal = n_optimal + str_count(x, pattern = "c....c.cg")
  n_optimal = n_optimal + str_count(x, pattern = "c.....c.cg")
  n_optimal = n_optimal + str_count(x, pattern = "c......c.cg")
  n_optimal = n_optimal + str_count(x, pattern = "c........c.cg")
  return(n_optimal)
}

TpA = function(x){
  x = unname(as.character(x))
  nT = str_count(x, "t")
  nA = str_count(x, "a")
  nTpA = str_count(x, "ta")
  N = nchar(x)
  ObsExpTpA = (nTpA/(nT * nA))*N
  return(ObsExpTpA)
}

TAcontent = function(x){
  x = unname(as.character(x))
  nT = str_count(x, "t")
  nA = str_count(x, "a")
  N = nchar(x)
  tacont = (nT + nA)/N
  return(tacont)
  
}

TpA_actual = function(x){
  x = unname(as.character(x))
  nTpA = str_count(x, "ta")
  return(nTpA)
}

seqs = seqinr::read.fasta("sequence_data/all_seqs.fasta",as.string = T)
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
       cpg = append(cpg, CpG(seqs[[j]][1]))
       tpa = append(tpa, TpA(seqs[[j]][1]))
       gc = append(gc, GCcontent(seqs[[j]][1]))
       ta = append(ta, TAcontent(seqs[[j]][1]))
       cpg_actual = append(cpg_actual, CpG_actual(seqs[[j]][1]))
       tpa_actual= append(tpa_actual, TpA_actual(seqs[[j]][1]))
       accessions = append(accessions, names(seqs)[j])
       clade = append(clade, metadata$Clade[metadata$Accession==names(seqs)[j]])
       #host_group = append(host_group, metadata$Group[metadata$Accession==names(seqs)[j]])
       ZAP_optimal_motifs = append(ZAP_optimal_motifs, ZAP_optimal(seqs[[j]][1]))
       ZAP_suboptimal_motifs = append(ZAP_suboptimal_motifs, ZAP_suboptimal(seqs[[j]][1]))
    }
  
df = data.frame(accessions, clade, cpg, gc, cpg_actual, 
                ZAP_optimal_motifs, ZAP_suboptimal_motifs, tpa, ta, tpa_actual)
write.csv(df, "output_data/N_CpG.csv")
df$clade = factor(df$clade, c("Cosmopolitan AF1b", "Cosmopolitan AM2a", "Arctic A", "Asian SEA2a", 
                            "Asian SEA2b", 
                            "Bats DR","Bats TB1", "Bats LC",
                           "Bats EF-E2","RAC-SK SCSK"))

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
                    labels = mylabels)+
  geom_image(
    data = tibble(clade = "Cosmopolitan AF1b", cpg = 0.67),
    aes(image = "phylopic_images/dog.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Cosmopolitan AM2a", cpg = 0.67),
    aes(image = "phylopic_images/mongoose.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Asian SEA2a", cpg = 0.67),
    aes(image = "phylopic_images/dog.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Asian SEA2b", cpg = 0.67),
    aes(image = "phylopic_images/cfb.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Arctic A", cpg = 0.67),
    aes(image = "phylopic_images/fox.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Bats LC", cpg = 0.67),
    aes(image = "phylopic_images/lasiurus.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Bats DR", cpg = 0.67),
    aes(image = "phylopic_images/desmodus.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Bats TB1", cpg = 0.67),
    aes(image = "phylopic_images/tadarida.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Bats EF-E2", cpg = 0.67),
    aes(image = "phylopic_images/eptesicus.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "RAC-SK SCSK", cpg = 0.67),
    aes(image = "phylopic_images/skunk.png"),
    size = 0.08)

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
  scale_x_discrete(labels = mylabels)+
  geom_image(
    data = tibble(clade = "Cosmopolitan AF1b", cpg_actual = 44),
    aes(image = "phylopic_images/dog.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Cosmopolitan AM2a", cpg_actual = 44),
    aes(image = "phylopic_images/mongoose.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Asian SEA2a", cpg_actual = 44),
    aes(image = "phylopic_images/dog.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Asian SEA2b", cpg_actual = 44),
    aes(image = "phylopic_images/cfb.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Arctic A", cpg_actual = 44),
    aes(image = "phylopic_images/fox.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Bats LC", cpg_actual = 44),
    aes(image = "phylopic_images/lasiurus.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Bats DR", cpg_actual = 44),
    aes(image = "phylopic_images/desmodus.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Bats TB1", cpg_actual = 44),
    aes(image = "phylopic_images/tadarida.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Bats EF-E2", cpg_actual = 44),
    aes(image = "phylopic_images/eptesicus.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "RAC-SK SCSK", cpg_actual = 44),
    aes(image = "phylopic_images/skunk.png"),
    size = 0.08)

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
                    labels = mylabels)+
  geom_image(
    data = tibble(clade = "Cosmopolitan AF1b", tpa = 0.73),
    aes(image = "phylopic_images/dog.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Cosmopolitan AM2a", tpa = 0.73),
    aes(image = "phylopic_images/mongoose.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Asian SEA2a", tpa = 0.73),
    aes(image = "phylopic_images/dog.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Asian SEA2b", tpa = 0.73),
    aes(image = "phylopic_images/cfb.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Arctic A", tpa = 0.73),
    aes(image = "phylopic_images/fox.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Bats LC", tpa = 0.73),
    aes(image = "phylopic_images/lasiurus.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Bats DR", tpa = 0.73),
    aes(image = "phylopic_images/desmodus.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Bats TB1", tpa = 0.73),
    aes(image = "phylopic_images/tadarida.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Bats EF-E2", tpa = 0.73),
    aes(image = "phylopic_images/eptesicus.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "RAC-SK SCSK", tpa = 0.73),
    aes(image = "phylopic_images/skunk.png"),
    size = 0.08)


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
  scale_x_discrete(labels = mylabels)+
  geom_image(
    data = tibble(clade = "Cosmopolitan AF1b", tpa_actual = 76),
    aes(image = "phylopic_images/dog.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Cosmopolitan AM2a", tpa_actual = 76),
    aes(image = "phylopic_images/mongoose.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Asian SEA2a", tpa_actual = 76),
    aes(image = "phylopic_images/dog.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Asian SEA2b", tpa_actual = 76),
    aes(image = "phylopic_images/cfb.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Arctic A", tpa_actual = 76),
    aes(image = "phylopic_images/fox.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Bats LC", tpa_actual = 76),
    aes(image = "phylopic_images/lasiurus.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Bats DR", tpa_actual = 76),
    aes(image = "phylopic_images/desmodus.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Bats TB1", tpa_actual = 76),
    aes(image = "phylopic_images/tadarida.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "Bats EF-E2", tpa_actual = 76),
    aes(image = "phylopic_images/eptesicus.png"),
    size = 0.08)+
  geom_image(
    data = tibble(clade = "RAC-SK SCSK", tpa_actual = 76),
    aes(image = "phylopic_images/skunk.png"),
    size = 0.08)

egg::ggarrange(p1,p6,p3,p8, labels = c("A", "B", "C", "D"))

png("plots/Figure 5.png", width = 8, height = 8, units = 'in', res = 600)
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
  scale_fill_manual(values = mypal, 
                     name = "Clade", guide = guide_legend(),
                     labels = mylabels) +
  scale_x_discrete(labels = mylabels) +
  scale_y_continuous(trans = "log", 
                     breaks = c(1/4,1/2,1,2,4, 8),
                     limits = c(-0.02,10))+
  guides(fill = "none")+
  geom_image(
    data = tibble(clade = "Cosmopolitan AF1b", ratio = 8),
    aes(image = "phylopic_images/dog.png"),
    size = 0.15)+
  geom_image(
    data = tibble(clade = "Cosmopolitan AM2a", ratio = 8),
    aes(image = "phylopic_images/mongoose.png"),
    size = 0.15)+
  geom_image(
    data = tibble(clade = "Asian SEA2a", ratio = 8),
    aes(image = "phylopic_images/dog.png"),
    size = 0.15)+
  geom_image(
    data = tibble(clade = "Asian SEA2b", ratio = 8),
    aes(image = "phylopic_images/cfb.png"),
    size = 0.15)+
  geom_image(
    data = tibble(clade = "Arctic A", ratio = 8),
    aes(image = "phylopic_images/fox.png"),
    size = 0.15)+
  geom_image(
    data = tibble(clade = "Bats LC", ratio = 8),
    aes(image = "phylopic_images/lasiurus.png"),
    size = 0.15)+
  geom_image(
    data = tibble(clade = "Bats DR", ratio = 8),
    aes(image = "phylopic_images/desmodus.png"),
    size = 0.15)+
  geom_image(
    data = tibble(clade = "Bats TB1", ratio = 8),
    aes(image = "phylopic_images/tadarida.png"),
    size = 0.15)+
  geom_image(
    data = tibble(clade = "Bats EF-E2", ratio = 8),
    aes(image = "phylopic_images/eptesicus.png"),
    size = 0.15)+
  geom_image(
    data = tibble(clade = "RAC-SK SCSK", ratio = 8),
    aes(image = "phylopic_images/skunk.png"),
    size = 0.15)
p4

library(DescTools)
MeanCI(df$cpg[df$host_group == "Carnivore"], conf.level = 0.95)
MeanCI(df$cpg[df$host_group == "Bat"], conf.level = 0.95)

MeanCI(df$tpa[df$host_group == "Carnivore"], conf.level = 0.95)
MeanCI(df$tpa[df$host_group == "Bat"], conf.level = 0.95)
# 
# t.test(data = df, cpg ~ host_group)
# t.test(data = df, tpa ~ host_group)
# 
# summary(lm(data = df, ZAP_optimal_motifs ~ cpg_actual))
# mean(ZAP_optimal_motifs/cpg_actual)
# 
# summary(lm(data= df, ZAP_optimal_motifs ~ cpg_actual))
# summary(lm(data= df, ZAP_suboptimal_motifs ~ cpg_actual))
# 
# p5 = ggplot(data = df, aes(y = ZAP_optimal_motifs, 
#                       x = cpg_actual, colour = clade))+
#   geom_jitter(alpha = 0.6,
#               size = 2, height = 0.3, width = 0.3) + 
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ylab("No. ZAP optimal motifs") + xlab("No. CpG dinucleotides") +
#   scale_colour_manual(values = mypal, 
#                     name = "Clade", guide = guide_legend(),
#                     labels = mylabels) +
#   geom_smooth(method = "lm", se = F, col = "black")
# 
# p6 = ggplot(data = df, aes(y = ZAP_suboptimal_motifs, 
#                       x = cpg_actual, colour = clade))+
#   geom_jitter(alpha = 0.6,
#               size = 2, height = 0.3, width = 0.3) + 
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ylab("No. ZAP suboptimal motifs") + xlab("No. CpG dinucleotides") +
#   scale_colour_manual(values = mypal, 
#                       name = "Clade", guide = guide_legend(),
#                       labels = mylabels) +
#   geom_smooth(method = "lm", se = F, col = "black")+
#   theme(legend.position = "bottom")
# 
# ggpubr::ggarrange(p5,p6, labels = c("A", "B"), 
#                   common.legend = T, legend = "bottom")
# 
# p7 = ggplot(data = df, aes(y = ZAP_suboptimal_motifs, 
#                       x = ZAP_optimal_motifs, colour = clade))+
#   geom_jitter(alpha = 0.6,
#               size = 2, height = 0.3, width = 0.3) + 
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ylab("No. ZAP suboptimal motifs") + 
#   xlab("No. ZAP optimal motifs") +
#   scale_colour_manual(values = mypal, 
#                       name = "Clade", guide = guide_legend(),
#                       labels = mylabels) +
#   geom_smooth(method = "lm", se = F, col = "black")+
#   theme(legend.position = "bottom")

source("code/ZAP_CpG_locs.R")
ggpubr::ggarrange(p4, g, nrow = 2, labels = "AUTO")

png("plots/Figure 7.png", width = 10, height = 7.5, units = 'in', res = 600)
ggpubr::ggarrange(p4, g, nrow = 2, labels = "AUTO")
dev.off()
