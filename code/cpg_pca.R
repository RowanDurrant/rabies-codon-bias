library(ggplot2)
library(ggpubr)
library(ggh4x)
library(boot)
library(readxl)

df_pca = read.csv("output_data/PCA_output.csv")
df_cpg = read.csv("output_data/N_CpG.csv")
df_purine = read_xlsx("output_data/Nucleotide_composition_N.xlsx")
colnames(df_purine)[1] = "Accession"

df = merge(df_pca, df_cpg, by.x = c("X", "clade"), by.y = c("accessions", "clade"))
df = merge(df_purine,df, by.x = "Accession", by.y = "X")

df$clade = factor(df$clade, c("Cosmopolitan AF1b", "Cosmopolitan AM2a", "Arctic A", "Asian SEA2a", 
                              "Asian SEA2b", 
                              "Bats DR","Bats TB1", "Bats LC",
                              "Bats EF-E2","RAC-SK SCSK"))

pc = c()
comp = c()
r2 = c()
for(i in c("PC1", "PC2", "PC3")){
  for(j in c("cpg","cpg_actual","gc","tpa","tpa_actual","ta","ZAP_optimal_motifs",
             "ZAP_suboptimal_motifs", colnames(df)[3:59])){
    pc = append(pc, i)
    comp = append(comp,j)
    r2 = append(r2, summary(lm(df[,i]~df[,j]))$adj.r.squared)
  }
}

r_df = as.data.frame(cbind(pc,comp,r2))

set.seed(37856)
tree = read.tree("sequence_data/tree/outgroup_removed.nwk")
ultra.tree = chronos(tree, 0)
Ainv.phylo<-inverseA(ultra.tree,nodes="TIPS")$Ainv

# colnames(df)[colnames(df) == "%G2+A2"] = "GA2"
# m1.phylo = MCMCglmm(PC1~GA2, 
#                     random=~Accession, 
#                     ginverse=list(Accession=Ainv.phylo), 
#                     data=df)
# summary(m1.phylo) #pMCMC = 0.868
# posterior.mode(m1.phylo$VCV[,1]/rowSums(m1.phylo$VCV)) #0.9860341 
# HPDinterval(m1.phylo$VCV[,1]/rowSums(m1.phylo$VCV)) #0.9779891 0.9931691
# 
# colnames(df)[colnames(df) == "%G1+C1"] = "GC1"
# m2.phylo = MCMCglmm(PC2~GC1, 
#                     random=~Accession, 
#                     ginverse=list(Accession=Ainv.phylo), 
#                     data=df)
# summary(m2.phylo) #pMCMC = <0.001
# posterior.mode(m2.phylo$VCV[,1]/rowSums(m2.phylo$VCV)) #0.9874719
# HPDinterval(m2.phylo$VCV[,1]/rowSums(m2.phylo$VCV)) #0.9769578 0.9921081

colnames(df)[colnames(df) == "%G3+T3"] = "GT3"
m3.phylo = MCMCglmm(PC3~GT3, 
                    random=~Accession, 
                    ginverse=list(Accession=Ainv.phylo), 
                    data=df)
summary(m3.phylo) #pMCMC = <0.001
posterior.mode(m3.phylo$VCV[,1]/rowSums(m3.phylo$VCV)) #0.9842323
HPDinterval(m3.phylo$VCV[,1]/rowSums(m3.phylo$VCV)) #0.9710975 0.9911913

my_pal = c("#332288","#88CCEE","#CCDDAA","#44AA99","#117733",  
          "#DDCC77","#999933", "#AA4499","#CC6677","#882255")
my_labels = c("Cosmo AF1b\n(dog)",
             "Cosmo AM2a\n(mongoose)",
             "Arctic A\n(arctic fox)",
             "Asian SEA2a\n(dog)",
             "Asian SEA2b\n(CFB)",
             "Bat DR\n(vampire bat)",
             "Bat TB1\n(Mexican free\n-tailed bat)",
             "Bat LC\n(hoary bat)",
             "Bat EF-E2\n(big brown bat)",
             "RAC-SK SCSK\n(skunk)")

# p1 = ggplot(data = df, aes(x = GA2/100, y = PC1))+
#   geom_smooth(method = "lm", se = F, colour = "black")+
#   geom_point(size = 2, aes(col = clade, shape = clade),
#              alpha = 0.8) +
#   scale_color_manual(values = my_pal, name = "Clade",
#                      labels = my_labels) +
#   scale_shape_manual(name = "Clade",
#                      labels = my_labels,
#                      values = c(17,17,17,17,17,16,16,16,16,17))+
#   theme_bw()+
#   xlab("GA2 content")+
#   annotate("text", x = max(df$GA2)/100, y = -7.5, 
#            label = paste(bquote("R^2 == "),
#                          round(summary(lm(data = df, 
#                                           GA2 ~ PC1))$adj.r.squared,3)), 
#            parse = T,hjust = 1)
# 
# p2 = ggplot(data = df, aes(x = GC1/100, y = PC2))+
#   geom_smooth(method = "lm", se = F, colour = "black")+
#   geom_point(size = 2, aes(col = clade, shape = clade),
#              alpha = 0.8) +
#   scale_color_manual(values = my_pal, name = "Clade",
#                      labels = my_labels) +
#   scale_shape_manual(name = "Clade",
#                      labels = my_labels,
#                      values = c(17,17,17,17,17,16,16,16,16,17))+
#   theme_bw()+
#   xlab("GC1 content")+
#   annotate("text", x = min(df$GC1)/100, y = min(df$PC2), 
#            label = paste(bquote("R^2 == "),round(summary(lm(data = df, GC1 ~ PC2))$adj.r.squared,3)), 
#            parse = T,hjust = 0)


p3 = ggplot(data = df, aes(x = GT3/100, y = PC3))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 0.8) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  theme_bw()+
  xlab("GT3 content")+
  annotate("text", x = min(df$GT3)/100, y = min(df$PC3), 
           label = paste(bquote("R^2 == "),
                         round(summary(lm(data = df, 
                                          GT3 ~ PC3))$adj.r.squared,3)), 
           parse = T,hjust = 0)


g1 = ggplot(data = df, aes(x = PC1, y = PC2))+ 
  geom_point(size = 2, alpha = 0.8, aes(col = clade, shape = clade)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  #ylim(-6,6) + xlim(-8,8)+
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  xlab("PC1 (22.2% explained var.)") + 
  ylab("PC2 (20.4% explained var.)")+
  theme_bw() + 
  coord_axes_inside(labels_inside = TRUE) +
  theme(legend.position = "none")

g2 = ggplot(data = df, aes(x = PC1, y = PC3))+ 
  geom_point(size = 2, alpha = 0.8, aes(col = clade, shape = clade)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  #ylim(-6,6) + xlim(-8,8)+
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  xlab("PC1 (22.2% explained var.)") + 
  ylab("PC3 (16.9% explained var.)")+
  theme_bw() + 
  coord_axes_inside(labels_inside = TRUE) +
  theme(legend.position = "none")

g3 = ggplot(data = df, aes(x = PC2, y = PC3))+ 
  geom_point(size = 2, alpha = 0.8, aes(col = clade, shape = clade)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  #ylim(-6,6) + xlim(-8,8)+
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  xlab("PC2 (20.4% explained var.)") + 
  ylab("PC3 (16.9% explained var.)")+
  theme_bw() + 
  coord_axes_inside(labels_inside = TRUE) +
  theme(legend.position = "none")

ggpubr::ggarrange(g1,g2,g3, p1,p2,p3, common.legend = T, legend = "bottom",
          labels = "AUTO", align = "v")

source("code/ENC-GC3_plot.R")
source("code/loadings_nucleotides.R")


ggpubr::ggarrange(g1,g2,g3, p5,p3,p4, common.legend = T, legend = "bottom",
                  labels = "AUTO", align = "v")

png("plots/Figure 3.png", width = 11, height = 7, units = 'in', res = 600)
ggpubr::ggarrange(g1,g2,g3, p5,p3,p4, common.legend = T, legend = "bottom",
                  labels = "AUTO", align = "v")

dev.off()
