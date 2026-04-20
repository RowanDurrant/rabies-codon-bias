set.seed(84298913)
library(ape)
library(ggplot2)
library(ggpubr)
library(ggh4x)
library(boot)
library(readxl)
library(MCMCglmm)

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
cc = c()
for(i in c("PC1", "PC2", "PC3")){
  for(j in c("cpg","cpg_actual","gc","tpa","tpa_actual","ta","ZAP_optimal_motifs",
             "ZAP_suboptimal_motifs", colnames(df)[3:43])){
    pc = append(pc, i)
    comp = append(comp,j)
    cc = append(cc, cor(df[,j]/100, df[,i]))
  }
}

r_df = as.data.frame(cbind(pc,comp,cc))
r_df$cc_abs = abs(as.numeric(r_df$cc))

r_df1 = r_df[r_df$pc == "PC1",]
r_df2 = r_df[r_df$pc == "PC2",]
r_df3 = r_df[r_df$pc == "PC3",]
top2 = rbind(r_df1[rev(order(r_df1$cc_abs)),][1:4,],
             r_df2[rev(order(r_df2$cc_abs)),][1:4,],
             r_df3[rev(order(r_df3$cc_abs)),][1:4,])

top2$upper = NA
top2$lower = NA
top2$label = NA

library(confintr)
for(i in 1:nrow(top2)){
  
  top2$lower[i] =ci_cor(y = df[,top2$pc[i]], x = df[,top2$comp[i]])$interval[1]
  top2$upper[i] =ci_cor(y = df[,top2$pc[i]], x = df[,top2$comp[i]])$interval[2]
  top2$label[i] = paste(top2$pc[i], "~", top2$comp[i])
}

top2$label = factor(top2$label, levels = rev(top2$label))

p4 = ggplot(data = top2, aes(y = as.numeric(cc), x = label))+
  geom_point(colour = "blue", size = 0.5)+
  geom_errorbar(aes(ymin=as.numeric(lower), ymax=as.numeric(upper)),
                width = 0.1)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_hline(yintercept = c(1,-1))+
  scale_x_discrete(labels = c("PC3~AC", "PC3~GT",
                              "PC3~GT3", "PC3~AC3",
                              "PC2~GA", "PC2~CT",
                              "PC2~GA3", "PC2~CT3",
                              "PC1~UpA", "PC1~AT1",
                              "PC1~GC1", 'PC1~C1'))+
  theme_bw()+ ylab("Correlation coefficient") + xlab("")+
  coord_flip()

tree = read.tree("sequence_data/all_outgroups/no_outgroup.nwk")
ultra.tree = chronos(tree, 0)
Ainv.phylo<-inverseA(ultra.tree,nodes="TIPS")$Ainv

prior1 <- list(
  G = list(G1 = list(V = 1, nu = 0.002)),
  R = list(V = 1, nu = 0.002)
)

colnames(df)[colnames(df) == "%C1"] = "pC1"
m1.phylo = MCMCglmm(PC1~pC1,
                    random=~Accession,
                    ginverse=list(Accession=Ainv.phylo),
                    data=df, prior=prior1)
plot(m1.phylo$Sol)
plot(m1.phylo$VCV)
summary(m1.phylo) #pMCMC = <0.001
posterior.mode(m1.phylo$VCV[,1]/rowSums(m1.phylo$VCV)) #0.9879975
HPDinterval(m1.phylo$VCV[,1]/rowSums(m1.phylo$VCV)) #0.9797923 0.9935779

colnames(df)[colnames(df) == "%G3+A3"] = "GA3"
m2.phylo = MCMCglmm(PC2~GA3,
                    random=~Accession,
                    ginverse=list(Accession=Ainv.phylo),
                    data=df, prior = prior1)
plot(m2.phylo$Sol)
plot(m2.phylo$VCV)
summary(m2.phylo) #pMCMC = <0.001
posterior.mode(m2.phylo$VCV[,1]/rowSums(m2.phylo$VCV)) #0.9825264 
HPDinterval(m2.phylo$VCV[,1]/rowSums(m2.phylo$VCV)) #0.9717638 0.9883312

colnames(df)[colnames(df) == "%G3+T3"] = "GT3"
m3.phylo = MCMCglmm(PC3~GT3, 
                    random=~Accession, 
                    ginverse=list(Accession=Ainv.phylo), 
                    data=df, prior = prior1)
plot(m3.phylo$Sol)
plot(m3.phylo$VCV)
summary(m3.phylo) #pMCMC = <0.001
posterior.mode(m3.phylo$VCV[,1]/rowSums(m3.phylo$VCV)) #0.9893406
HPDinterval(m3.phylo$VCV[,1]/rowSums(m3.phylo$VCV)) #0.9794927 0.9945502

summary(lm(data = df, PC1~pC1))$adj.r.squared #0.6656295
bootstrap = boot(df,function(data,indices)
  summary(lm(PC1~pC1,data[indices,]))$adj.r.squared,R=10000)
quantile(bootstrap$t,c(0.025,0.975)) #0.6200261 0.7084308 

summary(lm(data = df, PC2~GA3))$adj.r.squared #0.7289015
bootstrap = boot(df,function(data,indices)
  summary(lm(PC2~GA3,data[indices,]))$adj.r.squared,R=10000)
quantile(bootstrap$t,c(0.025,0.975)) #0.6878962 0.7678734 

summary(lm(data = df, PC3~GT3))$adj.r.squared #0.8399556
bootstrap = boot(df,function(data,indices)
  summary(lm(PC3~GT3,data[indices,]))$adj.r.squared,R=10000)
quantile(bootstrap$t,c(0.025,0.975)) #0.8141142 0.8640103


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

p1 = ggplot(data = df, aes(x = pC1/100, y = PC1))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 0.8) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  theme_bw()+
  xlab("C1 content")+
  annotate("text", x = max(df$pC1)/100, y = -7.5,
           label = paste(bquote("R^2 == "),
                         round(summary(lm(data = df,
                                          pC1 ~ PC1))$adj.r.squared,3)),
           parse = T,hjust = 1)

p2 = ggplot(data = df, aes(x = GA3/100, y = PC2))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 0.8) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  theme_bw()+
  xlab("GA3 content")+
  annotate("text", x = min(df$GA3)/100, y = min(df$PC2),
           label = paste(bquote("R^2 == "),round(summary(lm(data = df, GA3 ~ PC2))$adj.r.squared,3)),
           parse = T,hjust = 0)

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
  xlab("PC1 (22.4% explained var.)") + 
  ylab("PC2 (20.6% explained var.)")+
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
  xlab("PC1 (22.4% explained var.)") + 
  ylab("PC3 (16.5% explained var.)")+
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
  xlab("PC2 (20.6% explained var.)") + 
  ylab("PC3 (16.5% explained var.)")+
  theme_bw() + 
  coord_axes_inside(labels_inside = TRUE) +
  theme(legend.position = "none")

ggpubr::ggarrange(g1,g2,g3, p1,p2,p3, common.legend = T, legend = "bottom",
          labels = "AUTO", align = "v")

source("code/ENC-GC3_plot.R")

ggpubr::ggarrange(ggpubr::ggarrange(g1,g2,g3,
                                    nrow = 3,
                                    labels = "AUTO"), 
                  ggpubr::ggarrange(p1,p2,p3, nrow = 3, 
                                    legend = "none", 
                                    labels = c("D", "E", "F")),
                  ggpubr::ggarrange(p4,p5, nrow = 2,
                                    
                                    labels = c("G", "H")),
                  
                  ncol = 3,
                  common.legend = T, legend = "bottom")

png("plots/Figure 3.png", width = 12, height = 8, units = 'in', res = 600)
ggpubr::ggarrange(ggpubr::ggarrange(g1,g2,g3,
                                    nrow = 3,
                                    labels = "AUTO"), 
                  ggpubr::ggarrange(p1,p2,p3, nrow = 3, 
                                    legend = "none", 
                                    labels = c("D", "E", "F")),
                  ggpubr::ggarrange(p4,p5, nrow = 2,
                                    labels = c("G", "H")),
                  
                  ncol = 3,
                  common.legend = T, legend = "bottom")

dev.off()
