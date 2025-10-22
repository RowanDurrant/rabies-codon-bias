library(ggplot2)
library(ggpubr)

df_pca = read.csv("output_data/PCA_output.csv")
df_cpg = read.csv("output_data/N_CpG.csv")
df_purine = read.csv("output_data/purines.csv")

df = merge(df_pca, df_cpg, by.x = c("Accession", "clade"), by.y = c("accessions", "clade"))
df = merge(df, df_purine, by.x = c("Accession", "clade"), by.y = c("Accession", "Clade"))

df$clade = factor(df$clade, c("Cosmo AF1b", "Cosmo AM2a", "Arctic A", "Asian SEA2a", 
                              "Asian SEA2b", 
                              "Bat DR","Bat TB1", "Bat LC",
                              "Bat EF-E2","RAC-SK SCSK"))


summary(lm(data = df, cpg ~ PC1))$adj.r.squared
summary(lm(data = df, cpg ~ PC2))$adj.r.squared

summary(lm(data = df, tpa ~ PC1))$adj.r.squared
summary(lm(data = df, tpa ~ PC2))$adj.r.squared

df$cpg.c = df$cpg - mean(df$cpg)
df$tpa.c = df$tpa - mean(df$tpa)

summary(lm(data = df, PC2 ~ cpg.c+tpa.c))$adj.r.squared
summary(lm(data = df, PC2 ~ cpg.c*tpa.c))$adj.r.squared

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


p1 = ggplot(data = df, aes(x = cpg, y = PC1))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  theme_bw()+
  xlab("Obs/Exp CpG")+
  annotate("text", x = 0.42, y = -8, label = bquote("R^2 == 0.01023"), 
           parse = T,hjust = 0)

p2 = ggplot(data = df, aes(x = cpg, y = PC2))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  theme_bw()+
  xlab("Obs/Exp CpG")+
  annotate("text", x = 0.42, y = -8, label = bquote("R^2 == 0.3529"), 
           parse = T,hjust = 0)

p5 = ggplot(data = df, aes(x = tpa, y = PC1))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  theme_bw()+
  xlab("Obs/Exp UpA")+
  annotate("text", x = 0.58, y = -8, label = bquote("R^2 == 0.01651"), 
           parse = T,hjust = 0)

p6 = ggplot(data = df, aes(x = tpa, y = PC2))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  theme_bw()+
  xlab("Obs/Exp UpA")+
  annotate("text", x = 0.58, y = -6, label = bquote("R^2 == 0.4998"), 
           parse = T,hjust = 0)

library(ggh4x)
g1 = ggplot(data = df, aes(x = PC1, y = PC2))+ 
  geom_point(size = 2, aes(col = clade, shape = clade)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  xlab("PC1 (24.2% explained var.)") + 
  ylab("PC2 (21.8% explained var.)")+
  theme_bw() + ylim(-10, 10) + xlim(-10, 10)+
  coord_axes_inside(labels_inside = TRUE) +
  theme(legend.position = "none")

g1

summary(lm(data = df, PC1~pyrimidines))
summary(lm(data = df, PC1~pyrimidines1))
summary(lm(data = df, PC1~pyrimidines2))
summary(lm(data = df, PC1~pyrimidines3))

p8 = ggplot(data = df, aes(x = pyrimidines3, y = PC1))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  theme_bw()+
  xlab("Pyrimidine content at 3rd position")+
  annotate("text", x = 0.5, y = -8, label = bquote("R^2 == 0.7293"), 
           parse = T,hjust = 0)


ggarrange(g1, ggarrange(p8,p6, labels = c("B","C"), common.legend = T, legend = "none"), 
          labels = c("A", "B"), nrow = 2, common.legend = T, legend = "bottom")

png("plots/Figure 4.png", width = 8, height = 8, units = 'in', res = 600)
ggarrange(g1, ggarrange(p8,p6, labels = c("B","C"), common.legend = T, legend = "none"), 
          labels = c("A", "B"), nrow = 2, common.legend = T, legend = "bottom")
dev.off()
