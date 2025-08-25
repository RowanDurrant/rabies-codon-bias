library(ggplot2)
library(ggpubr)

df_pca = read.csv("output_data/PCA_output_corrected.csv")
df_cpg = read.csv("output_data/N_CpG.csv")
df_vol = read.csv("output_data/sequence_volatility_sp_normalised.csv")
meta = read.csv("sequence_data/metadata.csv")

df = merge(df_pca, df_cpg, by.x = c("Accession", "clade"), by.y = c("accessions", "clade"))
df = merge(df, df_vol, by = "Accession")
df = merge(df, meta, by.x = c("Accession", "clade"), by.y = c("Accession", "Clade"))

summary(lm(data = df, cpg ~ PC1))$adj.r.squared
summary(lm(data = df, cpg ~ PC2))$adj.r.squared

summary(lm(data = df, tpa ~ PC1))$adj.r.squared
summary(lm(data = df, tpa ~ PC2))$adj.r.squared

df$clade = factor(df$clade, c("Cosmo AF1b", "Cosmo AM2a", "Arctic A", "Asian SEA2a", 
                              "Asian SEA2b", 
                              "Bat DR","Bat TB1", "Bat LC",
                              "Bat EF-E2","RAC-SK SCSK"))

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
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlab("Obs/Exp CpG")+
  annotate("text", x = 0.42, y = -8, label = bquote("R^2 == 0.05332057"), 
           parse = T,hjust = 0)+
  ylim(-9,9)

p2 = ggplot(data = df, aes(x = cpg, y = PC2))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlab("Obs/Exp CpG")+
  annotate("text", x = 0.42, y = -8, label = bquote("R^2 == 0.2964098"), 
           parse = T,hjust = 0)+
  ylim(-9,9)

p5 = ggplot(data = df, aes(x = tpa, y = PC1))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlab("Obs/Exp UpA")+
  annotate("text", x = 0.58, y = -8, label = bquote("R^2 == 0.08100637"), 
           parse = T,hjust = 0)+
  ylim(-9,9)

p6 = ggplot(data = df, aes(x = tpa, y = PC2))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlab("Obs/Exp UpA")+
  annotate("text", x = 0.58, y = -8, label = bquote("R^2 == 0.4222357"), 
           parse = T,hjust = 0)+
  ylim(-9,9)

ggarrange(p1,p2,p5,p6, common.legend = T, legend = "bottom")

summary(lm(data = df, Hamming~ PC1))$adj.r.squared
summary(lm(data = df, Hamming ~ PC2))$adj.r.squared
summary(lm(data = df, Miyata ~ PC1))$adj.r.squared
summary(lm(data = df, Miyata ~ PC2))$adj.r.squared
summary(lm(data = df, Hamming_stop0~ PC1))$adj.r.squared
summary(lm(data = df, Hamming_stop0 ~ PC2))$adj.r.squared
summary(lm(data = df, Miyata_stop0 ~ PC1))$adj.r.squared
summary(lm(data = df, Miyata_stop0 ~ PC2))$adj.r.squared
summary(lm(data = df, Hamming_norm ~ PC1))$adj.r.squared
summary(lm(data = df, Hamming_norm ~ PC2))$adj.r.squared
summary(lm(data = df, Miyata_norm ~ PC1))$adj.r.squared
summary(lm(data = df, Miyata_norm ~ PC2))$adj.r.squared
summary(lm(data = df, Hamming_stop0_norm~ PC1))$adj.r.squared
summary(lm(data = df, Hamming_stop0_norm ~ PC2))$adj.r.squared
summary(lm(data = df, Miyata_stop0_norm ~ PC1))$adj.r.squared
summary(lm(data = df, Miyata_stop0_norm ~ PC2))$adj.r.squared

p7 = ggplot(data = df, aes(x = Hamming, y = PC1))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlab("Mean Hamming Volatility")+
  annotate("text", x = min(df$Hamming), y = -8, label = bquote("R^2 == 0.05895453"), 
           parse = T,hjust = 0)+
  ylim(-9,9)

p8 = ggplot(data = df, aes(x = Hamming, y = PC2))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlab("Mean Hamming Volatility")+
  annotate("text", x = min(df$Hamming), y = -8, label = bquote("R^2 == 0.4798058"), 
           parse = T,hjust = 0)+
  ylim(-9,9)
p9 = ggplot(data = df, aes(x = Miyata, y = PC1))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlab("Mean Miyata Volatility")+
  annotate("text", x = min(df$Miyata), y = -8, label = bquote("R^2 == -0.00224301"), 
           parse = T,hjust = 0)+
  ylim(-9,9)

p10 = ggplot(data = df, aes(x = Miyata, y = PC2))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlab("Mean Miyata Volatility")+
  annotate("text", x = min(df$Miyata), y = -8, label = bquote("R^2 == 0.2745893"), 
           parse = T,hjust = 0)+
  ylim(-9,9)

p11 = ggplot(data = df, aes(x = Hamming_stop0, y = PC1))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlab("Mean Hamming Volatility\n(stop = 0)")+
  annotate("text", x = min(df$Hamming_stop0), y = -8, label = bquote("R^2 == 0.161566"), 
           parse = T,hjust = 0)+
  ylim(-9,9)

p12 = ggplot(data = df, aes(x = Hamming_stop0, y = PC2))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlab("Mean Hamming Volatility\n(stop = 0)")+
  annotate("text", x = min(df$Hamming_stop0), y = -8, label = bquote("R^2 == 0.1641607"), 
           parse = T,hjust = 0)+
  ylim(-9,9)
p13 = ggplot(data = df, aes(x = Miyata_stop0, y = PC1))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlab("Mean Miyata Volatility\n(stop = 0)")+
  annotate("text", x = min(df$Miyata_stop0), y = -8, label = bquote("R^2 == 0.3380395"), 
           parse = T,hjust = 0)+
  ylim(-9,9)

p14 = ggplot(data = df, aes(x = Miyata_stop0, y = PC2))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlab("Mean Miyata Volatility\n(stop = 0)")+
  annotate("text", x = min(df$Miyata_stop0), y = -8, label = bquote("R^2 == 0.2026575"), 
           parse = T,hjust = 0)+
  ylim(-9,9)

p15 = ggplot(data = df, aes(x = Hamming_norm, y = PC1))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlab("Normalised Hamming Volatility")+
  annotate("text", x = min(df$Hamming_norm), y = -8, label = bquote("R^2 == 0.03867999"), 
           parse = T,hjust = 0)+
  ylim(-9,9)

p16 = ggplot(data = df, aes(x = Hamming_norm, y = PC2))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlab("Normalised Hamming Volatility")+
  annotate("text", x = min(df$Hamming_norm), y = -8, label = bquote("R^2 == 0.3479759"), 
           parse = T,hjust = 0)+
  ylim(-9,9)
p17 = ggplot(data = df, aes(x = Miyata_norm, y = PC1))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlab("Normalised Miyata Volatility")+
  annotate("text", x = min(df$Miyata_norm), y = -8, label = bquote("R^2 == -0.001898993"), 
           parse = T,hjust = 0)+
  ylim(-9,9)

p18 = ggplot(data = df, aes(x = Miyata_norm, y = PC2))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlab("Normalised Miyata Volatility")+
  annotate("text", x = min(df$Miyata_norm), y = -8, label = bquote("R^2 == 0.3316058"), 
           parse = T,hjust = 0)+
  ylim(-9,9)

p19 = ggplot(data = df, aes(x = Hamming_stop0_norm, y = PC1))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlab("Normalised Hamming Volatility\n(stop = 0)")+
  annotate("text", x = min(df$Hamming_stop0_norm), y = -8, label = bquote("R^2 == 0.002885119"), 
           parse = T,hjust = 0)+
  ylim(-9,9)

p20 = ggplot(data = df, aes(x = Hamming_stop0_norm, y = PC2))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlab("Normalised Hamming Volatility\n(stop = 0)")+
  annotate("text", x = min(df$Hamming_stop0_norm), y = -8, label = bquote("R^2 == 0.1970466"), 
           parse = T,hjust = 0)+
  ylim(-9,9)
p21 = ggplot(data = df, aes(x = Miyata_stop0_norm, y = PC1))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlab("Normalised Miyata Volatility (stop = 0)")+
  annotate("text", x = min(df$Miyata_stop0_norm), y = -8, label = bquote("R^2 == 0.2857238"), 
           parse = T,hjust = 0)+
  ylim(-9,9)

p22 = ggplot(data = df, aes(x = Miyata_stop0_norm, y = PC2))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 1) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlab("Normalised Miyata Volatility\n(stop = 0)")+
  annotate("text", x = min(df$Miyata_stop0_norm), y = -8, label = bquote("R^2 == 0.002129578"), 
           parse = T,hjust = 0)+
  ylim(-9,9)

ggarrange(p7,p8,p9,p10,
          p11, p12, p13, p14,
          p15, p16, p17, p18,
          p19,p20,p21,p22,
          nrow = 4, ncol = 4,
         common.legend = T, legend = "bottom")

library(ggh4x)
g1 = ggplot(data = df, aes(x = PC1, y = PC2))+ 
  geom_point(size = 2, aes(col = clade, shape = clade)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  xlab("PC1 (24.0% explained var.)") + 
  ylab("PC2 (21.7% explained var.)")+
  theme_bw() + ylim(-10, 10) + xlim(-10, 10)+
  coord_axes_inside(labels_inside = TRUE) +
  theme(legend.position = "none")

g1

ggarrange(g1, ggarrange(p6,p8, common.legend = T, legend = "bottom"), nrow = 2)

png("plots/Figure 5.png", width = 7, height = 7.5, units = 'in', res = 600)
ggarrange(g1, ggarrange(p6,p8, common.legend = T, legend = "bottom"), nrow = 2)
dev.off()
  