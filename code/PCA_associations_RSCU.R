library(ggplot2)
library(ggpubr)
library(ggh4x)
library(boot)

df_pca = read.csv("output_data/PCA_output_RSCU.csv")
df_cpg = read.csv("output_data/N_CpG.csv")
df_purine = read.csv("output_data/purines.csv")

df = merge(df_pca, df_cpg, by.x = c("Accession", "clade"), by.y = c("accessions", "clade"))
df = merge(df, df_purine, by.x = c("Accession", "clade"), by.y = c("Accession", "Clade"))

df$clade = factor(df$clade, c("Cosmo AF1b", "Cosmo AM2a", "Arctic A", "Asian SEA2a", 
                              "Asian SEA2b", 
                              "Bat DR","Bat TB1", "Bat LC",
                              "Bat EF-E2","RAC-SK SCSK"))

# summary(lm(data = df, cpg ~ PC1))$adj.r.squared
# summary(lm(data = df, cpg ~ PC2))$adj.r.squared
# 
# bootstrap = boot(df,function(data,indices)
#   summary(lm(cpg ~ PC2,data[indices,]))$adj.r.squared,R=10000)
# quantile(bootstrap$t,c(0.025,0.975))

summary(lm(data = df, tpa ~ PC1))$adj.r.squared
summary(lm(data = df, tpa ~ PC2))$adj.r.squared

bootstrap = boot(df,function(data,indices)
  summary(lm(tpa ~ PC2,data[indices,]))$adj.r.squared,R=10000)
quantile(bootstrap$t,c(0.025,0.975))

df$cpg.c = df$cpg - mean(df$cpg)
df$tpa.c = df$tpa - mean(df$tpa)

summary(lm(data = df, PC2 ~ cpg.c+tpa.c))$adj.r.squared
# summary(lm(data = df, PC2 ~ cpg.c*tpa.c))$adj.r.squared
# 
# summary(lm(data = df, PC2 ~ cpg + tpa))$adj.r.squared
# summary(lm(data = df, PC2 ~ cpg*tpa))$adj.r.squared
# 
# bootstrap = boot(df,function(data,indices)
#   summary(lm(PC2 ~ cpg + tpa,data[indices,]))$adj.r.squared,R=10000)
# quantile(bootstrap$t,c(0.025,0.975))

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


p6 = ggplot(data = df, aes(x = tpa, y = PC2))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 0.8) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  theme_bw()+
  xlab("Obs/Exp UpA")+
  annotate("text", x = 0.58, y = -6, label = bquote("R^2 == 0.442"), 
           parse = T,hjust = 0)

g1 = ggplot(data = df, aes(x = PC1, y = PC2))+ 
  geom_point(size = 2, alpha = 0.8, aes(col = clade, shape = clade)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  ylim(-8,8) + xlim(-10,10)+
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  xlab("PC1 (25.0% explained var.)") + 
  ylab("PC2 (20.7% explained var.)")+
  theme_bw() + 
  coord_axes_inside(labels_inside = TRUE) +
  theme(legend.position = "none", aspect.ratio = 0.9)

g1
# 
# summary(lm(data = df, PC1~purines))
# summary(lm(data = df, PC1~purines1))
# summary(lm(data = df, PC1~purines2))
summary(lm(data = df, PC1~purines3))

bootstrap = boot(df,function(data,indices)
  summary(lm(PC1~purines3,data[indices,]))$adj.r.squared,R=10000)
quantile(bootstrap$t,c(0.025,0.975))

p8 = ggplot(data = df, aes(x = purines3, y = PC1))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 0.8) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  theme_bw()+
  xlab("Purine content at 3rd position")+
  annotate("text", x = 0.49, y = -8, label = bquote("R^2 == 0.631"), 
           parse = T,hjust = 0)


ggarrange(g1, ggarrange(p8,p6, labels = c("B","C"), common.legend = T, legend = "none"), 
          labels = c("A", "B"), nrow = 2, common.legend = T, legend = "bottom",
          heights = c(2,1))

# png("plots/Figure 4.png", width = 8, height = 8, units = 'in', res = 600)
# ggarrange(g1, ggarrange(p8,p6, labels = c("B","C"), common.legend = T, legend = "none"), 
#           labels = c("A", "B"), nrow = 2, common.legend = T, legend = "bottom",
#           heights = c(2,1))
# 
# dev.off()
