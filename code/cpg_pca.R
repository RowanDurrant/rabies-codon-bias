library(ggplot2)
library(ggpubr)
library(ggh4x)
library(boot)

df_pca = read.csv("output_data/PCA_output.csv")
df_cpg = read.csv("output_data/N_CpG.csv")
df_purine = read.csv("output_data/purines.csv")

df = merge(df_pca, df_cpg, by.x = c("X", "clade"), by.y = c("accessions", "clade"))
df = merge(df_purine,df, by.x = c("Accession", "Clade"), by.y = c("X", "clade"))

df$clade = factor(df$Clade, c("Cosmopolitan AF1b", "Cosmopolitan AM2a", "Arctic A", "Asian SEA2a", 
                              "Asian SEA2b", 
                              "Bats DR","Bats TB1", "Bats LC",
                              "Bats EF-E2","RAC-SK SCSK"))

summary(lm(data = df, gc ~ PC1))$adj.r.squared #0.223
summary(lm(data = df, gc ~ PC2))$adj.r.squared #0.292
summary(lm(data = df, gc ~ PC3))$adj.r.squared #0

summary(lm(data = df, cpg_actual ~ PC1))$adj.r.squared #0.219
summary(lm(data = df, cpg_actual ~ PC2))$adj.r.squared #0.240
summary(lm(data = df, cpg_actual ~ PC3))$adj.r.squared #0.152

summary(lm(data = df, cpg ~ PC1))$adj.r.squared #0.157
summary(lm(data = df, cpg ~ PC2))$adj.r.squared #0.179
summary(lm(data = df, cpg ~ PC3))$adj.r.squared #0.182

# bootstrap = boot(df,function(data,indices)
#   summary(lm(cpg ~ PC2,data[indices,]))$adj.r.squared,R=10000)
# quantile(bootstrap$t,c(0.025,0.975))

summary(lm(data = df, tpa ~ PC1))$adj.r.squared #0.309
summary(lm(data = df, tpa ~ PC2))$adj.r.squared #0.313
summary(lm(data = df, tpa ~ PC3))$adj.r.squared #0.014

summary(lm(data = df, tpa_actual ~ PC1))$adj.r.squared #0.178
summary(lm(data = df, tpa_actual ~ PC2))$adj.r.squared #0.141
summary(lm(data = df, tpa_actual ~ PC3))$adj.r.squared #0.015

# bootstrap = boot(df,function(data,indices)
#   summary(lm(tpa ~ PC2,data[indices,]))$adj.r.squared,R=10000)
# quantile(bootstrap$t,c(0.025,0.975))

summary(lm(data = df, PC1~purines))$adj.r.squared #0.314
summary(lm(data = df, PC1~purines1))$adj.r.squared #0.170
summary(lm(data = df, PC1~purines2))$adj.r.squared #0.440
summary(lm(data = df, PC1~purines3))$adj.r.squared #0.331

# bootstrap = boot(df,function(data,indices)
#   summary(lm(PC1~purines3,data[indices,]))$adj.r.squared,R=10000)
# quantile(bootstrap$t,c(0.025,0.975))

# df$cpg.c = df$cpg - mean(df$cpg)
# df$tpa.c = df$tpa - mean(df$tpa)
# 
# summary(lm(data = df, PC2 ~ cpg.c+tpa.c))$adj.r.squared
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
  annotate("text", x = 0.58, y = -5, 
           label = paste(bquote("R^2 == "),round(summary(lm(data = df, tpa ~ PC2))$adj.r.squared,3)), 
           parse = T,hjust = 0)


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
  annotate("text", x = 0.465, y = -5, 
           label = paste(bquote("R^2 == "),round(summary(lm(data = df, purines3 ~ PC1))$adj.r.squared,3)), 
           parse = T,hjust = 0)

g1 = ggplot(data = df, aes(x = PC1, y = PC2))+ 
  geom_point(size = 2, alpha = 0.8, aes(col = clade, shape = clade)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  ylim(-5,8) + xlim(-8,8)+
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  xlab("PC1 (24.2% explained var.)") + 
  ylab("PC2 (21.8% explained var.)")+
  theme_bw() + 
  coord_axes_inside(labels_inside = TRUE) +
  theme(legend.position = "none", aspect.ratio = 0.9)

g1

ggpubr::ggarrange(g1, ggpubr::ggarrange(p8,p6, labels = c("B","C"), common.legend = T, legend = "none"), 
          labels = c("A", "B"), nrow = 2, common.legend = T, legend = "bottom",
          heights = c(2,1))

png("plots/Figure 4.png", width = 8, height = 8, units = 'in', res = 600)
ggarrange(g1, ggarrange(p8,p6, labels = c("B","C"), common.legend = T, legend = "none"), 
          labels = c("A", "B"), nrow = 2, common.legend = T, legend = "bottom",
          heights = c(2,1))

dev.off()
