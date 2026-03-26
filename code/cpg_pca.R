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


p1 = ggplot(data = df, aes(x = `%G2+A2`, y = PC1))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 0.8) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  theme_bw()+
  xlab("%G2+A2")+
  annotate("text", x = max(df$`%G2+A2`), y = -7.5, 
           label = paste(bquote("R^2 == "),
                         round(summary(lm(data = df, 
                                          `%G2+A2` ~ PC1))$adj.r.squared,3)), 
           parse = T,hjust = 1)

p2 = ggplot(data = df, aes(x = `%G1+C1`, y = PC2))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 0.8) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  theme_bw()+
  xlab("%G1+C1")+
  annotate("text", x = min(df$`%G1+C1`), y = min(df$PC2), 
           label = paste(bquote("R^2 == "),round(summary(lm(data = df, `%G1+C1` ~ PC2))$adj.r.squared,3)), 
           parse = T,hjust = 0)


p3 = ggplot(data = df, aes(x = `%G3+T3`, y = PC3))+
  geom_smooth(method = "lm", se = F, colour = "black")+
  geom_point(size = 2, aes(col = clade, shape = clade),
             alpha = 0.8) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  theme_bw()+
  xlab("%G3+T3")+
  annotate("text", x = min(df$`%G3+T3`), y = min(df$PC3), 
           label = paste(bquote("R^2 == "),
                         round(summary(lm(data = df, 
                                          `%G3+T3` ~ PC3))$adj.r.squared,3)), 
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

png("plots/Figure 3.png", width = 11, height = 7, units = 'in', res = 600)
ggpubr::ggarrange(g1,g2,g3, p1,p2,p3, common.legend = T, legend = "bottom",
                  labels = "AUTO", align = "v")

dev.off()
