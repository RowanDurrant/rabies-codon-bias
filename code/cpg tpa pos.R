library(ggplot2)
library(ggpubr)

df = read.csv("output_data/myseqs_RSDUc_all.tsv", sep='\t', row.names = 1, header= TRUE)
df$Accessions = rownames(df)
meta = read.csv("sequence_data/metadata.csv")

df = merge(df, meta, by.x = "Accessions", by.y = "Accession")
df$Clade = factor(df$Clade, c("Cosmo AF1b", "Cosmo AM2a", "Arctic A", "Asian SEA2a", 
                              "Asian SEA2b", 
                              "Bat TB1",
                              "Bat DR", "Bat EF-E2","RAC-SK SCSK", "Bat LC"))

my_pal <- c("#332288","#88CCEE","#CCDDAA","#44AA99","#117733",  
            "#999933", "#DDCC77","#CC6677","#882255","#AA4499")
my_labels = c("Cosmo AF1b\n(dog)",
              "Cosmo AM2a\n(mongoose)",
              "Arctic A\n(arctic fox)",
              "Asian SEA2a\n(dog)",
              "Asian SEA2b\n(CFB)",
              "Bat TB1\n(Mexican free\n-tailed bat)",
              "Bat DR\n(vampire bat)",
              "Bat EF-E2\n(big brown bat)",
              "RAC-SK SCSK\n(skunk)",
              "Bat LC\n(hoary bat)")

p1 = ggplot(data = df, aes(x = CpGpos1, y = CpGpos2))+
  geom_abline(slope = 1, linetype = "dashed", colour = "grey")+
  geom_point(size = 2, aes(col = Clade, shape = Clade),
             alpha = 1)+
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlim(0,0.4) + ylim(0,0.4)+
  geom_smooth(se = F, method = "lm", col = "black")

p2 = ggplot(data = df, aes(x = CpGpos1, y = CpGbridge))+
  geom_abline(slope = 1, linetype = "dashed", colour = "grey")+
  geom_point(size = 2, aes(col = Clade, shape = Clade),
             alpha = 1)+
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlim(0,0.4) + ylim(0,0.4)+
  geom_smooth(se = F, method = "lm", col = "black")

p3 = ggplot(data = df, aes(x = CpGpos2, y = CpGbridge))+
  geom_abline(slope = 1, linetype = "dashed", colour = "grey")+
  geom_point(size = 2, aes(col = Clade, shape = Clade),
             alpha = 1)+
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlim(0,0.4) + ylim(0,0.4)+
  geom_smooth(se = F, method = "lm", col = "black")

p4 = ggplot(data = df, aes(x = TpApos2, y = TpAbridge))+
  geom_abline(slope = 1, linetype = "dashed", colour = "grey")+
  geom_point(size = 2, aes(col = Clade, shape = Clade),
             alpha = 1)+
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  theme_bw()+
  xlim(0,0.4) + ylim(0,0.4)+
  geom_smooth(se = F, method = "lm", col = "black")

ggarrange(p1,p2,p3,p4, common.legend = T, legend = "bottom")


cor.test(df$CpGpos1, df$CpGpos2) #0.2196182  (0.1277046, 0.3077918)
cor.test(df$CpGpos1, df$CpGbridge) #0.1999115 (0.1073755, 0.2890136)
cor.test(df$CpGbridge, df$CpGpos2) #-0.125178 (-0.21717309, -0.03097882)
cor.test(df$TpApos2, df$TpAbridge) #-0.139859 (-0.23136491 -0.04590021)
