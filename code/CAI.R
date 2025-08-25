library("readxl")
library("ggplot2")
library("ggtext")
library(ggpubr)

df = read.csv("output_data/HIVE-CUTS_CAI.csv")
metadata = read.csv("sequence_data/metadata.csv")

df$clade = NA
for(i in 1:nrow(df)){
  df$clade[i] = metadata$Clade[metadata$Accession==df$Name[i]]
}

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

library("tidyr")
df2 <- pivot_longer(df, cols=c(Canis_familiaris, Vulpes_lagopus, Eptesicus_fuscus, Desmodus_rotundus)) 

canis_familiaris_eCAI = 0.785
desmodus_rotundus_eCAI = 0.753
eptesicus_fuscus_eCAI = 0.755
vulpes_lagopus_eCAI = 0.778

df2$normalised = NA
df2$normalised[df2$name == "Canis_familiaris"] = df2$value[df2$name == "Canis_familiaris"]/canis_familiaris_eCAI
df2$normalised[df2$name == "Eptesicus_fuscus"] = df2$value[df2$name == "Eptesicus_fuscus"]/eptesicus_fuscus_eCAI
df2$normalised[df2$name == "Desmodus_rotundus"] = df2$value[df2$name == "Desmodus_rotundus"]/desmodus_rotundus_eCAI
df2$normalised[df2$name == "Vulpes_lagopus"] = df2$value[df2$name == "Vulpes_lagopus"]/vulpes_lagopus_eCAI

p = ggplot(data = df2, aes(x = name, y = normalised, fill = clade))+
  geom_boxplot(alpha = 0.9)+ 
  theme_bw()+ 
scale_fill_manual(values = mypal, labels = mylabels,
                   name = "Clade")+
scale_x_discrete(labels=c(expression(italic("Canis familiaris")), expression(italic("Vulpes lagopus")),
                          expression(italic("Desmodus rotundus")), expression(italic("Eptesicus fuscus"))))+
  ylab("nCAI") +
  xlab("Reference host") +
  guides(fill = "none")

p

source("code/compare cpg cai.R")

p2 = ggarrange(p, p1, common.legend = T, legend = "right")
p2

png("plots/Figure 6.png", width = 12, height = 6, units = 'in', res = 600)
p2
dev.off()
