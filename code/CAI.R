library("readxl")
library("ggplot2")
library("ggtext")

df = read.csv("output_data/HIVE-CUTS_CAI.csv")
metadata = read.csv("sequence_data/metadata.csv")

df$clade = NA
for(i in 1:nrow(df)){
  df$clade[i] = metadata$Clade[metadata$Accession==df$Name[i]]
}

df$clade = factor(df$clade, c("Cosmo AF1b", "Cosmo AM2a", "Arctic A", "Asian SEA2a", 
                                "Asian SEA2b", 
                                "Bat TB1",
                                "Bat DR", "Bat EF-E2","RAC-SK SCSK", "Bat LC"))

my_pal <- c("#332288","#88CCEE","#CCDDAA","#44AA99","#117733",  
            "#999933", "#DDCC77","#CC6677","#882255","#AA4499")

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
scale_fill_manual(values = my_pal, labels = c("Cosmo AF1b\n(dog)",
                                               "Cosmo AM2a\n(mongoose)",
                                              "Arctic A\n(arctic fox)",
                                               "Asian SEA2a\n(dog)",
                                               "Asian SEA2b\n(CFB)",
                                               "Bat TB1\n(Mexican free\n-tailed bat)",
                                               "Bat DR\n(vampire bat)",
                                               "Bat EF-E2\n(big brown bat)",
                                               "RAC-SK SCSK\n(skunk)",
                                               "Bat LC\n(hoary bat)"),
                   name = "Clade")+
scale_x_discrete(labels=c("Canis familiaris", "Vulpes lagopus",
                                               "Desmodus rotundus", "Eptesicus fuscus"))+
  ylab("Normalised codon adapatation index") +
  xlab("Reference host") +
  theme(legend.position = "bottom")

p

png("plots/Figure 6.png", width = 7.5, height = 6, units = 'in', res = 600)
p
dev.off()
