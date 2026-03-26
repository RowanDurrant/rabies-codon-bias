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

df$clade = factor(df$clade, c("Cosmopolitan AF1b", "Cosmopolitan AM2a", "Arctic A", "Asian SEA2a", 
                              "Asian SEA2b", 
                              "Bats DR","Bats TB1", "Bats LC",
                              "Bats EF-E2","RAC-SK SCSK"))

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

canis_familiaris_eCAI = 0.784
desmodus_rotundus_eCAI = 0.754
eptesicus_fuscus_eCAI = 0.753
vulpes_lagopus_eCAI = 0.774

df$normalised = NA
df$normalised[df$Reference == "Canis familiaris"] = df$CAI[df$Reference == "Canis familiaris"]/canis_familiaris_eCAI
df$normalised[df$Reference == "Eptesicus fuscus"] = df$CAI[df$Reference == "Eptesicus fuscus"]/eptesicus_fuscus_eCAI
df$normalised[df$Reference == "Desmodus rotundus"] = df$CAI[df$Reference == "Desmodus rotundus"]/desmodus_rotundus_eCAI
df$normalised[df$Reference == "Vulpes lagopus"] = df$CAI[df$Reference == "Vulpes lagopus"]/vulpes_lagopus_eCAI

p = ggplot(data = df, aes(x = Reference, y = normalised, fill = clade))+
  geom_boxplot(alpha = 0.9)+ 
  theme_bw()+ 
scale_fill_manual(values = mypal, labels = mylabels,
                   name = "Clade")+
scale_x_discrete(labels=c(expression(italic("C. familiaris")), expression(italic("V. lagopus")),
                          expression(italic("D. rotundus")), expression(italic("E. fuscus"))))+
  ylab("nCAI") +
  xlab("Reference host") +
  guides(fill = "none")

p

source("code/compare_cpg_cai.R")

p2 = ggpubr::ggarrange(p, p1, common.legend = T, legend = "right", labels = "AUTO")
p2

png("plots/Figure 4.png", width = 10, height = 5, units = 'in', res = 600)
p2
dev.off()
