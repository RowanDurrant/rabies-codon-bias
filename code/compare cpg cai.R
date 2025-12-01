library(tidyr)
library(ggplot2)
library(boot)

df = read.csv("output_data/HIVE-CUTS_CAI.csv")
cpg = read.csv("output_data/N_CpG.csv")

df = as.data.frame(cbind(df,cpg))

df = pivot_longer(df, cols= 2:5, names_to = "Reference_species", values_to = "CAI")

canis_familiaris_eCAI = 0.785
desmodus_rotundus_eCAI = 0.753
eptesicus_fuscus_eCAI = 0.755
vulpes_lagopus_eCAI = 0.778

df$nCAI = NA
df$nCAI[df$Reference_species == "Canis_familiaris"] = df$CAI[df$Reference_species == "Canis_familiaris"]/canis_familiaris_eCAI
df$nCAI[df$Reference_species == "Eptesicus_fuscus"] = df$CAI[df$Reference_species == "Eptesicus_fuscus"]/eptesicus_fuscus_eCAI
df$nCAI[df$Reference_species == "Desmodus_rotundus"] = df$CAI[df$Reference_species == "Desmodus_rotundus"]/desmodus_rotundus_eCAI
df$nCAI[df$Reference_species == "Vulpes_lagopus"] = df$CAI[df$Reference_species == "Vulpes_lagopus"]/vulpes_lagopus_eCAI

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

p1 = ggplot(data = df, aes(x = cpg, y = nCAI))+
  geom_point(aes(shape = Reference_species,
                 colour = clade))+
  geom_smooth(se = F, method = "lm", col = "black")+
  theme_bw()+ 
  scale_colour_manual(values = mypal, labels = mylabels,
                    name = "Clade")+
  scale_shape_manual(name = "Reference host",
                     values = c(15,16,17,3),
                     labels = c(expression(italic("Canis familiaris")),
                                expression(italic("Desmodus rotundus")),
                                expression(italic("Eptesicus fuscus")),
                                expression(italic("Vulpes lagopus"))))+
  xlab("Obs/Exp CpG content")

summary(lm(data = df, nCAI ~ cpg))
bootstrap = boot(df,function(data,indices)
  summary(lm(nCAI ~ cpg,data[indices,]))$adj.r.squared,R=10000)
quantile(bootstrap$t,c(0.025,0.975))
