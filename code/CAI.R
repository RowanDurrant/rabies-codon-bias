library("readxl")
library("ggplot2")
library("ggtext")

df = read.csv("HIVE-CUTS_CAI.csv")

reference = as.data.frame(read_excel("CAI_N.xlsx"))

reference$Clade[reference$Clade == "CFB"] = "Chinese ferret badger"
reference$Clade[reference$Clade == "Dog SEA2a"] = "Dog (SEA2a)"
reference$Clade[reference$Clade == "VB"] = "Vampire bat"
reference$Clade[reference$Clade == "Dog AF1b"] = "Dog (AF1b)"
reference$Clade[reference$Clade == "EF"] = "Big brown bat"
reference$Clade[reference$Clade == "SCSC"] = "Skunk"
reference$Clade[reference$Clade == "TB"] = "Free-tailed bat"


df$clade = reference$Clade



df$clade = factor(df$clade, c("Dog (AF1b)", "Mongoose","Dog (SEA2a)", "Chinese ferret badger",
                              "Free-tailed bat",
                              "Vampire bat", "Big brown bat","Skunk", "Hoary bat"))

my_pal <- c("#332288","#88CCEE","#44AA99","#117733","#999933",
            "#DDCC77","#CC6677","#882255","#AA4499")

library("tidyr")
df2 <- pivot_longer(df, cols=c(Canis_familiaris, Eptesicus_fuscus, Desmodus_rotundus)) 

ggplot(data = df2, aes(x = clade, y = value, color = name))+
  geom_boxplot()+ 
  # geom_jitter(alpha  = 0.4, 
  #       width = 0.4, height = 0) +
  theme_bw()+ 
  scale_color_manual(values = my_pal, labels = c("Canis familiaris",
                     "Desmodus rotundus", "Eptesicus fuscus"),
                     name = "Reference host")+
  scale_x_discrete(labels=c("Cosmo AF1b\n(dog)",
                            "Cosmo AM2a\n(mongoose)",
                            "Asian SEA2a\n(dog)",
                            "Asian SEA2b\n(CFB)",
                            "Bat TB1\n(Mexican free\n-tailed bat)",
                            "Bat DR\n(vampire bat)",
                            "Bat EF-E2\n(big brown bat)",
                            "RAC-SK SCSC\n(skunk)",
                            "Bat LC\n(hoary bat)"))+ 
  ylab("Codon Adapatation Index") +
  xlab("Clade") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = ggtext::element_markdown(), legend.position = "bottom")


canis_familiaris_eCAI = 0.785
desmodus_rotundus_eCAI = 0.753
eptesicus_fuscus_eCAI = 0.755

df2$normalised = NA
df2$normalised[df2$name == "Canis_familiaris"] = df2$value[df2$name == "Canis_familiaris"]/canis_familiaris_eCAI
df2$normalised[df2$name == "Eptesicus_fuscus"] = df2$value[df2$name == "Eptesicus_fuscus"]/eptesicus_fuscus_eCAI
df2$normalised[df2$name == "Desmodus_rotundus"] = df2$value[df2$name == "Desmodus_rotundus"]/desmodus_rotundus_eCAI

ggplot(data = df2, aes(x = name, y = normalised, fill = clade))+
  geom_boxplot(alpha = 0.9)+ 
  theme_bw()+ 
scale_fill_manual(values = my_pal, labels = c("Cosmo AF1b\n(dog)",
                                               "Cosmo AM2a\n(mongoose)",
                                               "Asian SEA2a\n(dog)",
                                               "Asian SEA2b\n(CFB)",
                                               "Bat TB1\n(Mexican free\n-tailed bat)",
                                               "Bat DR\n(vampire bat)",
                                               "Bat EF-E2\n(big brown bat)",
                                               "RAC-SK SCSC\n(skunk)",
                                               "Bat LC\n(hoary bat)"),
                   name = "Clade")+
scale_x_discrete(labels=c("Canis familiaris",
                                               "Desmodus rotundus", "Eptesicus fuscus"))+
  ylab("Normalised codon adapatation index") +
  xlab("Reference host") +
  theme(legend.position = "bottom")

