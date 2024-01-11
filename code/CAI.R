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



df$clade = factor(df$clade, c("Dog (AF1b)","Dog (SEA2a)", 
                              "Chinese ferret badger",
                             "Mongoose", "Skunk",
                             "Big brown bat", "Hoary bat", "Free-tailed bat",
                            "Vampire bat"))

my_pal <- c("#648FFF", "#DC267F", "#FFB000")

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
                            "Asian SEA2a\n(dog)",
                            "Asian SEA2b\n(CFB)",
                            "Cosmo AM2a\n(mongoose)",
                            "RAC-SK SCSC\n(skunk)",
                            "Bat EF-E2\n(big brown bat)",
                            "Bat LC\n(hoary bat)",
                            "Bat TB1\n(Mexican free-tailed bat)",
                            "Bat DR\n(vampire bat)"))+ 
  ylab("Codon Adapatation Index") +
  xlab("Clade") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = ggtext::element_markdown(), legend.position = "bottom")
