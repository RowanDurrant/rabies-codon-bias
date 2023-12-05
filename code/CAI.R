library("readxl")
library("ggplot2")
library("ggtext")

df = as.data.frame(read_excel("CAI_N.xlsx"))

df$Clade[df$Clade == "CFB"] = "Chinese ferret badger"
df$Clade[df$Clade == "Dog SEA2a"] = "Dog (SEA2a)"
df$Clade[df$Clade == "VB"] = "Vampire bat"
df$Clade[df$Clade == "Dog AF1b"] = "Dog (AF1b)"
df$Clade[df$Clade == "EF"] = "Big brown bat"
df$Clade[df$Clade == "SCSC"] = "Skunk"
df$Clade[df$Clade == "TB"] = "Free-tailed bat"

df$Clade = factor(df$Clade, c("Dog (AF1b)","Dog (SEA2a)", 
                              "Chinese ferret badger",
                             "Mongoose", "Skunk",
                             "Big brown bat", "Hoary bat", "Free-tailed bat",
                            "Vampire bat"))

my_pal <- c("#b66dff","#490092","#006ddb","#6db6ff","#b6dbff",
            "#004949","#009292","#ff6db6","#ffb6db"
            
            )

ggplot(data = df, aes(x = Clade, y = CAI))+
  geom_boxplot()+ 
  geom_jitter(aes(color = Clade), alpha  = 0.4, 
        width = 0.4, height = 0) +
  theme_bw()+ scale_color_manual(values = my_pal, labels = c("Cosmo AF1b (dog)",
                                                             "Asian SEA2a (dog)",
                                                             "Asian SEA2b (CFB)",
                                                             "Cosmo AM2a (mongoose)",
                                                             "RAC-SK SCSC (skunk)",
                                                             "Bat EF-E2 (big brown bat)",
                                                             "Bat LC (hoary bat)",
                                                             "Bat TB1 (Mexican free-tailed bat)",
                                                             "Bat DR (vampire bat)")) +
  scale_x_discrete(labels=c("Cosmo AF1b\n(dog)",
                            "Asian SEA2a\n(dog)",
                            "Asian SEA2b\n(CFB)",
                            "Cosmo AM2a\n(mongoose)",
                            "RAC-SK SCSC\n(skunk)",
                            "Bat EF-E2\n(big brown bat)",
                            "Bat LC\n(hoary bat)",
                            "Bat TB1\n(Mexican free-tailed bat)",
                            "Bat DR\n(vampire bat)"))+ 
  ylab("CAI (*Canis familiaris* reference)") +
  xlab("Host species") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none",
        axis.title.y = ggtext::element_markdown())
