###need to run enc-gc3 script first to get df2!

library(readxl)

df = as.data.frame(read_excel("data/Nucleotide_composition_N.xlsx"))
df$GC3s = df$`%G3+C3`/100
df$bat = NA
df$bat[df$Host %in% c("Big brown bat", "Hoary bat", "Free-tailed bat",
                      "Vampire bat")] = "Bats"
df$bat[df$Host %in% c("Dog (SEA2a)", "Chinese ferret badger",
                      "Dog (AF1b)", "Mongoose", "Skunk")] = "Carnivores"
df$Host = factor(df$Host, c("Dog (AF1b)", "Mongoose","Dog (SEA2a)", "Chinese ferret badger",
                            "Free-tailed bat",
                            "Vampire bat", "Big brown bat","Skunk", "Hoary bat"))

library(ggpubr)
library(ggplot2)

f1 = function(x){
  2+x+29/(x^2+(1-x)^2)
}

my_pal <- c("#332288","#88CCEE","#44AA99","#117733","#999933",
            "#DDCC77","#CC6677","#882255","#AA4499")

df$GC12s = (df2$`%G1+C1` + df2$`%G2+C2`)/200
df$GC3s = df2$`%G3+C3`/100

lm(data=df[df$bat == "Bats",], GC12s ~ GC3s)
lm(data=df[df$bat == "Carnivores",], GC12s ~ GC3s)
lm(data=df, GC12s ~ GC3s)

p2 = ggplot(data = df, aes(x = GC3s, y = GC12s, col = Host, 
                           group = bat, linetype = bat)) + 
  theme_bw() + geom_point(position = position_jitter(width = .0005, height = 0.0005), 
                          size = 2, 
                          alpha = 0.8, aes(shape = bat)) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = c("Cosmo AF1b\n(dog)",
                                "Cosmo AM2a\n(mongoose)",
                                "Asian SEA2a\n(dog)",
                                "Asian SEA2b\n(CFB)",
                                "Bat TB1\n(Mexican free\n-tailed bat)",
                                "Bat DR\n(vampire bat)",
                                "Bat EF-E2\n(big brown bat)",
                                "RAC-SK SCSK\n(skunk)",
                                "Bat LC\n(hoary bat)")) +
  scale_shape_manual(values = c(16,17), name = "Host group") +
  stat_smooth(method = "lm", col = "black", se = F)+  
  stat_regline_equation(col = "black")+
  scale_linetype_manual(name="", 
                          breaks=c("Bats", "Carnivores"), 
                          labels = c("Bats", "Carnivores"),
                          values = c("dotted", "solid"))


#run ENC-GC3 script to get p1
ggarrange(p1, p2, labels = c("A", "B"), common.legend = T, legend = "bottom")


png("plots/Figure 4.png", width = 9, height = 5, units = 'in', res = 600)
ggarrange(p1, p2, labels = c("A", "B"), common.legend = T, legend = "bottom")


dev.off()
