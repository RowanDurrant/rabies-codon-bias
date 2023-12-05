df2 = as.data.frame(read_excel("Nucleotide_composition_N.xlsx"))
df$GC3s = df2$`%G3+C3`/100
df$bat = NA
df$bat[df$Host %in% c("Big brown bat", "Hoary bat", "Free-tailed bat",
                      "Vampire bat")] = "Bats"
df$bat[df$Host %in% c("Dog (SEA2a)", "Chinese ferret badger",
                      "Dog (AF1b)", "Mongoose", "Skunk")] = "Carnivores"

f1 = function(x){
  2+x+29/(x^2+(1-x)^2)
}

my_pal <- c("#004949","#009292","#ff6db6","#ffb6db",
            "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
            "#920000","#924900","#db6d00","#24ff24","#ffff6d",
            "#555555", "#999999")

df$GC12s = (df2$`%G1+C1` + df2$`%G2+C2`)/200
df$GC3s = df2$`%G3+C3`/100

p5 = ggplot(data = df, aes(x = GC3s, y = GC12s, col = Host)) + 
  theme_bw() + geom_point(size = 2, aes(shape = bat)) +
  ylim(0.42, 0.445)+ xlim(0.425, 0.515)+
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = c("Bat EF-E2\n(big brown bat)",
                                "Bat LC\n(hoary bat)",
                                "Bat TB1\n(Mexican free-tailed bat)",
                                "Bat DR\n(vampire bat)",
                                "Asian SEA2a\n(dog)",
                                "Asian SEA2b\n(CFB)",
                                "Cosmo AF1b\n(dog)",
                                "Cosmo AM2a\n(mongoose)",
                                "RAC-SK SCSC\n(skunk)",
                                "Gannoruwa bat lyssavirus"
                     )) +
  scale_shape_manual(values = c(16,17), name = "Host group") +
  stat_smooth(method = "lm", col = "black", se = F)+  
  stat_regline_equation(#label.y = 0.4425,label.x= 0.44,
    col = "black"
  )

clade_names <- c(
  `Big brown bat` = "Bat EF-E2\n(big brown bat)",
  `Hoary bat` = "Bat LC\n(hoary bat)",
  `Free-tailed bat` = "Bat TB1\n(MFTB)",
  `Vampire bat` = "Bat DR\n(vampire bat)",
  `Dog (SEA2a)` = "Asian SEA2a\n(dog)",
  `Chinese ferret badger` = "Asian SEA2b\n(CFB)",
  `Dog (AF1b)` = "Cosmo AF1b\n(dog)",
  `Mongoose` = "Cosmo AM2a\n(mongoose)",
  `Skunk` = "RAC-SK SCSC\n(skunk)"
)

p6 = ggplot(data = df, aes(x = GC3s, y = GC12s, col = Host)) + 
  theme_bw() + geom_point(aes(shape = bat)) +
  ylim(0.42, 0.443)+ xlim(0.428, 0.505)+
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = c("Bat EF-E2\n(big brown bat)",
                                "Bat LC\n(hoary bat)",
                                "Bat TB1\n(Mexican free-tailed bat)",
                                "Bat DR\n(vampire bat)",
                                "Asian SEA2a\n(dog)",
                                "Asian SEA2b\n(CFB)",
                                "Cosmo AF1b\n(dog)",
                                "Cosmo AM2a\n(mongoose)",
                                "RAC-SK SCSC\n(skunk)",
                                "Gannoruwa bat lyssavirus"
                     )) +
  scale_shape_manual(values = c(16,17), name = "Host group") +
  stat_smooth(method = "lm", col = "black", se = F)+ 
  stat_regline_equation(label.y = 0.4415,label.x= 0.43,
                        col = "black", size = 3 )+
  facet_wrap(~Host, labeller = as_labeller(clade_names))

ggarrange(p5,p6, labels = c("A","B"), common.legend = T, legend = "bottom")
