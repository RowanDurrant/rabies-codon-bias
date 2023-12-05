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

p1 = ggplot(data = df, aes(x = GC3s, y = ENC, colour = Host))+
  geom_point(size = 2, aes(shape = bat)) + theme_bw() +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = c("Bat EF-E2 (big brown bat)",
                                "Bat LC (hoary bat)",
                                "Bat TB1 (Mexican free-tailed bat)",
                                "Bat DR (vampire bat)",
                                "Asian SEA2a (dog)",
                                "Asian SEA2b (CFB)",
                                "Cosmo AF1b (dog)",
                                "Cosmo AM2a (mongoose)",
                                "RAC-SK SCSC (skunk)")) +
  scale_shape_manual(values = c(16,17), name = "Host group") +
  stat_function(fun=f1, col = "black") +
  xlim(0,1) + ylim(0,62)
p1

