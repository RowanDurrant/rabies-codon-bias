#PR2 plot
#AU-bias and GC-bias at 3rd codon position
#this should only be for 4-codon amino acids
#F4 = A, G, P, T, V

library("readxl")
df = as.data.frame(read_excel("Codon_usage.xlsx"))
for(i in 3:ncol(df)){
  colnames(df)[i] = paste0(as.character(df[1,i]),"_",colnames(df)[i])
}
colnames(df)[2] = "Accession no."
df = df[2:nrow(df),]


df$Host = factor(df$Host, c("Big brown bat", "Hoary bat", "Free-tailed bat",
                             "Vampire bat", "Dog (SEA2a)", "Chinese ferret badger",
                             "Dog (AF1b)", "Mongoose", "Skunk"))
df$bat = NA
df$bat[df$Host %in% c("Big brown bat", "Hoary bat", "Free-tailed bat",
                        "Vampire bat")] = "Bats"
df$bat[df$Host %in% c("Dog (SEA2a)", "Chinese ferret badger",
                        "Dog (AF1b)", "Mongoose", "Skunk")] = "Carnivores"

df$A3 = NA
df$C3 = NA
df$T3 = NA
df$G3 = NA

for(j in 1:nrow(df)){
  df$A3[j] = sum(as.numeric(df$A_GCA[j]), as.numeric(df$G_GGA[j]), as.numeric(df$P_CCA[j]), 
                 as.numeric(df$T_ACA[j]), as.numeric(df$V_GTA[j]))
  df$C3[j] = sum(as.numeric(df$A_GCC[j]), as.numeric(df$G_GGC[j]), as.numeric(df$P_CCC[j]), 
                 as.numeric(df$T_ACC[j]), as.numeric(df$V_GTC[j]))
  df$T3[j] = sum(as.numeric(df$A_GCT[j]), as.numeric(df$G_GGT[j]), as.numeric(df$P_CCT[j]), 
                 as.numeric(df$T_ACT[j]), as.numeric(df$V_GTT[j]))
  df$G3[j] = sum(as.numeric(df$A_GCG[j]), as.numeric(df$G_GGG[j]), as.numeric(df$P_CCG[j]), 
                 as.numeric(df$T_ACG[j]), as.numeric(df$V_GTG[j]))
}

df$A3prop = df$A3/(df$A3 + df$T3)
df$G3prop = df$G3/(df$G3 + df$C3)

my_pal <- c("#004949","#009292","#ff6db6","#ffb6db",
            "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
            "#920000","#924900","#db6d00","#24ff24","#ffff6d",
            "#555555", "#999999")

p3 = ggplot(data = df, aes(x = G3prop, y = A3prop, col = Host))+
  geom_point(size = 2, aes(shape = bat)) + theme_bw() + xlim(0,1) + ylim(0,1)+
  scale_color_manual(values = my_pal, name = "Host species") +
  scale_shape_manual(values = c(16,17), name = "Host group") +
  geom_hline(yintercept = 0.5) + geom_vline(xintercept = 0.5)+
  xlab("G3/(G3+C3)") + ylab("A3/(A3+U3)")#+
  # theme(legend.position = "bottom", legend.box = "vertical")

p4 = ggplot(data = df, aes(x = G3prop, y = A3prop, col = Host))+
  geom_point(size = 2, aes(shape = bat)) + theme_bw() +
  scale_color_manual(values = my_pal, name = "Host species") +
  scale_shape_manual(values = c(16,17), name = "Host group") +
  xlab("G3/(G3+C3)") + ylab("A3/(A3+U3)")

ggarrange(p3, p4, labels = c("A", "B"), common.legend = T, legend = "bottom")
