library(readxl)
library(ggpubr)
library(ggplot2)

df = as.data.frame(read_excel("output_data/Codon_usage_N.xlsx"))
for(i in 3:ncol(df)){
  colnames(df)[i] = paste0(as.character(df[1,i]),"_",colnames(df)[i])
}
colnames(df)[1] = "Accession no."
df = df[2:nrow(df),]

F_function = function(codons){
  n = sum(codons)
  pis = (codons/n)^2
  f = (n*sum(pis)-1)/(n-1)
  return(f)
}

df$ENC = NA
for(k in 1:nrow(df)){
  F_codons = as.numeric(df[k,2:3])
  L_codons = as.numeric(df[k,4:9])
  I_codons = as.numeric(df[k,10:12])
  V_codons = as.numeric(df[k,13:16])
  S_codons = as.numeric(df[k,17:22])
  P_codons = as.numeric(df[k,23:26])
  T_codons = as.numeric(df[k,27:30])
  A_codons = as.numeric(df[k,31:34])
  Y_codons = as.numeric(df[k,35:36])
  H_codons = as.numeric(df[k,37:38])
  Q_codons = as.numeric(df[k,39:40])
  N_codons = as.numeric(df[k,41:42])
  K_codons = as.numeric(df[k,43:44])
  D_codons = as.numeric(df[k,45:46])
  E_codons = as.numeric(df[k,47:48])
  C_codons = as.numeric(df[k,49:50])
  R_codons = as.numeric(df[k,52:56])
  G_codons = as.numeric(df[k,57:60])
  
  #F2 = C, D, E, F, H, K, N, Q, Y
  #F3 = I
  #F4 = A, G, P, T, V
  #F6 = R, L, S
  
  F2s = c(F_function(C_codons), F_function(D_codons), F_function(E_codons), F_function(F_codons),
          F_function(H_codons), F_function(K_codons), F_function(N_codons), F_function(Q_codons),
          F_function(Y_codons))
  F3s = c(F_function(I_codons))
  F4s = c(F_function(A_codons), F_function(G_codons), F_function(P_codons), F_function(T_codons),
          F_function(V_codons))
  F6s = c(F_function(L_codons), F_function(S_codons), F_function(R_codons))
  
  ENC = 2 + (9/mean(F2s)) + (1/mean(F3s) + (5/mean(F4s)) + (3/mean(F6s)))
  df$ENC[k] = ENC
}

df2 = as.data.frame(read_excel("output_data/Nucleotide_composition_N.xlsx"))
df$GC3s = df2$`%G3+C3`/100

metadata = read.csv("sequence_data/metadata.csv")

df$clade = NA
for(i in 1:nrow(df)){
  df$clade[i] = metadata$Clade[metadata$Accession==df$`Accession no.`[i]]
}

df$clade = factor(df$clade, c("Cosmo AF1b", "Cosmo AM2a", "Arctic A", "Asian SEA2a", 
                              "Asian SEA2b", 
                              "Bat TB1",
                              "Bat DR", "Bat EF-E2","RAC-SK SCSK", "Bat LC"))

f1 = function(x){
  2+x+29/(x^2+(1-x)^2)
}

my_pal <- c("#332288","#88CCEE","#CCDDAA","#44AA99","#117733",  
            "#999933", "#DDCC77","#CC6677","#882255","#AA4499")

p1 = ggplot(data = df, aes(x = GC3s, y = ENC))+
  geom_point(size = 2, aes(colour = clade, shape = clade)) + theme_bw() +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = c("Cosmo AF1b\n(dog)",
                                "Cosmo AM2a\n(mongoose)",
                                "Arctic A\n(arctic fox)",
                                "Asian SEA2a\n(dog)",
                                "Asian SEA2b\n(CFB)",
                                "Bat TB1\n(Mexican free\n-tailed bat)",
                                "Bat DR\n(vampire bat)",
                                "Bat EF-E2\n(big brown bat)",
                                "RAC-SK SCSK\n(skunk)",
                                "Bat LC\n(hoary bat)")) +
  scale_shape_manual(name = "Clade",
                     labels = c("Cosmo AF1b\n(dog)",
                                "Cosmo AM2a\n(mongoose)",
                                "Arctic A\n(arctic fox)",
                                "Asian SEA2a\n(dog)",
                                "Asian SEA2b\n(CFB)",
                                "Bat TB1\n(Mexican free\n-tailed bat)",
                                "Bat DR\n(vampire bat)",
                                "Bat EF-E2\n(big brown bat)",
                                "RAC-SK SCSK\n(skunk)",
                                "Bat LC\n(hoary bat)"),
                     values = c(17,17,17,17,17,16,16,16,17,16))+
  stat_function(fun=f1, col = "black") +
  xlim(0.35,0.6) + ylim(48,62)
p1

df = as.data.frame(read_excel("output_data/Nucleotide_composition_N.xlsx"))
df$GC3s = df$`%G3+C3`/100

df$clade = NA
for(i in 1:nrow(df)){
  df$clade[i] = metadata$Clade[metadata$Accession==df$`Accession no.`[i]]
}

df$bat = NA
df$bat[df$clade %in% c("Bat TB1",
                       "Bat DR", "Bat EF-E2","Bat LC")] = "Bats"
df$bat[df$clade %in% c("Cosmo AF1b", "RAC-SK SCSK", "Cosmo AM2a", "Arctic A", "Asian SEA2a", 
                       "Asian SEA2b")] = "Carnivores"
df$clade = factor(df$clade, c("Cosmo AF1b", "Cosmo AM2a", "Arctic A", "Asian SEA2a", 
                              "Asian SEA2b", 
                              "Bat TB1",
                              "Bat DR", "Bat EF-E2","RAC-SK SCSK", "Bat LC"))

df$GC12s = (df2$`%G1+C1` + df2$`%G2+C2`)/200
df$GC3s = df2$`%G3+C3`/100

lm(data=df[df$bat == "Bats",], GC12s ~ GC3s)
lm(data=df[df$bat == "Carnivores",], GC12s ~ GC3s)
lm(data=df, GC12s ~ GC3s)

p2 = ggplot(data = df, aes(x = GC3s, y = GC12s,  
                           linetype = bat)) + 
  theme_bw() + 
  geom_point(position = position_jitter(width = .0005, height = 0.0005), 
             size = 2, 
             alpha = 0.8, aes(col = clade, shape = clade)) +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = c("Cosmo AF1b\n(dog)",
                                "Cosmo AM2a\n(mongoose)",
                                "Arctic A\n(arctic fox)",
                                "Asian SEA2a\n(dog)",
                                "Asian SEA2b\n(CFB)",
                                "Bat TB1\n(Mexican free\n-tailed bat)",
                                "Bat DR\n(vampire bat)",
                                "Bat EF-E2\n(big brown bat)",
                                "RAC-SK SCSK\n(skunk)",
                                "Bat LC\n(hoary bat)")) +
  scale_shape_manual(name = "Clade",
                     labels = c("Cosmo AF1b\n(dog)",
                                "Cosmo AM2a\n(mongoose)",
                                "Arctic A\n(arctic fox)",
                                "Asian SEA2a\n(dog)",
                                "Asian SEA2b\n(CFB)",
                                "Bat TB1\n(Mexican free\n-tailed bat)",
                                "Bat DR\n(vampire bat)",
                                "Bat EF-E2\n(big brown bat)",
                                "RAC-SK SCSK\n(skunk)",
                                "Bat LC\n(hoary bat)"),
                     values = c(17,17,17,17,17,16,16,16,17,16))+
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

