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

df$clade = factor(df$clade, c("Cosmopolitan AF1b", "Cosmopolitan AM2a", "Arctic A", "Asian SEA2a", 
                              "Asian SEA2b", 
                              "Bats DR","Bats TB1", "Bats LC",
                              "Bats EF-E2","RAC-SK SCSK"))

f1 = function(x){
  2+x+(29/(x^2+(1-x)^2))
}

my_pal = c("#332288","#88CCEE","#CCDDAA","#44AA99","#117733",  
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

p = ggplot(data = df, aes(x = GC3s, y = ENC))+
  geom_point(size = 2, aes(colour = clade, shape = clade)) + theme_bw() +
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = mylabels) +
  scale_shape_manual(name = "Clade",
                     labels = mylabels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  stat_function(fun=f1, col = "black") +
  xlim(0.3,0.7) + ylim(48,61)+
  xlab("GC3 content")
p

df_pca = read.csv("output_data/PCA_output_RSCU.csv")
df_pca$clade = factor(df$clade, c("Cosmopolitan AF1b", "Cosmopolitan AM2a", "Arctic A", "Asian SEA2a", 
                              "Asian SEA2b", 
                              "Bats DR","Bats TB1", "Bats LC",
                              "Bats EF-E2","RAC-SK SCSK"))
g1 = ggplot(data = df_pca, aes(x = PC1, y = PC2))+ 
  geom_point(size = 2, alpha = 0.8, aes(col = clade, shape = clade)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  ylim(-8,8) + xlim(-10,10)+
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  xlab("PC1 (22.4% explained var.)") + 
  ylab("PC2 (20.1% explained var.)")+
  theme_bw() + 
  coord_axes_inside(labels_inside = TRUE)

ggpubr::ggarrange(p, g1, common.legend = T, legend = "bottom",
                  align = "h")

png("plots/Supplementary Figure 2.png", width = 9, height = 6, units = 'in', res = 600)
ggpubr::ggarrange(p, g1, common.legend = T, legend = "bottom",
                  align = "h")
dev.off()
