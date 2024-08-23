library(ggplot2)
library(readxl)

df = as.data.frame(read_excel("output_data/Codon_usage_N.xlsx"))
for(i in 2:ncol(df)){
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

metadata = read.csv("sequence_data/metadata.csv")

df$clade = NA
for(i in 1:nrow(df)){
  df$clade[i] = metadata$Clade[metadata$Accession==df$`Accession no.`[i]]
}

df$clade = factor(df$clade, c("Cosmo AF1b", "Cosmo AM2a", "Arctic A", "Asian SEA2a", 
                              "Asian SEA2b", 
                              "Bat TB1",
                              "Bat DR", "Bat EF-E2","RAC-SK SCSK", "Bat LC"))


p = ggplot(data = df, aes(x = clade, y = ENC))+
  geom_boxplot()+ 
  geom_jitter(aes(color = clade), size  = 0.5, 
              width = 0.4, height = 0) + theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")+ 
  xlab("Clade") +
  scale_color_manual(values = c("#332288","#88CCEE","#CCDDAA","#44AA99","#117733",  
                                "#999933", "#DDCC77","#CC6677","#882255","#AA4499"), name = "Clade", guide = guide_legend(),
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
  scale_x_discrete(labels = c("Cosmo AF1b\n(dog)",
                              "Cosmo AM2a\n(mongoose)",
                              "Arctic A\n(arctic fox)",
                              "Asian SEA2a\n(dog)",
                              "Asian SEA2b\n(CFB)",
                              "Bat TB1\n(Mexican free\n-tailed bat)",
                              "Bat DR\n(vampire bat)",
                              "Bat EF-E2\n(big brown bat)",
                              "RAC-SK SCSK\n(skunk)",
                              "Bat LC\n(hoary bat)"))

p

png("plots/Figure 2.png", width = 7.5, height = 5, units = 'in', res = 600)
p
dev.off()

mean(df$ENC)
sd(df$ENC)

mean(df$ENC[df$clade == "Cosmo AF1b"])
mean(df$ENC[df$clade == "Asian SEA2b"])
