library("readxl")
df = as.data.frame(read_excel("Codon_usage.xlsx"))
for(i in 3:ncol(df)){
  colnames(df)[i] = paste0(as.character(df[1,i]),"_",colnames(df)[i])
}
colnames(df)[2] = "Accession no."
df = df[2:nrow(df),]

F_function = function(codons){
  n = sum(codons)
  pis = (codons/n)^2
  f = (n*sum(pis)-1)/(n-1)
  return(f)
}

df$ENC = NA
for(k in 1:nrow(df)){
  F_codons = as.numeric(df[k,3:4])
  L_codons = as.numeric(df[k,5:10])
  I_codons = as.numeric(df[k,11:13])
  V_codons = as.numeric(df[k,14:17])
  S_codons = as.numeric(df[k,18:23])
  P_codons = as.numeric(df[k,24:27])
  T_codons = as.numeric(df[k,28:31])
  A_codons = as.numeric(df[k,32:35])
  Y_codons = as.numeric(df[k,36:37])
  H_codons = as.numeric(df[k,38:39])
  Q_codons = as.numeric(df[k,40:41])
  N_codons = as.numeric(df[k,42:43])
  K_codons = as.numeric(df[k,44:45])
  D_codons = as.numeric(df[k,46:47])
  E_codons = as.numeric(df[k,48:49])
  C_codons = as.numeric(df[k,50:51])
  R_codons = as.numeric(df[k,52:57])
  G_codons = as.numeric(df[k,58:61])
  
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

df$Host = factor(df$Host, c("Big brown bat", "Hoary bat", "Free-tailed bat",
                              "Vampire bat", "Dog (SEA2a)", "Chinese ferret badger",
                              "Dog (AF1b)", "Mongoose", "Skunk"))


library(ggplot2)
ggplot(data = df, aes(x = Host, y = ENC))+
  geom_boxplot()+ theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


df2 = as.data.frame(read_excel("Nucleotide_composition.xlsx"))
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
  scale_color_manual(values = my_pal, name = "Host species") +
  scale_shape_manual(values = c(16,17), name = "Host group") +
  stat_function(fun=f1, col = "black") +
  xlim(0,1) + ylim(0,62)
# p2 = ggplot(data = df, aes(x = GC3s, y = ENC, colour = Host))+
#   geom_point(size = 2, shape = bat) + theme_bw() +
#   scale_color_manual(values = my_pal, name = "Host species") +
#   scale_shape_manual(values = c(16,17), name = "Host group") +
#   stat_function(fun=f1, col = "black")
# library(ggpubr)
# ggarrange(p1, p2, labels = c("A", "B"), common.legend = T, legend = "bottom")
p1

df$GC12s = (df2$`%G1+C1` + df2$`%G2+C2`)/200
df$GC3s = df2$`%G3+C3`/100

p5 = ggplot(data = df, aes(x = GC3s, y = GC12s, col = Host)) + 
  theme_bw() + geom_point(size = 2, aes(shape = bat)) +
  ylim(0.42, 0.445)+ xlim(0.425, 0.515)+
  scale_color_manual(values = my_pal, name = "Host species") +
  scale_shape_manual(values = c(16,17), name = "Host group") +
  stat_smooth(method = "lm", col = "black", se = F)+  
  stat_regline_equation(#label.y = 0.4425,label.x= 0.44,
    col = "black"
    )

p6 = ggplot(data = df, aes(x = GC3s, y = GC12s, col = Host)) + 
  theme_bw() + geom_point(aes(shape = bat)) +
  ylim(0.42, 0.443)+ xlim(0.428, 0.505)+
  scale_color_manual(values = my_pal, name = "Host species") +
  scale_shape_manual(values = c(16,17), name = "Host group") +
  stat_smooth(method = "lm", col = "black", se = F)+ 
  stat_regline_equation(label.y = 0.4415,label.x= 0.43,
                        col = "black", size = 3 )+
  facet_wrap(~Host)

ggarrange(p5,p6, labels = c("A","B"), common.legend = T, legend = "bottom")
