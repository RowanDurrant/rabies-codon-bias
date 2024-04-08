library(readxl)
library(ggforce)
library(janitor)
library(devtools)
library(ggbiplot)
library(ggh4x)
library(ggplot2)

df = as.data.frame(read_excel("output_data/Codon_usage_N.xlsx"))
for(i in 2:ncol(df)){
  colnames(df)[i] = paste0(as.character(df[1,i]),"_",colnames(df)[i])
}
colnames(df)[1] = "Accession no."
df = df[2:nrow(df),]
rownames = df$`Accession no.`
df = df[,2:ncol(df)]

rownames(df) = rownames

df$F_codons = NA
df$L_codons = NA
df$I_codons = NA
df$V_codons = NA
df$S_codons = NA
df$P_codons = NA
df$T_codons = NA
df$A_codons = NA
df$Y_codons = NA
df$H_codons = NA
df$Q_codons = NA
df$N_codons = NA
df$K_codons = NA
df$D_codons = NA
df$E_codons = NA
df$C_codons = NA
df$R_codons = NA
df$G_codons = NA

for(k in 1:nrow(df)){
  df$F_codons[k] = sum(as.numeric(df[k,1:2]))
  df$L_codons[k] = sum(as.numeric(df[k,3:8]))
  df$I_codons[k] = sum(as.numeric(df[k,9:11]))
  df$V_codons[k] = sum(as.numeric(df[k,12:15]))
  df$S_codons[k] = sum(as.numeric(df[k,16:21]))
  df$P_codons[k] = sum(as.numeric(df[k,22:25]))
  df$T_codons[k] = sum(as.numeric(df[k,26:29]))
  df$A_codons[k] = sum(as.numeric(df[k,30:33]))
  df$Y_codons[k] = sum(as.numeric(df[k,34:35]))
  df$H_codons[k] = sum(as.numeric(df[k,36:37]))
  df$Q_codons[k] = sum(as.numeric(df[k,38:39]))
  df$N_codons[k] = sum(as.numeric(df[k,40:41]))
  df$K_codons[k] = sum(as.numeric(df[k,42:43]))
  df$D_codons[k] = sum(as.numeric(df[k,44:45]))
  df$E_codons[k] = sum(as.numeric(df[k,46:47]))
  df$C_codons[k] = sum(as.numeric(df[k,48:49]))
  df$R_codons[k] = sum(as.numeric(df[k,51:55]))
  df$G_codons[k] = sum(as.numeric(df[k,56:59]))
}


df <- sapply(df, as.numeric )
df = df[,60:ncol(df)]

df = remove_constant(df)

pc <- prcomp(df,
             center = TRUE,
             scale. = TRUE)
attributes(pc)

metadata = read.csv("sequence_data/metadata.csv")

g <- ggbiplot(pc,
              obs.scale = 1,
              var.scale = 1,
              groups = metadata$Clade,
              ellipse = F,
              circle = F)
g <- g + scale_color_discrete(name = '') + theme_bw() +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top')
print(g)

df2 = as.data.frame(pc$x)
df2$clade = metadata$Clade
df2$Accession = metadata$Accession

df2$clade = factor(df2$clade, c("Cosmo AF1b", "Cosmo AM2a", "Arctic A", "Asian SEA2a", 
                                "Asian SEA2b", 
                                "Bat TB1",
                                "Bat DR", "Bat EF-E2","RAC-SK SCSK", "Bat LC"))


my_pal <- c("#332288","#88CCEE","#CCDDAA","#44AA99","#117733",  
            "#999933", "#DDCC77","#CC6677","#882255","#AA4499")

g1 = ggplot(data = df2, aes(x = PC1, y = PC2))+ 
  geom_point(size = 2, aes(col = clade, shape = clade)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
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
  xlab("PC1 (33.1% explained var.)") + 
  ylab("PC2 (15.0% explained var.)")+
  theme_bw() + ylim(-10, 10) + xlim(-10, 10)+
  coord_axes_inside(labels_inside = TRUE) +
  theme(legend.position = "bottom", legend.box = "vertical")

g1
