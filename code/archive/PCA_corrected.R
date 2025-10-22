library(readxl)
library(ggforce)
library(janitor)
library(devtools)
library(ggbiplot)
library(ggh4x)
library(ggplot2)
library(stringr)
library(ggpubr)

source("code/codon_means.R")

df = as.data.frame(read_excel("output_data/Codon_usage_N.xlsx"))
for(i in 2:ncol(df)){
  colnames(df)[i] = paste0(as.character(df[1,i]),"_",colnames(df)[i])
}
colnames(df)[1] = "Accession no."
df = df[2:nrow(df),]
rownames = df$`Accession no.`
df = df[,2:ncol(df)]

df <- sapply(df, as.numeric )
df = remove_constant(df)

rownames(df) = rownames

metadata = read.csv("sequence_data/metadata.csv")

df = df[rownames(df) %in% metadata$Accession,]

#remove amino acids with only one codon (M + W)
df = df[,1:(ncol(df)-2)]

for(i in 1:nrow(df)){
  df[i,] = codon_mean(df[i,])
  
}

pc <- prcomp(df,
             center = TRUE,
             scale. = TRUE)
attributes(pc)

var_explained = pc$sdev^2 / sum(pc$sdev^2)

qplot(c(1:10), var_explained[1:10]) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 0.3)

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
  xlab("PC1 (24.0% explained var.)") + 
  ylab("PC2 (21.7% explained var.)")+
  theme_bw() + ylim(-10, 10) + xlim(-10, 10)+
  coord_axes_inside(labels_inside = TRUE) +
  theme(legend.position = "bottom", legend.box = "vertical")

g1

write.csv(df2, "output_data/PCA_output_corrected.csv")
