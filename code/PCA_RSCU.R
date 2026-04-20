library(readxl)
library(ggforce)
library(janitor)
library(devtools)
library(ggbiplot)
library(ggh4x)
library(ggplot2)
library(stringr)
library(ggpubr)

df = as.data.frame(read_excel("output_data/RSCU_N.xlsx"))
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

pc <- prcomp(df,
             center = TRUE,
             scale. = TRUE)
summary(pc)

library(factoextra)
fviz_eig(pc, 
         addlabels = TRUE)

var_explained = pc$sdev^2 / sum(pc$sdev^2)

qplot(c(1:10), var_explained[1:10]) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("") +
  ylim(0, 0.3)+
  theme_bw()+
  scale_x_continuous(breaks = 1:10)+
  geom_text(aes(label = signif(var_explained[1:10],3)), hjust=0.25, vjust=-0.5)

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
df2$clade = NA
for(i in 1:nrow(df2)){
  df2$clade[i] = metadata$Clade[metadata$Accession==rownames(df2[i,])]
}

df2$clade = factor(df2$clade, c("Cosmopolitan AF1b", "Cosmopolitan AM2a", "Arctic A", "Asian SEA2a", 
                                "Asian SEA2b", 
                                "Bats TB1",
                                "Bats DR", "Bats EF-E2","RAC-SK SCSK", "Bats LC"))


my_pal <- c("#332288","#88CCEE","#CCDDAA","#44AA99","#117733",  
            "#999933", "#DDCC77","#CC6677","#882255","#AA4499")
my_labels = c("Cosmo AF1b\n(dog)",
              "Cosmo AM2a\n(mongoose)",
              "Arctic A\n(arctic fox)",
              "Asian SEA2a\n(dog)",
              "Asian SEA2b\n(CFB)",
              "Bat DR\n(vampire bat)",
              "Bat TB1\n(Mexican free\n-tailed bat)",
              "Bat LC\n(hoary bat)",
              "Bat EF-E2\n(big brown bat)",
              "RAC-SK SCSK\n(skunk)")

g1 = ggplot(data = df2, aes(x = PC1, y = PC2))+ 
  geom_point(size = 2, alpha = 0.8, aes(col = clade, shape = clade)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  #ylim(-6,6) + xlim(-8,8)+
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  xlab("PC1 (22.4% explained var.)") + 
  ylab("PC2 (20.1% explained var.)")+
  theme_bw() + 
  coord_axes_inside(labels_inside = TRUE) +
  theme(legend.position = "none")

g2 = ggplot(data = df2, aes(x = PC1, y = PC3))+ 
  geom_point(size = 2, alpha = 0.8, aes(col = clade, shape = clade)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  #ylim(-6,6) + xlim(-8,8)+
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  xlab("PC1 (22.4% explained var.)") + 
  ylab("PC3 (17.1% explained var.)")+
  theme_bw() + 
  coord_axes_inside(labels_inside = TRUE) +
  theme(legend.position = "none")

g3 = ggplot(data = df2, aes(x = PC2, y = PC3,col = clade, shape = clade))+ 
  geom_point(size = 2, alpha = 0.8) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  #ylim(-6,6) + xlim(-8,8)+
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = my_labels) +
  scale_shape_manual(name = "Clade",
                     labels = my_labels,
                     values = c(17,17,17,17,17,16,16,16,16,17))+
  xlab("PC2 (20.1% explained var.)") + 
  ylab("PC3 (17.1% explained var.)")+
  theme_bw() + 
  coord_axes_inside(labels_inside = TRUE) +
  theme(legend.position = "none")

ggpubr::ggarrange(g1,g2,g3, nrow = 1, common.legend = T, legend = "bottom",
                  labels = "AUTO", align = "v")

write.csv(df2, "output_data/PCA_output_RSCU.csv")
