library(readxl)
df = as.data.frame(read_excel("data/Codon_usage_N.xlsx"))
for(i in 3:ncol(df)){
  colnames(df)[i] = paste0(as.character(df[1,i]),"_",colnames(df)[i])
}
colnames(df)[2] = "Accession no."
df = df[2:nrow(df),]
rownames = df$`Accession no.`
df = df[,3:ncol(df)]

df <- sapply(df, as.numeric )
rownames(df) = rownames

library(janitor)
df = remove_constant(df)

pc <- prcomp(df,
             center = TRUE,
             scale. = TRUE)
attributes(pc)

df = as.data.frame(read_excel("data/Codon_usage_N.xlsx"))
for(i in 3:ncol(df)){
  colnames(df)[i] = paste0(as.character(df[1,i]),"_",colnames(df)[i])
}
colnames(df)[2] = "Accession no."
df = df[2:nrow(df),]

library(devtools)
library(ggbiplot)
g <- ggbiplot(pc,
              obs.scale = 1,
              var.scale = 1,
              groups = df$Host,
              ellipse = F,
              circle = F)
g <- g + scale_color_discrete(name = '') + theme_bw() +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top')
print(g)

df2 = as.data.frame(pc$x)
df2$Host = df$Host
df2$Accession = df$`Accession no.`

df2$Host = factor(df$Host, c("Dog (AF1b)", "Mongoose","Dog (SEA2a)", "Chinese ferret badger",
                             "Free-tailed bat",
                             "Vampire bat", "Big brown bat","Skunk", "Hoary bat"))
df2$bat = NA
df2$bat[df2$Host %in% c("Big brown bat", "Hoary bat", "Free-tailed bat",
                      "Vampire bat")] = "Bats"
df2$bat[df2$Host %in% c("Dog (SEA2a)", "Chinese ferret badger",
                      "Dog (AF1b)", "Mongoose", "Skunk")] = "Carnivores"
library(ggh4x)

my_pal <- c("#332288","#88CCEE","#44AA99","#117733","#999933",
            "#DDCC77","#CC6677","#882255","#AA4499")

library(ggforce)

g1 = ggplot(data = df2, aes(x = PC1, y = PC2, colour = Host))+ 
  geom_point(size = 2, aes(shape = bat)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  scale_color_manual(values = my_pal, name = "Clade",
                     labels = c("Cosmo AF1b\n(dog)",
                                "Cosmo AM2a\n(mongoose)",
                                "Asian SEA2a\n(dog)",
                                "Asian SEA2b\n(CFB)",
                                "Bat TB1\n(Mexican free\n-tailed bat)",
                                "Bat DR\n(vampire bat)",
                                "Bat EF-E2\n(big brown bat)",
                                "RAC-SK SCSK\n(skunk)",
                                "Bat LC\n(hoary bat)")) +
  scale_fill_manual(values = my_pal, name = "Host species") +
  scale_shape_manual(values = c(16,17), name = "Host group") +
  xlab("PC1 (25.8% explained var.)") + 
  ylab("PC2 (22.6% explained var.)")+
  theme_bw() + ylim(-10, 10) + xlim(-10, 10)+
  coord_axes_inside(labels_inside = TRUE) +
  theme(legend.position = "bottom", legend.box = "vertical")

g1

write.csv(df2, "PCA_output.csv")

# Helper function 
#::::::::::::::::::::::::::::::::::::::::
var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
}
# Compute Coordinates
#::::::::::::::::::::::::::::::::::::::::
loadings <- pc$rotation
sdev <- pc$sdev
var.coord <- t(apply(loadings, 1, var_coord_func, sdev))
var.cos2 <- var.coord^2
comp.cos2 <- apply(var.cos2, 2, sum)
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))

pc2 = var.contrib[,"PC2"]
barplot(pc2)
max(pc2)
names(pc2[pc2>5])

pc1 = var.contrib[,"PC1"]
barplot(pc1)
max(pc1)
names(pc1[pc1>4.5])
