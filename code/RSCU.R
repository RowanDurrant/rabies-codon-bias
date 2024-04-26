library(readxl)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(colorspace)

df = as.data.frame(read_excel("output_data/RSCU_N.xlsx"))
for(i in 2:ncol(df)){
  colnames(df)[i] = paste0(as.character(df[1,i]),"_",colnames(df)[i])
  df[,i] = as.numeric(df[,i])
}
colnames(df)[1] = "Accession no."
df = df[2:nrow(df),]

metadata = read.csv("sequence_data/metadata.csv")
df$clade = NA
for(i in 1:nrow(df)){
  df$clade[i] = metadata$Clade[metadata$Accession==df$`Accession no.`[i]]
}

agg = aggregate(df[, 2:60], list(df$clade), mean)

agg$Group.1 = c("Arctic A\n(arctic fox)",
                "Asian SEA2a\n(dog)",
                "Asian SEA2b\n(CFB)",
                "Bat DR\n(vampire bat)",
                "Bat EF-E2\n(big brown bat)",
                "Bat LC\n(hoary bat)",
                "Bat TB1\n(Mexican free\n-tailed bat)",
                "Cosmo AF1b\n(dog)",
                "Cosmo AM2a\n(mongoose)",
                "RAC-SK SCSK\n(skunk)"
)

melt_data <- melt(agg, na.rm = FALSE, value.name = "rscu", id = "Group.1")

melt_data$Group.1 = factor(melt_data$Group.1, c("Bat LC\n(hoary bat)",
                                                "RAC-SK SCSK\n(skunk)",
                                                "Bat EF-E2\n(big brown bat)",
                                                "Bat DR\n(vampire bat)",
                                                "Bat TB1\n(Mexican free\n-tailed bat)",
                                                "Asian SEA2b\n(CFB)",
                                                "Asian SEA2a\n(dog)",
                                                "Arctic A\n(arctic fox)",
                                                "Cosmo AM2a\n(mongoose)",
                                                "Cosmo AF1b\n(dog)"))


p = ggplot(melt_data, aes(x = variable, y= Group.1, fill= rscu)) + 
  geom_tile() + xlab("Codon") + ylab("Clade") +
  scale_fill_continuous_divergingx(palette = 'RdBu', rev = T, mid = 1,
                                   l3 = 0, p3 = .8, p4 = .6,
                                   name = "RSCU") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  

p

png("plots/Figure 3.png", width = 9, height = 5, units = 'in', res = 600)
p
dev.off()
