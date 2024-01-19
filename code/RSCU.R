library("readxl")
df = as.data.frame(read_excel("data/RSCU_N.xlsx"))
for(i in 3:ncol(df)){
  colnames(df)[i] = paste0(as.character(df[1,i]),"_",colnames(df)[i])
  df[,i] = as.numeric(df[,i])
}
colnames(df)[2] = "Accession no."
df = df[2:nrow(df),]

agg = aggregate(df[, 3:61], list(df$Host), mean)

agg$Group.1 = c("Bat EF-E2\n(big brown bat)","Asian SEA2b\n(CFB)",
                "Cosmo AF1b\n(dog)","Asian SEA2a\n(dog)",
                "Bat TB1\n(Mexican free-tailed bat)","Bat LC\n(hoary bat)",
                "Cosmo AM2a\n(mongoose)",
                "RAC-SK SCSK\n(skunk)","Bat DR\n(vampire bat)"
)

library(reshape2)
melt_data <- melt(agg, na.rm = FALSE, value.name = "rscu", id = "Group.1")

melt_data$Group.1 = factor(melt_data$Group.1, c("Bat LC\n(hoary bat)",
                                                "RAC-SK SCSK\n(skunk)",
                                                "Bat EF-E2\n(big brown bat)",
                                                "Bat DR\n(vampire bat)",
                                                "Bat TB1\n(Mexican free-tailed bat)",
                                                "Asian SEA2b\n(CFB)",
                                                "Asian SEA2a\n(dog)",
                                                "Cosmo AM2a\n(mongoose)",
                                                "Cosmo AF1b\n(dog)"))

library(ggplot2)
library(RColorBrewer)
library(colorspace)
ggplot(melt_data, aes(x = variable, y= Group.1, fill= rscu)) + 
  geom_tile() + xlab("Codon") + ylab("Clade") +
  scale_fill_continuous_divergingx(palette = 'RdBu', rev = T, mid = 1,
                                   l3 = 0, p3 = .8, p4 = .6,
                                   name = "RSCU") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
