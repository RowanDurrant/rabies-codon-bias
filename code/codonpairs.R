install.packages("remotes")
remotes::install_github("alex-sbu/CPBias")
library(Biostrings)
library("CPBias")
library("readxl")

seqs = readDNAStringSet("sequences/N_gene/all_seqs.fasta")
seqs = seqs[1:413]
ref = read.csv("data/N_CpG.csv") #for matching accessions to clades

ID = c()
codonpair = c()
cps = c()
host = c()

for(i in 1:length(seqs)){
  bias = CPBtable(seqs[[i]])
  codonpair = append(codonpair, bias[[1]]$`codon pair`)
  cps = append(cps, bias[[1]]$cps)
  ID = append(ID, rep(seqs@ranges@NAMES[i], length(bias[[1]]$`codon pair`)))
  host = append(host, rep(ref$hosts[ref$accessions == seqs@ranges@NAMES[i]], 
                          length(bias[[1]]$`codon pair`)))
}

df = data.frame(ID, codonpair, cps, host)

library(tidyr)
df2 = drop_na(df)

df2$host = factor(df2$host, c("Dog (AF1b)", "Mongoose","Dog (SEA2a)", "Chinese ferret badger",
                              "Free-tailed bat",
                              "Vampire bat", "Big brown bat","Skunk", "Hoary bat" 
))

library(ggplot2)
library(RColorBrewer)
library(colorspace)
p = ggplot(df2, aes(x = codonpair, y= host, fill= cps)) + 
  geom_tile() + xlab("Codon pair") + ylab("Clade") +
  scale_fill_continuous_divergingx(palette = 'RdBu', rev = T, mid = 1,
                                   l3 = 0, p3 = .8, p4 = .6,
                                   name = "CPS") +
  scale_y_discrete(labels = c("Cosmo AF1b\n(dog)",
                              "Cosmo AM2a\n(mongoose)",
                              "Asian SEA2a\n(dog)",
                              "Asian SEA2b\n(CFB)",
                              "Bat TB1\n(Mexican free\n-tailed bat)",
                              "Bat DR\n(vampire bat)",
                              "Bat EF-E2\n(big brown bat)",
                              "RAC-SK SCSK\n(skunk)",
                              "Bat LC\n(hoary bat)"
  ))+
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p

mean(df2$cps[df2$host %in% c("Dog (AF1b)", "Dog (SEA2a)") &
               df2$cps > 0])
mean(df2$cps[df2$host %in% c("Dog (AF1b)", "Dog (SEA2a)") &
               df2$cps < 0])

mean(df2$cps[df2$host %in% c("Free-tailed bat",
                                 "Vampire bat", "Big brown bat", 
                                 "Hoary bat") &
               df2$cps > 0])
mean(df2$cps[df2$host %in% c("Free-tailed bat",
                             "Vampire bat", "Big brown bat", 
                             "Hoary bat") &
               df2$cps < 0])

ggplot(data = df2, aes(x = host, y = cps))+
  geom_jitter(alpha = 0.3)

hostspecies = unique(df2$host)
nopairs = c()
for(i in hostspecies){
  nopairs = append(nopairs, length(unique(df2$codonpair[df2$host == i])))
}

barplot(nopairs ~ hostspecies,
        ylim = c(0,1000))
