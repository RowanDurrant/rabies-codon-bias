library(stringr)
ZAP_optimal = function(x){
  x = unname(as.character(x))
  location = c()
  location = append(location, unlist(gregexpr(pattern = "C........G.CG", x)))
  location = append(location, unlist(gregexpr(pattern = "C.......G.CG", x)))
  location = append(location, unlist(gregexpr(pattern = "C......G.CG", x)))
  location = append(location, unlist(gregexpr(pattern = "C.....G.CG", x)))
  location = append(location, unlist(gregexpr(pattern = "C....G.CG", x)))
  return(location)
}


library("Biostrings")

seqs = readDNAStringSet("sequence_data/all_seqs.fasta")
seqs = seqs[1:length(seqs)-1]
metadata = read.csv("sequence_data/metadata.csv")

locations = c()
accession = c()
clade = c()

for(j in 1:length(seqs)){
  locations = append(locations, ZAP_optimal(seqs[j]))
  accession = append(accession, rep(names(seqs[j]), length(ZAP_optimal(seqs[j]))))
  clade = append(clade, rep(metadata$Clade[metadata$Accession==names(seqs[j])], length(ZAP_optimal(seqs[j]))))
}

df = as.data.frame(cbind(accession, locations, clade))
df = df[df$locations != -1,]
df$clade = factor(df$clade, c("Cosmo AF1b", "Cosmo AM2a", "Arctic A", "Asian SEA2a", 
                              "Asian SEA2b", 
                              "Bat TB1",
                              "Bat DR", "Bat EF-E2","RAC-SK SCSK", "Bat LC"))

library(ggplot2)
ggplot(data = df, aes(y = clade, x = as.numeric(locations)))+
  geom_vline(xintercept = as.numeric(unique(df$locations)),
             colour = "lightgray", linewidth = 0.1)+
  geom_count(aes(color = clade, group = clade, size = ..prop..), alpha = 0.8) + theme_classic()+ 
  xlim(1,1353)+
  xlab("position of ZAP optimal motifs (C(n7)G(n)CG)") + ylab("Clade") +
  scale_color_manual(values = c("#332288","#88CCEE","#CCDDAA","#44AA99","#117733",  
                                "#999933", "#DDCC77","#CC6677","#882255","#AA4499"), 
                     name = "Clade", guide = guide_legend(),
                     limits = c("Cosmo AF1b", "Cosmo AM2a", "Arctic A", "Asian SEA2a", 
                                  "Asian SEA2b", 
                                  "Bat TB1",
                                  "Bat DR", "Bat EF-E2","RAC-SK SCSK", "Bat LC"),
                     labels = c("Cosmo AF1b\n(dog)",
                                "Cosmo AM2a\n(mongoose)",
                                "Arctic A\n(arctic fox)",
                                "Asian SEA2a\n(dog)",
                                "Asian SEA2b\n(CFB)",
                                "Bat TB1\n(Mexican free\n-tailed bat)",
                                "Bat DR\n(vampire bat)",
                                "Bat EF-E2\n(big brown bat)",
                                "RAC-SK SCSK\n(skunk)",
                                "Bat LC\n(hoary bat)"
                     )) +
  scale_y_discrete(limits= c("Cosmo AF1b", "Cosmo AM2a", "Arctic A", "Asian SEA2a", 
                             "Asian SEA2b", 
                             "Bat TB1",
                             "Bat DR", "Bat EF-E2","RAC-SK SCSK", "Bat LC"), 
                   labels = c("Cosmo AF1b\n(dog)",
                              "Cosmo AM2a\n(mongoose)",
                              "Arctic A\n(arctic fox)",
                              "Asian SEA2a\n(dog)",
                              "Asian SEA2b\n(CFB)",
                              "Bat TB1\n(Mexican free\n-tailed bat)",
                              "Bat DR\n(vampire bat)",
                              "Bat EF-E2\n(big brown bat)",
                              "RAC-SK SCSK\n(skunk)",
                              "Bat LC\n(hoary bat)"
  ))+
  guides(colour="none")+
  scale_size_continuous(breaks = c(1,25,100), name = "no. sequences")

