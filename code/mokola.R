library(stringr)
library(Biostrings)
library(ggplot2)
library(ggpubr)
library(viridis)

CpG = function(x){
  x = unname(as.character(x))
  nG = str_count(x, "G")
  nC = str_count(x, "C")
  nCpG = str_count(x, "CG")
  N = nchar(x)
  ObsExpCpG = (nCpG/(nC * nG))*N
  return(ObsExpCpG)
}

GCcontent = function(x){
  x = unname(as.character(x))
  nG = str_count(x, "G")
  nC = str_count(x, "C")
  N = nchar(x)
  gccont = (nC + nG)/N
  return(gccont)
  
}

CpG_actual = function(x){
  x = unname(as.character(x))
  nCpG = str_count(x, "CG")
  return(nCpG)
}

seqs = readDNAStringSet("mokola/mokola.fasta")
meta = read.csv("mokola/mokola_metadata.csv")

cpg = c()
gc = c()
cpg_actual = c()
accessions = c()
seq_length = c()

for(j in 1:length(seqs)){
  cpg = append(cpg, CpG(seqs[j]))
  gc = append(gc, GCcontent(seqs[j]))
  cpg_actual = append(cpg_actual, CpG_actual(seqs[j]))
  accessions = append(accessions, names(seqs[j]))
  seq_length = append(seq_length, seqs[j]@ranges@width)
}

df = as.data.frame(cbind(meta, cpg, gc, cpg_actual, seq_length))

p1 = ggplot(data = df, aes(x = Country, y = as.numeric(cpg)))+
  geom_violin(scale = "width", draw_quantiles = 0.5) +
  geom_jitter(width = 0.2, height = 0, aes(colour = as.numeric(seq_length))) +
  ylab("Obs/Exp CpG content") + theme_bw() + scale_colour_viridis_c(direction = -1,
                                                                    name = "Sequence length")

p2 = ggplot(data = df, aes(x = Host, y = as.numeric(cpg)))+
  geom_violin(scale = "width", draw_quantiles = 0.5)+
  geom_jitter(width = 0.2, height = 0, aes(colour = as.numeric(seq_length))) +
    ylab("Obs/Exp CpG content") + theme_bw() + scale_colour_viridis_c(direction = -1,
                                                                      name = "Sequence length")

ggarrange(p1,p2, common.legend = T, legend = "bottom")

# 
# plot(df$cpg~df$Year)+
#   abline(lm(df$cpg~df$Year))