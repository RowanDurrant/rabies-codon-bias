dinuc = function(x, di){
  x = unname(as.character(x))
  n1 = str_count(x, str_split(di, "")[[1]][1])
  n2 = str_count(x, str_split(di, "")[[1]][2])
  ndi = str_count(x, di)
  N = nchar(x)
  ObsExp = (ndi/(n1 * n2))*N
  return(ObsExp)
}

dinucleotides = c("AA", "AC", "AG", "AT",
                  "CA", "CC", "CG", "CT",
                  "GA", "GC", "GG", "GT",
                  "TA", "TC", "TG", "TT")

library("Biostrings")

seqs = readDNAStringSet("sequence_data/all_seqs.fasta")
seqs = seqs[1:length(seqs)-1]
metadata = read.csv("sequence_data/metadata.csv")

accessions = c()
clade = c()
host_group = c()
dinuc_list = list()

for(j in 1:length(seqs)){
  for(i in dinucleotides){
    dinuc_list[[i]][j] = dinuc(seqs[j], i)
  }

  accessions = append(accessions, names(seqs[j]))
  clade = append(clade, metadata$Clade[metadata$Accession==names(seqs[j])])
  host_group = append(host_group, metadata$Group[metadata$Accession==names(seqs[j])])

}

library(reshape2)
df = melt(dinuc_list)
ggplot(df,aes(x = L1, y = value)) + 
  geom_jitter(alpha = 0.3, height = 0)+
  theme_bw()
