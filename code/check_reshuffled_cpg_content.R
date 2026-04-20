library(seqinr)
library(stringr)

CpG = function(x){
  x = unname(as.character(x))
  nG = str_count(x, "g")
  nC = str_count(x, "c")
  nCpG = str_count(x, "cg")
  N = nchar(x)
  ObsExpCpG = (nCpG/(nC * nG))*N
  return(ObsExpCpG)
}

reshuffled = seqinr::read.fasta("sequence_data/reshuffled.fasta", as.string= T)
reshuffled_AY352488 = reshuffled[1:100]
seqs = seqinr::read.fasta("sequence_data/all_seqs.fasta", as.string = T)

CpG(seqs["AY352488"])
cpg = c()
for(i in 1:length(reshuffled_AY352488)){
  cpg = append(cpg, CpG(reshuffled_AY352488[i]))
}
unique(cpg)
