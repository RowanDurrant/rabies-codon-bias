library(stringr)
library(Biostrings)

seqs = readDNAStringSet("sequence_data/all_seqs.fasta")
seqs = seqs[1:length(seqs)-1]
metadata = read.csv("sequence_data/metadata.csv")
metadata$pyrimidines = NA
metadata$pyrimidines1 = NA
metadata$pyrimidines2 = NA
metadata$pyrimidines3 = NA
for(i in 1:length(seqs)){
  split = str_split(as.character(seqs[i]), "")[[1]]
  metadata$pyrimidines[i] = length(split[split == "C" | split == "T"])/length(split)
  third_only = split[(1:(length(split)/3))*3]
  metadata$pyrimidines3[i] = length(third_only[third_only == "C" | third_only == "T"])/length(third_only)
  second_only = split[((1:(length(split)/3))*3)-1]
  metadata$pyrimidines2[i] = length(second_only[second_only == "C" | second_only == "T"])/length(second_only)
  first_only = split[((1:(length(split)/3))*3)-2]
  metadata$pyrimidines1[i] = length(first_only[first_only == "C" | first_only == "T"])/length(first_only)
}

write.csv(metadata, "output_data/purines.csv")
