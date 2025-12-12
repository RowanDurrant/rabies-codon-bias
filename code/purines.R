library(stringr)
library(seqinr)

seqs = read.fasta("sequence_data/all_seqs.fasta",as.string = T)
seqs = seqs[1:length(seqs)-1]
metadata = read.csv("sequence_data/metadata.csv")
metadata$purines = NA
metadata$purines1 = NA
metadata$purines2 = NA
metadata$purines3 = NA
for(i in 1:length(seqs)){
  split = str_split(as.character(seqs[i]), "")[[1]]
  metadata$purines[i] = length(split[split == "a" | split == "g"])/length(split)
  third_only = split[(1:(length(split)/3))*3]
  metadata$purines3[i] = length(third_only[third_only == "a" | third_only == "g"])/length(third_only)
  second_only = split[((1:(length(split)/3))*3)-1]
  metadata$purines2[i] = length(second_only[second_only == "a" | second_only == "g"])/length(second_only)
  first_only = split[((1:(length(split)/3))*3)-2]
  metadata$purines1[i] = length(first_only[first_only == "a" | first_only == "g"])/length(first_only)
}

write.csv(metadata, "output_data/purines.csv")
