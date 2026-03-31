library(uShuffleR)
shuffled = uShuffleR::shuffle("sequence_data/all_seqs.fasta", k = 2, n = 100)

writeXStringSet(shuffled, "sequence_data/reshuffled.fasta")
