library(ape)
library(Biostrings)

G_gene = readDNAStringSet("sequences/G_gene/AL_Bats_EF-E2_G_nucleotide_alignment.fasta")
N_gene = readDNAStringSet("sequences/N_gene/Eptesicus fuscus_Bats_EF-E2_N_nucleotide_alignment.fasta")

concatenated = DNAStringSet()
all_G = DNAStringSet()

all_G = c(all_G, G_gene)
for(i in 1:length(N_gene)){
  if(names(N_gene[i]) %in% names(G_gene)){
    conc = xscat(N_gene[i], G_gene[names(G_gene) == names(N_gene[i])])
    names(conc) = names(N_gene[i])
    if(length(grep("N", conc)) == 0 & length(grep("R", conc)) == 0 &
       length(grep("Y", conc)) == 0){
      concatenated = c(concatenated, conc)
      print("BBB")
    }
  }
  
}

G_gene = readDNAStringSet("sequences/G_gene/AL_Bats_LC_G_nucleotide_alignment.fasta")
N_gene = readDNAStringSet("sequences/N_gene/Hoary_bat_Bats_LC_N_nucleotide_alignment.fasta")
all_G = c(all_G, G_gene)
for(i in 1:length(N_gene)){
  if(names(N_gene[i]) %in% names(G_gene)){
    conc = xscat(N_gene[i], G_gene[names(G_gene) == names(N_gene[i])])
    names(conc) = names(N_gene[i])
    if(length(grep("N", conc)) == 0 & length(grep("R", conc)) == 0 &
       length(grep("Y", conc)) == 0){
      concatenated = c(concatenated, conc)
      print("Hoary")
    }
  }
  
}

G_gene = readDNAStringSet("sequences/G_gene/AL_Bats_TB1_G_nucleotide_alignment.fasta")
N_gene = readDNAStringSet("sequences/N_gene/Tadarida brasiliensis_Bats_TB1_N_nucleotide_alignment.fasta")
all_G = c(all_G, G_gene)
for(i in 1:length(N_gene)){
  if(names(N_gene[i]) %in% names(G_gene)){
    conc = xscat(N_gene[i], G_gene[names(G_gene) == names(N_gene[i])])
    names(conc) = names(N_gene[i])
    if(length(grep("N", conc)) == 0 & length(grep("R", conc)) == 0 &
       length(grep("Y", conc)) == 0){
      concatenated = c(concatenated, conc)
      print("Freetailed")
    }
  }
  
}

G_gene = readDNAStringSet("sequences/G_gene/cfb_Asian_SEA2b_G_nucleotide_alignment.fasta")
N_gene = readDNAStringSet("sequences/N_gene/CFB_SEA2b_N_nucleotide_alignment.fasta")
all_G = c(all_G, G_gene)
for(i in 1:length(N_gene)){
  if(names(N_gene[i]) %in% names(G_gene)){
    conc = xscat(N_gene[i], G_gene[names(G_gene) == names(N_gene[i])])
    names(conc) = names(N_gene[i])
    if(length(grep("N", conc)) == 0 & length(grep("R", conc)) == 0 &
       length(grep("Y", conc)) == 0){
      concatenated = c(concatenated, conc)
      print("CFB")
    }
  }
  
}

G_gene = readDNAStringSet("sequences/G_gene/Desmodus_rotundus_Bats_DR_G_nucleotide_alignment.fasta")
N_gene = readDNAStringSet("sequences/N_gene/Desmodus_rotundus_Bats_DR_N_nucleotide_alignment.fasta")
all_G = c(all_G, G_gene)
for(i in 1:length(N_gene)){
  if(names(N_gene[i]) %in% names(G_gene)){
    conc = xscat(N_gene[i], G_gene[names(G_gene) == names(N_gene[i])])
    names(conc) = names(N_gene[i])
    if(length(grep("N", conc)) == 0 & length(grep("R", conc)) == 0 &
       length(grep("Y", conc)) == 0){
      concatenated = c(concatenated, conc)
      print("Vamp")
    }
  }
  
}

G_gene = readDNAStringSet("sequences/G_gene/dog_Asian_SEA2a_G_nucleotide_alignment.fasta")
N_gene = readDNAStringSet("sequences/N_gene/Canis_familiaris_SEA2a_N_nucleotide_alignment.fasta")
all_G = c(all_G, G_gene)
for(i in 1:length(N_gene)){
  if(names(N_gene[i]) %in% names(G_gene)){
    conc = xscat(N_gene[i], G_gene[names(G_gene) == names(N_gene[i])])
    names(conc) = names(N_gene[i])
    if(length(grep("N", conc)) == 0 & length(grep("R", conc)) == 0 &
       length(grep("Y", conc)) == 0){
      concatenated = c(concatenated, conc)
      print("SEA2a")
    }
  }
  
}

G_gene = readDNAStringSet("sequences/G_gene/dog_Cosmopolitan_AF1b_G_nucleotide_alignment.fasta")
N_gene = readDNAStringSet("sequences/N_gene/Dog_AF1b_N_nucleotide_alignment.fasta")
all_G = c(all_G, G_gene)
for(i in 1:length(N_gene)){
  if(names(N_gene[i]) %in% names(G_gene)){
    conc = xscat(N_gene[i], G_gene[names(G_gene) == names(N_gene[i])])
    names(conc) = names(N_gene[i])
    if(length(grep("N", conc)) == 0 & length(grep("R", conc)) == 0 &
       length(grep("Y", conc)) == 0){
      concatenated = c(concatenated, conc)
      print("AF1b")
    }
  }
  
}

G_gene = readDNAStringSet("sequences/G_gene/scsc_RAC-SK_G_nucleotide_alignment.fasta")
N_gene = readDNAStringSet("sequences/N_gene/Skunk-SK_N_nucleotide_alignment.fasta")

for(i in 1:length(N_gene)){
  if(names(N_gene[i]) %in% names(G_gene)){
    conc = xscat(N_gene[i], G_gene[names(G_gene) == names(N_gene[i])])
    names(conc) = names(N_gene[i])
    if(length(grep("N", conc)) == 0 & length(grep("R", conc)) == 0 &
       length(grep("Y", conc)) == 0){
      concatenated = c(concatenated, conc)
      all_G = c(all_G, G_gene[names(G_gene) == names(N_gene[i])])
      print("Skunk")
    }
  }
  
}

concatenated_unique = unique(concatenated)
writeXStringSet(concatenated_unique, "sequences/N_G_conc.fasta")
all_G_unique = unique(all_G)
writeXStringSet(all_G_unique, "sequences/G_gene/all_G_unique.fasta")

all_N = readDNAStringSet("sequences/N_gene/all_seqs.fasta")
all_N_unique = unique(all_N)
writeXStringSet(all_N_unique, "sequences/N_gene/all_N_unique.fasta")
