#plot idea stolen from https://www.biorxiv.org/content/10.1101/2023.01.23.525187v1.full

tab=read.table('trees/Hoary_bat_Bats_LC_N_nucleotide_alignment.fasta.state',header=TRUE)
seqs = readDNAStringSet("sequences/Hoary_bat_Bats_LC_N_nucleotide_alignment.fasta")

nodes = c("Node1", "Node2", "Node3", "Node4", "Node5", "Node6", "Node7", "Node8",
          "Node9", "Node10", "Node11", "Node12", "Node13", "Node14")
descendents = list(c(paste0(tab$State[tab$Node == "Node2"], collapse = ""), #node 1
                     as.character(seqs$AF351845)),
                   c(paste0(tab$State[tab$Node == "Node3"], collapse = ""), #node 2
                     as.character(seqs$JQ685947)),
                   c(paste0(tab$State[tab$Node == "Node4"], collapse = ""), #node 3
                     paste0(tab$State[tab$Node == "Node9"], collapse = "")),
                   c(paste0(tab$State[tab$Node == "Node5"], collapse = ""), #node 4
                     paste0(tab$State[tab$Node == "Node8"], collapse = "")),
                   c(paste0(tab$State[tab$Node == "Node6"], collapse = ""), #node5
                     paste0(tab$State[tab$Node == "Node7"], collapse = "")),
                   c(as.character(seqs$AF351846), #node6
                     as.character(seqs$GU644717)),
                   c(as.character(seqs$GU644715), #node7
                     as.character(seqs$GU644721)),
                   c(as.character(seqs$GU644712), #node8
                     as.character(seqs$GU644716)),
                   c(paste0(tab$State[tab$Node == "Node10"], collapse = ""), #node9
                     as.character(seqs$GU644714)),
                   c(as.character(seqs$AF351858), #node10
                     as.character(seqs$GU644718)),
                   c(paste0(tab$State[tab$Node == "Node12"], collapse = ""), #node11
                     paste0(tab$State[tab$Node == "Node13"], collapse = "")),
                   c(as.character(seqs$AF394883), #node2
                     as.character(seqs$AF394884)),
                   c(paste0(tab$State[tab$Node == "Node14"], collapse = ""), #node13
                     as.character(seqs$GU644719)),
                   c(as.character(seqs$GU644713), #node14
                     as.character(seqs$GU644720))
                   )

position = c()
nucleotide = c()

for(i in 1:length(nodes)){
  Nodetab = tab[tab$Node == nodes[i],]
  for(k in 1:length(descendents[[i]])){
    descendent = descendents[[i]][k]
    for(j in 1:nrow(Nodetab)){
      if(Nodetab$State[j] == "C"){
        if(substring(descendent, j, j) == "T" & j >= 4){
          position = append(position, c("-3", "-2", "-1", "0", "1", "2", "3"))
          nucleotide = append(nucleotide, c(Nodetab$State[j-3], Nodetab$State[j-2],
                                            Nodetab$State[j-1], Nodetab$State[j],
                                            Nodetab$State[j+1], Nodetab$State[j+2],
                                            Nodetab$State[j+3]))
        
      }
    }
    
  }
    
    }
  }

df = data.frame(position, nucleotide)
df$position = factor(df$position, levels = c("-3", "-2", "-1", "0", "1", "2", "3"))

library(dplyr)
ggplot(df %>% count(position, nucleotide),  # Calculate label positions
       aes(position, n*100/17, fill=nucleotide)) +
  geom_bar(stat="identity") + ylim(0,100) + ylab("%")+
  theme_bw()
