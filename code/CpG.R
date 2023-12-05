##CpG 
# Observed to Expected CpG is calculated as below : 
# Obs/Exp CpG = Number of CpG * N / (Number of C * Number of G) 
# where N = length of sequence.

library(stringr)

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

library(Biostrings)

 files = c("sequences/N_gene/Dog_AF1b_N_nucleotide_alignment.fasta",
           "sequences/N_gene/Canis_familiaris_SEA2a_N_nucleotide_alignment.fasta",
           "sequences/N_gene/CFB_SEA2b_N_nucleotide_alignment.fasta",
           "sequences/N_gene/Desmodus_rotundus_Bats_DR_N_nucleotide_alignment.fasta",
           "sequences/N_gene/Eptesicus fuscus_Bats_EF-E2_N_nucleotide_alignment.fasta",
           "sequences/N_gene/Hoary_bat_Bats_LC_N_nucleotide_alignment.fasta",
           "sequences/N_gene/Mongoose_AM2a_N_nucleotide_alignment.fasta",
           "sequences/N_gene/Skunk-SK_N_nucleotide_alignment.fasta",
           "sequences/N_gene/Tadarida brasiliensis_Bats_TB1_N_nucleotide_alignment.fasta")

 Host = c("Dog (AF1b)","Dog (SEA2a)", "Chinese ferret badger","Vampire bat","Big brown bat",
          "Hoary bat", "Mongoose", "Skunk", "Free-tailed bat")

# files = c("sequences/WGS/BBB_Bats_EF-E2_whole_genome_nucleotide_alignment.fasta",
#           "sequences/WGS/CFB_Asian_SEA2b_whole_genome_nucleotide_alignment.fasta",
#           "sequences/WGS/dog_Asian_SEA2a_whole_genome_nucleotide_alignment.fasta",
#           "sequences/WGS/dog_Cosmopolitan_AF1b_whole_genome_nucleotide_alignment.fasta",
#           "sequences/WGS/free-tailed_Bats_TB1_whole_genome_nucleotide_alignment.fasta",
#           "sequences/WGS/hoary_Bats_LC_whole_genome_nucleotide_alignment.fasta",
#           "sequences/WGS/skunk_RAC-SK_whole_genome_nucleotide_alignment.fasta")
#Host = c("Big brown bat","Chinese ferret badger","Dog (SEA2a)","Dog (AF1b)", 
#         "Free-tailed bat","Hoary bat", "Skunk")
 
 # files = c("sequences/G_gene/dog_Cosmopolitan_AF1b_G_nucleotide_alignment.fasta",
 #           "sequences/G_gene/dog_Asian_SEA2a_G_nucleotide_alignment.fasta",
 #           "sequences/G_gene/cfb_Asian_SEA2b_G_nucleotide_alignment.fasta",
 #           "sequences/G_gene/Desmodus_rotundus_Bats_DR_G_nucleotide_alignment.fasta",
 #           "sequences/G_gene/AL_Bats_EF-E2_G_nucleotide_alignment.fasta",
 #           "sequences/G_gene/AL_Bats_LC_G_nucleotide_alignment.fasta",
 #           "sequences/G_gene/scsc_RAC-SK_G_nucleotide_alignment.fasta",
 #           "sequences/G_gene/AL_Bats_TB1_G_nucleotide_alignment.fasta")
 # 
 # Host = c("Dog (AF1b)","Dog (SEA2a)", "Chinese ferret badger","Vampire bat","Big brown bat",
 #          "Hoary bat", "Skunk", "Free-tailed bat")

cpg = c()
gc = c()
cpg_actual = c()
accessions = c()
hosts = c()

for(i in 1:length(files)){
  seqs = readDNAStringSet(files[i])
  for(j in 1:length(seqs)){
       cpg = append(cpg, CpG(seqs[j]))
       gc = append(gc, GCcontent(seqs[j]))
       cpg_actual = append(cpg_actual, CpG_actual(seqs[j]))
       accessions = append(accessions, names(seqs[j]))
       hosts = append(hosts, Host[i])
       
    }
  
}
df = data.frame(accessions, hosts, cpg, gc, cpg_actual)
write.csv(df, "N_CpG.csv")
df$hosts = factor(df$hosts, c("Big brown bat", "Hoary bat", "Free-tailed bat",
                              "Vampire bat", "Dog (SEA2a)", "Chinese ferret badger",
                              "Dog (AF1b)", "Mongoose", "Skunk"))

df$host_group = NA
df$host_group[df$hosts %in% c("Big brown bat", "Hoary bat", "Free-tailed bat",
                              "Vampire bat")] = "Bats"
df$host_group[df$hosts %in% c("Dog (SEA2a)", "Chinese ferret badger",
                              "Dog (AF1b)", "Mongoose", "Skunk")] = "Carnivores"

my_pal <- c("#b66dff","#490092","#006ddb","#6db6ff","#b6dbff",
            "#004949","#009292","#ff6db6","#ffb6db"
            
)
#write.csv(df, "cpg.csv")

library(ggplot2)
p1= ggplot(data = df, aes(x = hosts, y = cpg))+
  geom_boxplot()+ 
  geom_jitter(aes(color = hosts), alpha  = 0.4, 
              width = 0.4, height = 0) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")+
  ylab("Obs/Exp CpG") + xlab("Clade") +
  scale_color_manual(values = c("#004949","#009292","#ff6db6","#ffb6db",
                                "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
                                "#920000","#924900"), name = "Clade", guide = guide_legend(),
                     labels = c("Bat EF-E2\n(big brown bat)",
                                "Bat LC\n(hoary bat)",
                                "Bat TB1\n(Mexican free-tailed bat)",
                                "Bat DR\n(vampire bat)",
                                "Asian SEA2a\n(dog)",
                                "Asian SEA2b\n(CFB)",
                                "Cosmo AF1b\n(dog)",
                                "Cosmo AM2a\n(mongoose)",
                                "RAC-SK SCSC\n(skunk)",
                                "Gannoruwa bat lyssavirus"
                     )) +
  scale_x_discrete(labels = c("Bat EF-E2\n(big brown bat)",
                              "Bat LC\n(hoary bat)",
                              "Bat TB1\n(Mexican free-tailed bat)",
                              "Bat DR\n(vampire bat)",
                              "Asian SEA2a\n(dog)",
                              "Asian SEA2b\n(CFB)",
                              "Cosmo AF1b\n(dog)",
                              "Cosmo AM2a\n(mongoose)",
                              "RAC-SK SCSC\n(skunk)",
                              "Gannoruwa bat lyssavirus"
  ))

p2= ggplot(data = df, aes(x = hosts, y = gc))+
  geom_boxplot()+ 
  geom_jitter(aes(color = hosts), alpha  = 0.4, 
              width = 0.4, height = 0) + theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")+
  ylab("GC content") + xlab("Clade") +
  scale_color_manual(values = c("#004949","#009292","#ff6db6","#ffb6db",
                                "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
                                "#920000","#924900"), name = "Clade", guide = guide_legend(),
                     labels = c("Bat EF-E2\n(big brown bat)",
                                "Bat LC\n(hoary bat)",
                                "Bat TB1\n(Mexican free-tailed bat)",
                                "Bat DR\n(vampire bat)",
                                "Asian SEA2a\n(dog)",
                                "Asian SEA2b\n(CFB)",
                                "Cosmo AF1b\n(dog)",
                                "Cosmo AM2a\n(mongoose)",
                                "RAC-SK SCSC\n(skunk)",
                                "Gannoruwa bat lyssavirus"
                     )) +
  scale_x_discrete(labels = c("Bat EF-E2\n(big brown bat)",
                              "Bat LC\n(hoary bat)",
                              "Bat TB1\n(Mexican free-tailed bat)",
                              "Bat DR\n(vampire bat)",
                              "Asian SEA2a\n(dog)",
                              "Asian SEA2b\n(CFB)",
                              "Cosmo AF1b\n(dog)",
                              "Cosmo AM2a\n(mongoose)",
                              "RAC-SK SCSC\n(skunk)",
                              "Gannoruwa bat lyssavirus"
  )) 

p3 = ggplot(data = df, aes(x = hosts, y = cpg_actual))+
  geom_boxplot()+ 
  geom_jitter(aes(color = hosts), alpha  = 0.4, 
              width = 0.4, height = 0) + theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")+
  ylab("No. CpGs") + xlab("Clade") +
  scale_color_manual(values = c("#004949","#009292","#ff6db6","#ffb6db",
                                "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
                                "#920000","#924900"), name = "Clade", guide = guide_legend(),
                     labels = c("Bat EF-E2\n(big brown bat)",
                                "Bat LC\n(hoary bat)",
                                "Bat TB1\n(Mexican free-tailed bat)",
                                "Bat DR\n(vampire bat)",
                                "Asian SEA2a\n(dog)",
                                "Asian SEA2b\n(CFB)",
                                "Cosmo AF1b\n(dog)",
                                "Cosmo AM2a\n(mongoose)",
                                "RAC-SK SCSC\n(skunk)",
                                "Gannoruwa bat lyssavirus"
                     )) +
  scale_x_discrete(labels = c("Bat EF-E2\n(big brown bat)",
                              "Bat LC\n(hoary bat)",
                              "Bat TB1\n(Mexican free-tailed bat)",
                              "Bat DR\n(vampire bat)",
                              "Asian SEA2a\n(dog)",
                              "Asian SEA2b\n(CFB)",
                              "Cosmo AF1b\n(dog)",
                              "Cosmo AM2a\n(mongoose)",
                              "RAC-SK SCSC\n(skunk)",
                              "Gannoruwa bat lyssavirus"
  ))
library(ggpubr)
ggarrange(p1, ggarrange(p2, p3, labels = c("B", "C"), ncol = 1, nrow = 2), labels = c("A"))

# t.test(df$cpg[df$host_group == "Bats"], df$cpg[df$host_group == "Carnivores"])
# t.test(df$gc[df$host_group == "Bats"], df$gc[df$host_group == "Carnivores"])
# t.test(df$cpg_actual[df$host_group == "Bats"], df$cpg_actual[df$host_group == "Carnivores"])
# 
# SDUc = read.table("myseqs_SDUc_bootstrapped.tsv", sep = "\t", header = T)
# 
# df$CpGpos1 = NA
# df$CpGpos2 = NA
# df$CpGbridge = NA
# df$CpGpos1_95 = NA
# df$CpGpos2_95 = NA
# df$CpGbridge_95 = NA
# for(i in df$accessions){
#   df$CpGpos1[df$accessions == i] = SDUc$CpGpos1[SDUc$acc == i]
#   df$CpGpos2[df$accessions == i] = SDUc$CpGpos2[SDUc$acc == i]
#   df$CpGbridge[df$accessions == i] = SDUc$CpGbridge[SDUc$acc == i]
#  
#   df$CpGpos1_95[df$accessions == i] = (SDUc$CpGpos1[SDUc$acc == i] > SDUc$CpGpos1_low95CI[SDUc$acc == i] &
#                                          SDUc$CpGpos1[SDUc$acc == i] < SDUc$CpGpos1_high95CI[SDUc$acc == i])
#   df$CpGpos2_95[df$accessions == i] = (SDUc$CpGpos2[SDUc$acc == i] > SDUc$CpGpos2_low95CI[SDUc$acc == i] &
#                                          SDUc$CpGpos2[SDUc$acc == i] < SDUc$CpGpos2_high95CI[SDUc$acc == i])
#   df$CpGbridge_95[df$accessions == i] = (SDUc$CpGbridge[SDUc$acc == i] > SDUc$CpGbridge_low95CI[SDUc$acc == i] &
#                                            SDUc$CpGbridge[SDUc$acc == i] < SDUc$CpGbridge_high95CI[SDUc$acc == i])
#   
# }
# p4= ggplot(data = df, aes(x = hosts, y = CpGpos1))+
#   geom_boxplot()+ theme_bw()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ylab("CpG position 1") + xlab("Host") 
# 
# p5= ggplot(data = df, aes(x = hosts, y = CpGpos2))+
#   geom_boxplot()+ theme_bw()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ylab("CpG position 2") + xlab("Host") 
# 
# p6 = ggplot(data = df, aes(x = hosts, y = CpGbridge))+
#   geom_boxplot()+ theme_bw()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ylab("CpG bridge") + xlab("Host") 
# 
# library(ggpubr)
# ggarrange(p4, p5, p6, labels = c("A", "B", "C"), nrow = 1)
