library(ggplot2)
library(ggpubr)
library(ggh4x)
library(boot)

df_pca = read.csv("output_data/PCA_output_RSCU.csv")
df_cpg = read.csv("output_data/N_CpG.csv")
df_purine = read_xlsx("output_data/Nucleotide_composition_N.xlsx")
colnames(df_purine)[1] = "Accession"

df = merge(df_pca, df_cpg, by.x = c("X", "clade"), by.y = c("accessions", "clade"))
df = merge(df_purine,df, by.x = "Accession", by.y = "X")

df$clade = factor(df$clade, c("Cosmopolitan AF1b", "Cosmopolitan AM2a", "Arctic A", "Asian SEA2a", 
                              "Asian SEA2b", 
                              "Bats DR","Bats TB1", "Bats LC",
                              "Bats EF-E2","RAC-SK SCSK"))

pc = c()
comp = c()
r2 = c()
for(i in c("PC1", "PC2", "PC3")){
  for(j in c("cpg","cpg_actual","gc","tpa","tpa_actual","ta","ZAP_optimal_motifs",
             "ZAP_suboptimal_motifs", colnames(df)[3:59])){
    pc = append(pc, i)
    comp = append(comp,j)
    r2 = append(r2, summary(lm(df[,i]~df[,j]))$adj.r.squared)
  }
}

r_df = as.data.frame(cbind(pc,comp,r2))
