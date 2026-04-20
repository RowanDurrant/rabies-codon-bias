set.seed(84298913)
library(ape)
library(ggplot2)
library(ggpubr)
library(ggh4x)
library(boot)
library(readxl)
library(MCMCglmm)

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
cc = c()
for(i in c("PC1", "PC2", "PC3")){
  for(j in c("cpg","cpg_actual","gc","tpa","tpa_actual","ta","ZAP_optimal_motifs",
             "ZAP_suboptimal_motifs", colnames(df)[3:43])){
    pc = append(pc, i)
    comp = append(comp,j)
    cc = append(cc, cor(df[,j]/100, df[,i]))
  }
}

r_df = as.data.frame(cbind(pc,comp,cc))
r_df$cc_abs = abs(as.numeric(r_df$cc))
