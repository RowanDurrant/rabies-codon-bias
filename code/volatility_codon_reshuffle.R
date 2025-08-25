##Codon volatility

## v(c) = (i =1:9) Σ d(acid(c.i), acid(c))
#d = amino acid distance measure: Hamming or Miyata
#Hamming much simpler- probably just use this for now
library(readxl)
source("code/codon_reshuffle.R")

`%!in%` = Negate(`%in%`)
nucleotides = c("A", "T", "C", "G")
stop_codons = c("TGA", "TAA", "TAG")
df = as.data.frame(read_excel("output_data/Codon_usage_N.xlsx"))
volatility_df = read.csv("output_data/codon_volatility_stop_punished.csv")
volatility_df_0 = read.csv("output_data/codon_volatility.csv")
codons = colnames(df)[2:ncol(df)]

df$mean_volatility_h = NA
df$mean_volatility_m = NA
df$normalised_volatility_h = NA
df$normalised_volatility_m = NA
df$mean_volatility_0_h = NA
df$mean_volatility_0_m = NA
df$normalised_volatility_0_h = NA
df$normalised_volatility_0_m = NA

for(a in 2:nrow(df)){
  total_volatility_h = 0
  total_volatility_m = 0
  total_volatility_0_h = 0
  total_volatility_0_m = 0
  total_volatility_h_norm = 0
  total_volatility_m_norm = 0
  total_volatility_0_h_norm = 0
  total_volatility_0_m_norm = 0
  
  codon_count = 0
  reshuffled = codon_reshuffle(df[a,2:62])
  for(b in codons){
    total_volatility_h = total_volatility_h + (as.numeric(df[a, b])*as.numeric(volatility_df$hamming[volatility_df$codons == b]))
    total_volatility_m = total_volatility_m + (as.numeric(df[a, b])*as.numeric(volatility_df$miyata[volatility_df$codons == b]))
    total_volatility_0_h = total_volatility_0_h + (as.numeric(df[a, b])*as.numeric(volatility_df_0$hamming[volatility_df_0$codons == b]))
    total_volatility_0_m = total_volatility_0_m + (as.numeric(df[a, b])*as.numeric(volatility_df_0$miyata[volatility_df_0$codons == b]))
    total_volatility_h_norm = total_volatility_h_norm + (as.numeric(reshuffled[1, b])*as.numeric(volatility_df$hamming[volatility_df$codons == b]))
    total_volatility_m_norm = total_volatility_m_norm + (as.numeric(reshuffled[1, b])*as.numeric(volatility_df$miyata[volatility_df$codons == b]))
    total_volatility_0_h_norm = total_volatility_0_h_norm + (as.numeric(reshuffled[1, b])*as.numeric(volatility_df_0$hamming[volatility_df_0$codons == b]))
    total_volatility_0_m_norm = total_volatility_0_m_norm + (as.numeric(reshuffled[1, b])*as.numeric(volatility_df_0$miyata[volatility_df_0$codons == b]))
    
    codon_count = codon_count + as.numeric(df[a, b])
  }
  
  df$mean_volatility_h[a] = total_volatility_h/codon_count
  df$mean_volatility_m[a] = total_volatility_m/codon_count
  df$mean_volatility_0_h[a] = total_volatility_0_h/codon_count
  df$mean_volatility_0_m[a] = total_volatility_0_m/codon_count
  df$normalised_volatility_h[a] = total_volatility_h_norm/codon_count
  df$normalised_volatility_m[a] = total_volatility_m_norm/codon_count
  df$normalised_volatility_0_h[a] = total_volatility_0_h_norm/codon_count
  df$normalised_volatility_0_m[a] = total_volatility_0_m_norm/codon_count
}  

df2 = data.frame(cbind(df$CODONS, 
                       df$mean_volatility_h, df$mean_volatility_m,
                       df$mean_volatility_0_h, df$mean_volatility_0_m,
                       df$normalised_volatility_h, df$normalised_volatility_m,
                       df$normalised_volatility_0_h, df$normalised_volatility_0_m))
df2 = df2[2:nrow(df2),]
colnames(df2) = c("Accession", "Hamming", "Miyata", "Hamming_stop0",
                  "Miyata_stop0", "Hamming_norm", "Miyata_norm", 
                  "Hamming_stop0_norm", "Miyata_stop0_norm")

write.csv(df2,"output_data/sequence_volatility_sp_normalised.csv")
