##Codon volatility

## v(c) = (i =1:9) Σ d(acid(c.i), acid(c))
#d = amino acid distance measure: Hamming or Miyata
#Hamming much simpler- probably just use this for now
library(readxl)

`%!in%` = Negate(`%in%`)
nucleotides = c("A", "T", "C", "G")
stop_codons = c("TGA", "TAA", "TAG")
df = as.data.frame(read_excel("output_data/Codon_usage_N.xlsx"))
volatility_df = read.csv("output_data/archive/codon_volatility_stop_punished.csv")
codons = colnames(df)[2:ncol(df)]

df$mean_volatility_h = NA

for(a in 2:nrow(df)){
  total_volatility_h = 0
  
  codon_count = 0
  for(b in codons){
    total_volatility_h = total_volatility_h + (as.numeric(df[a, b])*as.numeric(volatility_df$hamming[volatility_df$codons == b]))

    codon_count = codon_count + as.numeric(df[a, b])
  }
  
  df$mean_volatility_h[a] = total_volatility_h/codon_count

}  

df2 = data.frame(cbind(df$CODONS, 
                       df$mean_volatility_h, df$normalised_volatility_h))
df2 = df2[2:nrow(df2),]
colnames(df2) = c("Accession", "Hamming")

write.csv(df2,"output_data/archive/sequence_volatility_sp_normalised.csv")
