##Codon volatility

## v(c) = (i =1:9) Σ d(acid(c.i), acid(c))
#d = amino acid distance measure: Hamming or Miyata
#Hamming much simpler- probably just use this for now
library(readxl)

`%!in%` = Negate(`%in%`)
nucleotides = c("A", "T", "C", "G")
stop_codons = c("TGA", "TAA", "TAG")
df = as.data.frame(read_excel("output_data/Codon_usage_N.xlsx"))

miyata_distances = read.csv("output_data/miyata_distances.csv", header = T, row.names = 1)

codons = colnames(df)[2:ncol(df)]
hamming = c()
miyata = c()

for(i in codons){
  amino_acid = df[1, i] #identify amino acid codon encodes
  cod_vol_h = 0
  cod_vol_m = 0
  
  for(j in 1:3){ #at each position:
    nuc = strsplit(i, "")[[1]][j] #identify nucleotide
    new_nucs = nucleotides[nucleotides %!in% nuc]  #possible mutations
    
    for(k in new_nucs){ #for each possible mutation:
      #substitute in mutation
      codon_new = i
      substr(codon_new,j,j) <- k
      
      if(codon_new %!in% stop_codons){
        
        #find new amino acid
        amino_acid_new = df[1, codon_new]
        
        if(amino_acid != amino_acid_new){
          cod_vol_h = cod_vol_h + 1
          cod_vol_m = cod_vol_m + miyata_distances[amino_acid, amino_acid_new]
          print(paste(amino_acid, "->", amino_acid_new))
        }
      }
    }
  }
  hamming = append(hamming, cod_vol_h)
  miyata = append(miyata, cod_vol_m)
}

volatility_df = data.frame(cbind(codons, hamming, miyata))
write.csv(volatility_df, "output_data/codon_volatility.csv")

df$mean_volatility_h = NA
df$mean_volatility_m = NA

for(a in 2:nrow(df)){
  total_volatility_h = 0
  total_volatility_m = 0
  codon_count = 0
  for(b in codons){
    total_volatility_h = total_volatility_h + (as.numeric(df[a, b])*as.numeric(volatility_df$hamming[volatility_df$codons == b]))
    total_volatility_m = total_volatility_m + (as.numeric(df[a, b])*as.numeric(volatility_df$miyata[volatility_df$codons == b]))
    codon_count = codon_count + as.numeric(df[a, b])
  }
  df$mean_volatility_h[a] = total_volatility_h/codon_count
  df$mean_volatility_m[a] = total_volatility_m/codon_count
}

df2 = data.frame(cbind(df$CODONS, df$mean_volatility_h, df$mean_volatility_m))
df2 = df2[2:nrow(df2),]
colnames(df2) = c("Accession", "Hamming", "Miyata")
write.csv(df2,"output_data/sequence_volatility.csv")
