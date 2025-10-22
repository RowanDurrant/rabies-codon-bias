##Codon volatility

## v(c) = (i =1:9) Σ d(acid(c.i), acid(c))
#d = amino acid distance measure: Hamming or Miyata
#Hamming much simpler- probably just use this for now
library(readxl)

`%!in%` = Negate(`%in%`)
nucleotides = c("A", "T", "C", "G")
stop_codons = c("TGA", "TAA", "TAG")
df = as.data.frame(read_excel("output_data/Codon_usage_N.xlsx"))

miyata_distances = read.csv("output_data/archive/miyata_distances.csv", header = T, row.names = 1)

codons = colnames(df)[2:ncol(df)]
hamming = c()
miyata = c()
totals = c()

for(i in codons){
  amino_acid = df[1, i] #identify amino acid codon encodes
  cod_vol_h = 0
  cod_vol_m = 0
  total = 0
  for(j in 1:3){ #at each position:
    nuc = strsplit(i, "")[[1]][j] #identify nucleotide
    new_nucs = nucleotides[nucleotides %!in% nuc]  #possible mutations
    
    for(k in new_nucs){ #for each possible mutation:
      #substitute in mutation
      codon_new = i
      substr(codon_new,j,j) <- k
      if(codon_new %!in% stop_codons){
        total = total+1
        #find new amino acid
        amino_acid_new = df[1, codon_new]
        total = total +1
        if(amino_acid != amino_acid_new){
          cod_vol_h = cod_vol_h + 1
          cod_vol_m = cod_vol_m + miyata_distances[amino_acid, amino_acid_new]
          print(paste(amino_acid, "->", amino_acid_new))
        }
      }
    }
  }
  hamming = append(hamming, cod_vol_h)
  totals = append(totals, total)
  miyata = append(miyata, cod_vol_m)
}

hamming = hamming/totals
miyata = miyata/(totals*max(miyata_distances))

volatility_df = data.frame(cbind(codons, hamming, miyata))
write.csv(volatility_df, "output_data/archive/codon_volatility.csv")
