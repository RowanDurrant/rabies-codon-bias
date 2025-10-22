codon_mean_calc = function(codon_vect){
  mean = sum(as.numeric(codon_vect))/length(codon_vect)
  freq = codon_vect - mean
  return(freq)
}

codon_mean = function(codon_freqs){
  codon_freqs[1:2] = codon_mean_calc(codon_freqs[1:2])
  codon_freqs[3:8] = codon_mean_calc(codon_freqs[3:8])
  codon_freqs[9:11] = codon_mean_calc(codon_freqs[9:11])
  codon_freqs[12:15] = codon_mean_calc(codon_freqs[12:15])
  codon_freqs[16:21] = codon_mean_calc(codon_freqs[16:21])
  codon_freqs[22:25] = codon_mean_calc(codon_freqs[22:25])
  codon_freqs[26:29] = codon_mean_calc(codon_freqs[26:29])
  codon_freqs[30:33] = codon_mean_calc(codon_freqs[30:33])
  codon_freqs[34:35] = codon_mean_calc(codon_freqs[34:35])
  codon_freqs[36:37] = codon_mean_calc(codon_freqs[36:37])
  codon_freqs[38:39] = codon_mean_calc(codon_freqs[38:39])
  codon_freqs[40:41] = codon_mean_calc(codon_freqs[40:41])
  codon_freqs[42:43] = codon_mean_calc(codon_freqs[42:43])
  codon_freqs[44:45] = codon_mean_calc(codon_freqs[44:45])
  codon_freqs[46:47] = codon_mean_calc(codon_freqs[46:47])
  codon_freqs[48:49] = codon_mean_calc(codon_freqs[48:49])
  codon_freqs[50:55] = codon_mean_calc(codon_freqs[50:55])
  codon_freqs[56:59] = codon_mean_calc(codon_freqs[56:59])
  
  return(codon_freqs)
}