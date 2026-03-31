df_loadings = read.csv("output_data/loadings.csv")
df_comp = read_xlsx("output_data/Nucleotide_comp_codons.xlsx")

df = cbind(df_loadings, df_comp)

pc = c()
comp = c()
r2 = c()
for(i in c("PC1", "PC2", "PC3")){
  for(j in c(colnames(df)[63:103])){
    pc = append(pc, i)
    comp = append(comp,j)
    r2 = append(r2, 
                summary(glm(formula = df[,j]/100~df[,i],
                            family = binomial))$coefficients[2,4])
  }
}

r_df = as.data.frame(cbind(pc,comp,r2))

plot(data = df, `%G3+T3` ~ PC3)
