df_loadings = read.csv("output_data/loadings.csv")
df_comp = read_xlsx("output_data/Nucleotide_comp_codons.xlsx")

df = cbind(df_loadings, df_comp)

pc = c()
comp = c()
cc = c()
for(i in c("PC1", "PC2", "PC3")){
  for(j in c(colnames(df)[63:103])){
    pc = append(pc, i)
    comp = append(comp,j)
   cc = append(cc, cor(df[,j]/100, df[,i]))
  }
}

r_df = as.data.frame(cbind(pc,comp,cc))
r_df$cc_abs = abs(as.numeric(r_df$cc))

#PC1 ~ %G3 r = -0.163281895540991
#PC2 ~ %G3+A3 r = 0.170171196
#PC3 ~ %G3+T3 r = 0.421705393

m1 = glm(formula = `%G3+T3`/100~PC3,
    family = binomial, 
    data = df)
summary(m1)
newdata1 = data.frame(PC3 = seq(-0.25,0.25,by=0.01))
newdata1$`%G3+T3`= predict(m1, newdata = newdata1, type = 'response')
p4 = ggplot(data = df, aes(y = `%G3+T3`/100, x =  PC3))+
  geom_point(colour = "blue", size = 4, alpha = 0.4, stroke = 0)+
  ylab("GT3 content of codon")+ xlab("PC3 loading")+
  geom_line(data = newdata1, aes(x = PC3, y = `%G3+T3`), col = 'red')+
  theme_bw()

m2 = glm(formula = `%G3`/100~PC1,
         family = binomial, 
         data = df)
summary(m2)

m3 = glm(formula = `%G3+A3`/100~PC2,
         family = binomial, 
         data = df)
summary(m3)
  
