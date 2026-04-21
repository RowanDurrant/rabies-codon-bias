library(readxl)
library(ggplot2)

df_loadings = read.csv("output_data/loadings.csv")
df_comp = read_xlsx("output_data/Nucleotide_comp_codons.xlsx")
df_cpg = read.csv("output_data/N_CpG.csv")

df = cbind(df_loadings, df_comp)
df = cbind(df, df_cpg[,4:12])

pc = c()
comp = c()
cc = c()
for(i in c("PC1", "PC2", "PC3")){
  for(j in c("cpg","cpg_actual","tpa","tpa_actual","ZAP_optimal_motifs",
             colnames(df)[62:102])){
    pc = append(pc, i)
    comp = append(comp,j)
   cc = append(cc, cor(df[,j]/100, df[,i]))
  }
}

r_df = as.data.frame(cbind(pc,comp,cc))
r_df$cc_abs = abs(as.numeric(r_df$cc))

r_df1 = r_df[r_df$pc == "PC1",]
r_df2 = r_df[r_df$pc == "PC2",]
r_df3 = r_df[r_df$pc == "PC3",]
top2 = rbind(r_df1[rev(order(r_df1$cc_abs)),][1:4,],
             r_df2[rev(order(r_df2$cc_abs)),][1:4,],
             r_df3[rev(order(r_df3$cc_abs)),][1:4,])

top2$upper = NA
top2$lower = NA
top2$label = NA

library(confintr)
for(i in 1:nrow(top2)){
  
  top2$lower[i] =ci_cor(y = df[,top2$pc[i]], x = df[,top2$comp[i]])$interval[1]
  top2$upper[i] =ci_cor(y = df[,top2$pc[i]], x = df[,top2$comp[i]])$interval[2]
  top2$label[i] = paste(top2$pc[i], "~", top2$comp[i])
}

top2$label = factor(top2$label, levels = rev(top2$label))

p4 = ggplot(data = top2, aes(y = as.numeric(cc), x = label))+
  geom_point(colour = "blue")+
  geom_errorbar(aes(ymin=as.numeric(lower), ymax=as.numeric(upper)),
                width = 0.1)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_x_discrete(labels = c("PC3~A3", "PC3~G3",
                              "PC3~GT3", "PC3~AC3",
                              "PC2~C3", "PC2~G3",
                              "PC2~GA3", "PC2~CT3",
                              "PC1~C", "PC1~AT1",
                              "PC1~GC1", 'PC1~C1'))+
  theme_bw()+ ylab("Correlation coefficient") + xlab("")+
  coord_flip()
p4  
