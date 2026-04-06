library(tidyr)
library(ggplot2)
library(boot)
library(MCMCglmm)
library(ape)

df = read.csv("output_data/HIVE-CUTS_CAI.csv")
cpg = read.csv("output_data/N_CpG.csv")

df = merge(df, cpg, by.x = "Name", by.y = "accessions")

canis_familiaris_eCAI = 0.784
desmodus_rotundus_eCAI = 0.754
eptesicus_fuscus_eCAI = 0.753
vulpes_lagopus_eCAI = 0.774

df$normalised = NA
df$normalised[df$Reference == "Canis familiaris"] = df$CAI[df$Reference == "Canis familiaris"]/canis_familiaris_eCAI
df$normalised[df$Reference == "Eptesicus fuscus"] = df$CAI[df$Reference == "Eptesicus fuscus"]/eptesicus_fuscus_eCAI
df$normalised[df$Reference == "Desmodus rotundus"] = df$CAI[df$Reference == "Desmodus rotundus"]/desmodus_rotundus_eCAI
df$normalised[df$Reference == "Vulpes lagopus"] = df$CAI[df$Reference == "Vulpes lagopus"]/vulpes_lagopus_eCAI

df$clade = factor(df$clade, c("Cosmopolitan AF1b", "Cosmopolitan AM2a", "Arctic A", "Asian SEA2a", 
                              "Asian SEA2b", 
                              "Bats DR","Bats TB1", "Bats LC",
                              "Bats EF-E2","RAC-SK SCSK"))

mypal = c("#332288","#88CCEE","#CCDDAA","#44AA99","#117733",  
          "#DDCC77","#999933", "#AA4499","#CC6677","#882255")
mylabels = c("Cosmo AF1b\n(dog)",
             "Cosmo AM2a\n(mongoose)",
             "Arctic A\n(arctic fox)",
             "Asian SEA2a\n(dog)",
             "Asian SEA2b\n(CFB)",
             "Bat DR\n(vampire bat)",
             "Bat TB1\n(Mexican free\n-tailed bat)",
             "Bat LC\n(hoary bat)",
             "Bat EF-E2\n(big brown bat)",
             "RAC-SK SCSK\n(skunk)")

p0 = ggplot(data = df, aes(x = cpg, y = normalised))+
  geom_point(aes(shape = Reference,
                 colour = clade))+
  geom_smooth(se = F, method = "lm", col = "black")+
  theme_bw()+ 
  scale_colour_manual(values = mypal, labels = mylabels,
                    name = "Clade")+
  scale_shape_manual(name = "Reference host",
                     values = c(15,16,17,3),
                     labels = c(expression(italic("Canis familiaris")),
                                expression(italic("Desmodus rotundus")),
                                expression(italic("Eptesicus fuscus")),
                                expression(italic("Vulpes lagopus"))))+
  xlab("Obs/Exp CpG content")

summary(lm(data = df, normalised ~ cpg))
bootstrap = boot(df,function(data,indices)
  summary(lm(normalised ~ cpg,data[indices,]))$adj.r.squared,R=10000)
quantile(bootstrap$t,c(0.025,0.975))

p1 = ggplot(data = df, aes(x = tpa, y = normalised))+
  geom_point(aes(shape = Reference,
                 colour = clade))+
  geom_smooth(se = F, method = "lm", col = "black")+
  theme_bw()+ 
  scale_colour_manual(values = mypal, labels = mylabels,
                      name = "Clade")+
  scale_shape_manual(name = "Reference host",
                     values = c(15,16,17,3),
                     labels = c(expression(italic("Canis familiaris")),
                                expression(italic("Desmodus rotundus")),
                                expression(italic("Eptesicus fuscus")),
                                expression(italic("Vulpes lagopus"))))+
  xlab("Obs/Exp UpA content") + ylab("nCAI")
summary(lm(data = df, normalised ~ tpa))
bootstrap = boot(df,function(data,indices)
  summary(lm(normalised ~ tpa,data[indices,]))$adj.r.squared,R=10000)
print(quantile(bootstrap$t,c(0.025,0.5,0.975)))

tree = read.tree("sequence_data/tree/outgroup_removed.nwk")
ultra.tree = chronos(tree, 0)
Ainv.phylo<-inverseA(ultra.tree,nodes="TIPS")$Ainv


#I think reference species should be included somewhere? Random effect??
m1.phylo<-MCMCglmm(normalised ~ tpa, random=~Name, 
                   ginverse=list(Name=Ainv.phylo), data=df)
summary(m1.phylo) #<0.001
posterior.mode(m1.phylo$VCV[,1]/rowSums(m1.phylo$VCV)) #0.8184092
HPDinterval(m1.phylo$VCV[,1]/rowSums(m1.phylo$VCV)) #0.7898131 0.8545874

m2.phylo<-MCMCglmm(normalised ~ cpg + Reference, random=~Name, 
                   ginverse=list(Name=Ainv.phylo), data=df)
summary(m2.phylo) #<0.001
posterior.mode(m2.phylo$VCV[,1]/rowSums(m2.phylo$VCV)) #0.8496247
HPDinterval(m2.phylo$VCV[,1]/rowSums(m2.phylo$VCV)) #0.816564 0.8736289
