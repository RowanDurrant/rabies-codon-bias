#works in R.4.2.1
library(PhylogeneticEM)
library(phytools)
library("readxl")
#https://pbastide.github.io/PhylogeneticEM/articles/tutorial.html

mytree <- read.tree('sequence_data/trees/rooted_N/all_seqs_noout.nwk')
plot(mytree, show.tip.label = TRUE)


#choose most likely model
#https://doi.org/10.1016/j.ympev.2013.02.008
#chronomodels <- c('relaxed', 'correlated', 'clock')

ultratree =chronos(mytree, lambda=1, model = 'relaxed')
#relaxed: log-Lik = -283.1224 PHIIC = 2663.58 --> BEST
#correlated: log-Lik = -297.1143, PHIIC = 2690.59
#clock: log-Lik = -289.2297, PHIIC = 1276.46 


mydata <- read.csv('output_data/myseqs_RSDUc.tsv', sep='\t', row.names = 1, header= TRUE)

#mydata <- mydata[,c(1,4,7)]
mydata <- mydata[1:(nrow(mydata)-1),] #remove gannoruwa bat lyssavirus
mydata <- t(mydata)


mydatamat = data.matrix(mydata)
  
plot(ultratree, show.tip.label = TRUE)
 # run this to determine alpha values, but >300 might not work on a laptop
alphas <- find_grid_alpha(ultratree)
  
myres <- PhyloEM(phylo = ultratree,
             Y_data = mydatamat,
             process = "scOU",                   ## scalar OU model
             random.root = TRUE,                 ## Root is stationary (true model)
             stationary.root = TRUE,
             alpha = tail(alphas), #get the 5 highest alpha values of the vector
             K_max = 30,                         ## Maximal number of shifts
             parallel_alpha = TRUE,              
             Ncores = 5)
  
myres
  
tips = ultratree$tip.label
colours = c()

df = as.data.frame(read_excel("output_data/Codon_usage_N.xlsx"))
df = df[2:nrow(df),]
for(i in tips){
  if(df$Host[df$CODONS == i] == "Dog (SEA2a)"){colours = append(colours, "#490092")}
  if(df$Host[df$CODONS == i] == "Dog (AF1b)"){colours = append(colours, "#b66dff")}
  if(df$Host[df$CODONS == i] == "Chinese ferret badger"){colours = append(colours, "#006ddb")}
  if(df$Host[df$CODONS == i] == "Mongoose"){colours = append(colours, "#6db6ff")}
  if(df$Host[df$CODONS == i] == "Vampire bat"){colours = append(colours, "#ffb6db")}
  if(df$Host[df$CODONS == i] == "Free-tailed bat"){colours = append(colours, "#ff6db6")}
  if(df$Host[df$CODONS == i] == "Big brown bat"){colours = append(colours, "#004949")}
  if(df$Host[df$CODONS == i] == "Skunk"){colours = append(colours, "#b6dbff")}
  if(df$Host[df$CODONS == i] == "Hoary bat"){colours = append(colours, "#009292")}
}

plot(myres)
plot(myres, automatic_colors = F, color_characters = colours)
plot_criterion(myres)
