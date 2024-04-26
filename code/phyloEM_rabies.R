#works in R.4.2.1
library(PhylogeneticEM)
library(phytools)
library("readxl")
#https://pbastide.github.io/PhylogeneticEM/articles/tutorial.html

mytree <- read.tree('sequence_data/trees/all_seqs_noout.nwk')
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
             K_max = 20,                         ## Maximal number of shifts
             parallel_alpha = TRUE,              
             Ncores = 5)
  
myres
  
tips = ultratree$tip.label
colours = c()
metadata = read.csv("sequence_data/metadata.csv")

for(i in tips){
  if(metadata$Clade[metadata$Accession == i] == "Asian SEA2a"){colours = append(colours, "#490092")}
  if(metadata$Clade[metadata$Accession == i] == "Cosmo AF1b"){colours = append(colours, "#b66dff")}
  if(metadata$Clade[metadata$Accession == i] == "Arctic A"){colours = append(colours, "red")}
  if(metadata$Clade[metadata$Accession == i] == "Asian SEA2b"){colours = append(colours, "#006ddb")}
  if(metadata$Clade[metadata$Accession == i] == "Cosmo AM2a"){colours = append(colours, "#6db6ff")}
  if(metadata$Clade[metadata$Accession == i] == "Bat DR"){colours = append(colours, "#ffb6db")}
  if(metadata$Clade[metadata$Accession == i] == "Bat TB1"){colours = append(colours, "#ff6db6")}
  if(metadata$Clade[metadata$Accession == i] == "Bat EF-E2"){colours = append(colours, "#004949")}
  if(metadata$Clade[metadata$Accession == i] == "RAC-SK SCSK"){colours = append(colours, "#b6dbff")}
  if(metadata$Clade[metadata$Accession == i] == "Bat LC"){colours = append(colours, "#009292")}
}

plot(myres)
plot(myres, automatic_colors = F, color_characters = colours)
plot_criterion(myres)
