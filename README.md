# Rabies codon bias

Repo for the manuscript "Differences in codon usage and CpG content between host species-specific rabies clades". 

## Data

All genetic data and metadata are stored in the sequence_data folder. Output from CAICal, DinuQ and my scripts are in the output_data folder.

## Code

This code has been written and run on R version 4.3.2. Each script loosely corresponds to one figure:

- Figure 1 is generated using IQTree and phylogeny.R.
- Figure 2 uses ENC.R.
- Figure 3 uses RSCU.R.
- Figure 4 uses both ENC-GC3.R and Neutrality plot.R.
- Figure 5 uses PCA.R.
- Figure 6 uses data from CoCoPUTS analysed with CAICal and plotted in CAI.R.
- Figure 7 uses CpG.R.
- Figure 8 is generated using the ancestral sequence reconstruction capability of IQtree and plotted in ancestral_state.R.
- The RSDUc values for Figure 9 are calculated using DinuQ in SDUc.ipynb and phylogeneticEM analysis is caried out in phyloEM_rabies.R.
- Figure 10 also used CpG.R.

