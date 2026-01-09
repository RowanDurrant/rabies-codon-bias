# Rabies codon bias

Repo for the manuscript "Differences in codon usage between host-species-specific rabies virus clades are driven by UpA and purine content". 

## Data

All genetic data and metadata are stored in the sequence_data folder. Output from CAICal, DinuQ and my scripts are in the output_data folder.

## Code

This code has been written and run on R version 4.3.2. Each script loosely corresponds to one figure:

- Figure 1 is generated using IQTree, phylogeny.R and ENC.R.
- Figure 2 uses RSCU.R.
- Figure 3 uses PCA.R and cpg_pca.R.
- Figure 4 uses data from CoCoPUTS analysed with CAICal and plotted in CAI.R and compare_cpg_cai.R.
- Figure 5 uses CpG.R.
- Figure 6 is generated using the ancestral sequence reconstruction capability of IQtree and plotted in ancestral_state.R.
- Figure 7 is plotted in CpG.R and ZAP_CpG_locs.R.

