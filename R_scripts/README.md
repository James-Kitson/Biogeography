# Kitson et al Cratopine Biogeography

The R scripts in this folder are as follows:

- [Plot_MCC_ggtree.R](https://github.com/James-Kitson/Biogeography/blob/master/R_scripts/Plot_MCC_ggtree.R): This script plots the BEAST MCC tree as presented in the supplementary material.
- [Supplementary_trees_R](https://github.com/James-Kitson/Biogeography/blob/master/R_scripts/Supplementary_trees.R): This script plots the COII tree for the definition of OTUs and the seperate locus trees.
- [BEAST_array_unstratified_DECJX.R](https://github.com/James-Kitson/Biogeography/blob/master/R_scripts/BEAST_array_unstratified_DECJX.R): This script performs the DEC+J+X analysis on a cluster when submitted by the submit_BEAST_DECJX_array.sh script. On a slurm job manager this is typically submitted as "sbatch --array=1-1000 submit_BEAST_DECJX_array.sh"

Each analysis folder contains:
- Input files
- Metadata files
- Output consensus or MCC tree files (users wishing to obtain full tree samples should clone the repository and rerun the analyses)

***********************************************************************************************************

The Rnotebook should be executed in the root of the cloned repository. All filepaths are correctly specified for the various analyses.