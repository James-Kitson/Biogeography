# Kitson et al Cratopine Biogeography

The R scripts in this folder are as follows:

- [Plot_MCC_ggtree.R](https://github.com/James-Kitson/Biogeography/blob/master/R_scripts/Plot_MCC_ggtree.R): This script plots the BEAST MCC tree as presented in the supplementary material.
- [Supplementary_trees.R](https://github.com/James-Kitson/Biogeography/blob/master/R_scripts/Supplementary_trees.R): This script plots the COII tree for the definition of OTUs and the seperate locus trees.
- [BEAST_array_unstratified_DECJX.R](https://github.com/James-Kitson/Biogeography/blob/master/R_scripts/BEAST_array_unstratified_DECJX.R): This script performs the DEC+J+X analysis on a cluster when submitted by the submit_BEAST_DECJX_array.sh script. On a slurm job manager this is typically submitted as "sbatch --array=1-1000 submit_BEAST_DECJX_array.sh"
- [sample_BEAST_trees.R](https://github.com/James-Kitson/Biogeography/blob/master/R_scripts/sample_BEAST_trees.R): This script takes logcombiner combined tree file sand draws out 1000 random trees.
-[Model_testing.R](https://github.com/James-Kitson/Biogeography/blob/master/R_scripts/Model_testing.R): Ascript based heavily on Nick Matzke's own scrips for model testing across biogeographic models in BioGeoBEARS.
