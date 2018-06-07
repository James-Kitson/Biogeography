# Kitson et al Cratopine Biogeography - R Scripts

The R scripts in this folder are as follows:

- [Plot_MCC_ggtree.R](https://github.com/James-Kitson/Biogeography/blob/master/R_scripts/Plot_MCC_ggtree.R): This script plots the BEAST MCC tree as presented in the supplementary material.
- [BEAST_array_unstratified_DECJX.R](https://github.com/James-Kitson/Biogeography/blob/master/R_scripts/BEAST_array_unstratified_DECJX.R): This script performs the DEC+J+X analysis on a cluster when submitted by the submit_BEAST_DECJX_array.sh script. On a slurm job manager this is typically submitted as "sbatch --array=1-1000 submit_BEAST_DECJX_array.sh"
- [sample_BEAST_trees.R](https://github.com/James-Kitson/Biogeography/blob/master/R_scripts/sample_BEAST_trees.R): This script takes logcombiner combined tree file sand draws out 1000 random trees.
- [Node_age_stats.R](https://github.com/James-Kitson/Biogeography/blob/master/R_scripts/Node_age_stats.R): This script takes the output from the BioGeoBEARs array run performed on a cluster, combines the output and performs the statistical tests required.
- [Megatree_plotting_Script.R](https://github.com/James-Kitson/Biogeography/blob/master/R_scripts/Megatree_plotting_Script.R): This script takes the output from RAxML and produces the species delimitation tree for the supplementary information.

### Note
Some of these analyses will need to be performed on a computing cluster as it will simply take far too long to perform the required BioGeoBEARs analyses sequentially on 1000 trees. As a result there is quite a bit of uploading and downloading between computers that I can't predict for all users. So be aware that you may need to alter file paths and data structures appropriately to get this to work on your machine.