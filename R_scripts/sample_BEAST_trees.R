########################################################################################################
######################## Script for subsetting BEAST Cratopus trees ############################
########################################################################################################

### Clear the workspace
rm(list=ls())

### open APE
library(ape)

### read in the combined post-burnin trees from logcombiner
combined_BEAST<-read.nexus("BEAST/Tree_and_log_files/Relaxed_clock_modeltest2_Good_age/combined_relaxed_trees.trees")

### pull 100 random trees out of the list
sampled_BEAST_trees<-sample(combined_BEAST, 1000, replace = FALSE)

### write these out in newick format for BioGeoBEARs to process
write.tree(sampled_BEAST_trees, file="BEAST/Tree_and_log_files/Relaxed_clock_modeltest2_Good_age/BEAST_newick.trees")
