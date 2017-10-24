########################################################################################################
######################## Script for plotting Cratopus trees ############################
########################################################################################################

### Clear the workspace
rm(list=ls())

### open APE
library(ape)
library(phytools)
library(plyr)
library(ggtree)

names<-read.csv("Metadata/Multilocus_names_final.csv", stringsAsFactors = FALSE)
clades<-read.csv("Metadata/Multilocus_clades.csv", stringsAsFactors = FALSE)

#############################################################################
#### check MrBayes vs BEAST topology
#############################################################################

BEASTMCC <- read.beast(file = "BEAST/Tree_and_log_files/Relaxed_clock_modeltest2_Good_age/Relaxed_clock_MCC.tree")

### spit out a newick tree to use with BioGeoBEARS for the relaxed BEAST analysis
write.tree(BEASTMCC@phylo, file = "BEAST/Tree_and_log_files/Relaxed_clock_modeltest2_Good_age/Relaxed_clock_MCC_newick.tree")

BEASTMCC@phylo$tip.label<-names$Alt_label[match(BEASTMCC@phylo$tip.label, names$Name)]

#############################################################################
### Plot the BEAST tree with different annotations and node bars 
#############################################################################

# Plot with posterior support
BEASTMCC_plot_support<-ggtree(BEASTMCC, ndigits=2) +
  #geom_range(range='height_0.95_HPD', color='blue', alpha=.6, size=1) +
  geom_text2(aes(label=posterior, subset=posterior>=0.95), vjust=-0.1,  color='red') +
  geom_text2(aes(label=posterior, subset=posterior<0.95), vjust=-0.1,  color='blue') +
  ggtitle("Support Values") +
  #ggplot2::xlim(0, 25) +
  theme_tree2()
BEASTMCC_plot_support<-BEASTMCC_plot_support %<+% clades + geom_tiplab(aes(colour= colour))

# Plot the node numbers for easy reference
BEASTMCC_plot_nodes<-ggtree(BEASTMCC, ndigits=2) +
  #geom_range(range='height_0.95_HPD', color='blue', alpha=.6, size=1) +
  geom_text2(aes(subset=!isTip, label=node), vjust=-0.1,  color='blue') +
  ggtitle("Node number") +
  #ggplot2::xlim(0, 25) +
  theme_tree2()
BEASTMCC_plot_nodes<-BEASTMCC_plot_nodes %<+% clades + geom_tiplab(aes(colour= colour))

# Plot the node ages
BEASTMCC_plot_nodeage<-ggtree(BEASTMCC, ndigits=2) +
  #geom_range(range='height_0.95_HPD', color='blue', alpha=.6, size=1) +
  geom_text2(aes(label=height, subset=!isTip), vjust=-0.1,  color='red') +
  ggtitle("Age values") +
  #ggplot2::xlim(0, 25) +
  theme_tree2()
BEASTMCC_plot_nodeage<-BEASTMCC_plot_nodeage %<+% clades + geom_tiplab(aes(colour= colour))

# Plot the 95% HPD intervals on node ages
BEASTMCC_plot_nodebars<-ggtree(BEASTMCC, ndigits=2) +
  geom_range(range='height_0.95_HPD', color='blue', alpha=.6, size=1) +
  #geom_text(aes(label=posterior), vjust=-0.1,  color='red') +
  ggtitle("Age range bars") +
  #ggplot2::xlim(0, 25) +
  theme_tree2()
BEASTMCC_plot_nodebars<-BEASTMCC_plot_nodebars %<+% clades + geom_tiplab(aes(colour= colour))

# Write the 95% HPD intervals on node ages
BEASTMCC_plot_nodeageHPD<-ggtree(BEASTMCC, ndigits=2) +
  #geom_range(range='height_0.95_HPD', color='blue', alpha=.6, size=1) +
  geom_text2(aes(label=height_0.95_HPD, subset=!isTip), vjust=-0.1,  color='red') +
  ggtitle("Age range values") +
  #ggplot2::xlim(0, 25) +
  theme_tree2()
BEASTMCC_plot_nodeageHPD<-BEASTMCC_plot_nodeageHPD %<+% clades + geom_tiplab(aes(colour= colour))


pdf(file = "Figures/BEAST_Relaxed_clock_MCC.pdf", 8.27,11.69)
revts(rotate(BEASTMCC_plot_support,59)) + ggplot2::xlim(-15, 20)
revts(rotate(BEASTMCC_plot_nodes,59)) + ggplot2::xlim(-15, 20)
revts(rotate(BEASTMCC_plot_nodeage,59)) + ggplot2::xlim(-15, 20)
revts(rotate(BEASTMCC_plot_nodebars,59)) + ggplot2::xlim(-15, 20)
revts(rotate(BEASTMCC_plot_nodeageHPD,59)) + ggplot2::xlim(-15, 20)
dev.off()

