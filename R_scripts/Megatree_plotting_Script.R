### A script that plots the COII megatree for appendix S1

### Clear the workspace
rm(list=ls())

### load libraries
library(ape)
library(phytools)

#megatree<-read.nexus("Data/trees/Supplementary_Trees/Supplementary_COII.nex")
megatree<-read.newick("RAxML/Outputs/RAxML_bestTree.Supplementary_COII")
megatree<-root(megatree, outgroup = "Outgroup")
  
name<-read.csv("Metadata/Supplementary_names.csv", stringsAsFactors = FALSE)

megatree$tip.colour<- name$colour[match(megatree$tip.label,name$Sample)]
megatree$tip.label <- name$name[match(megatree$tip.label,name$Sample)]

megatree<-ladderize(megatree)


pdf(file = "~/Desktop/Random Cratopus stuff to dump after review/RAxML for SI tree/fig_S1.1_v3.pdf",width = 25, height = 110)
plot(megatree, tip.color = megatree$tip.colour, cex = 0.8)
dev.off()