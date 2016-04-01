########################################################################################################
######################## Script for plotting Cratopus trees ############################
########################################################################################################

## @knitr MLtreeprocess

### Clear the workspace
rm(list=ls())

### open APE
library(ape)
library(phytools)
library(plyr)

### read in the trees
my.tree<-read.nexus("Data/trees/Cratopus_MCC_tree.nex")
my.tree2<-read.nexus("Data/resolved_tree.nex")

### read in the list of names
name<-read.csv("Data/all_names.csv")
### read.csv turns text into factors, this gets messy later when plotting
### so make it character data
name<-data.frame(lapply(name, as.character), stringsAsFactors=FALSE)

### the next line uses match to perform the same function as vlookup in excel
my.tree$tip.label <- (name$Alt_label[match(my.tree$tip.label,name$Name)])
my.tree2$tip.label <- (name$Alt_label[match(my.tree2$tip.label,name$Name)])

### use cophylo to rotate the MCC tree relative to the consensus tree
obj<-cophylo(my.tree,my.tree2)

### plot and look at the tip associations
plot(obj)