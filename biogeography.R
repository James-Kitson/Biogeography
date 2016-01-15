########################################################################################################
######################## Script for plotting and biogeographic analysis of Cratopus trees ############################
########################################################################################################

### Clear the workspace
rm(list=ls())

### open APE
library(ape)
library(phytools)
library(plyr)
library(RColorBrewer)
library(BioGeoBEARS)
library(optimx)
library(FD)
library(snow)
library(parallel)

### read in the tree
my.trees<-read.nexus("Data/All_dating_mcorrected.nex.con.tre")
my.tree<-my.trees[[1]]

test