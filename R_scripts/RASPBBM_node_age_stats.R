########################################################################################################
######################## Script for performng stats on Cratopus trees ############################
########################################################################################################

## @knitr processdata1

### Clear the workspace
rm(list=ls())

### open APE
library(ape)
library(phytools)
library(plyr)

### read in the trees
my.trees<-read.nexus("Data/trees/concatenated_loci/All_dating_mcorrected.nex.con.tre")

### process the trees
my.tree<-my.trees[[1]] ### raw outpout from MrBayes with polytomies

########################################################################################################

### Calculate the branch depths for each node (distance from tip) and the node number
node.depths<-as.numeric(branching.times(my.tree))
node.depths<-round(node.depths,7)
nodes<-seq(length(my.tree$tip.label)+1, length(my.tree$tip.label)+length(my.tree$node.label),1)
node.depths<-cbind(nodes,node.depths)
node.depths<-as.data.frame(node.depths)

### read in the node ages and heights from the vstat file
HPD<-read.csv("Data/bipartition_ages.csv")

### add the node ages and credible intervals to the node depths
node.depths$Median_age <- round((HPD$Median[match(node.depths$node.depths,HPD$Median_height)]),2)
node.depths$Median_high <- round((HPD$CredInt_Upper[match(node.depths$node.depths,HPD$Median_height)]),2)
node.depths$Median_low <- round((HPD$CredInt_Lower[match(node.depths$node.depths,HPD$Median_height)]),2)

########################################################################################################

### define node numbers for the groups of nodes we will test.
### Nodes with Flight loss
fl.loss.nodes<-c(65,87,92,110)
### Colonisations of Réunion from Mauritius
Reu.col.nodes<-c(65,68,73,74,75,92,95,99)
### internal speciation events on Réunion
Reu.insitu.nodes<-c(76,77,78,87,88,89,96,109)
### Colonisations of Mauritius from Réunion
Mau.col.nodes<-c(86,97,110)
### internal speciation event on Mauritius
Mau.insitu.nodes<-c(60,61,62,63,64,66,67,69,72,76,87,88,89,90,91,93,94,98)

########################################################################################################

### set the number of bootstrap replicates
iter<-1000

########################################################################################################
#### testing age of flight loss nodes vs normal nodes
########################################################################################################

### create a df of just the non-flight loss nodes by subsetting node.depths to everything except the flight loss nodes.
fl<-node.depths[!node.depths$nodes %in% fl.loss.nodes, ]

### create a vector to fill with bootstrap values
mean.fl<-numeric(iter)

## @knitr BBMflightloss

### repeatedly sample the non-flightloss node ages to get a vector of mean node ages where the sample size is the same as the flight loss nodes
for(i in 1:iter){
  mean.fl[i]<-mean(sample(fl$Median_age, size =length(fl.loss.nodes),
                          replace=FALSE))
}

### plot a histogram and label the mean flightloss node age
node.samp<-hist(mean.fl,breaks=100,
                main="Histogram of flightloss vs normal node ages",
                xlab="Mean age of resampled flight capable nodes",
                xlim=c(0,4))
abline(v=mean(node.depths$Median_age[node.depths$nodes %in% fl.loss.nodes]),
       col="red")


## @knitr BBMflightlossttest
### test the bootstrapped non-flightloss node ages against the mean age of flight loss nodes
t.test(mean.fl, mu=mean(node.depths$Median_age[node.depths$nodes %in% fl.loss.nodes]))
BBMfl.test<-t.test(mean.fl, mu=mean(node.depths$Median_age[node.depths$nodes %in% fl.loss.nodes]))

########################################################################################################
#### testing age of Mauritius colonisation nodes vs Réunion colonisation nodes
########################################################################################################

## @knitr processdata2

### create a df of just the Réunion colonistion nodes
Reu.col<-node.depths[node.depths$nodes %in% Reu.col.nodes, ]

### create a vector to fill with bootstrap values
mean.Reu.col<-numeric(iter)

### loop across the Réunion colonisation node ages to get a vector of mean node ages where the sample size is the same as the Mauritius nodes
for(i in 1:iter){
  mean.Reu.col[i]<-mean(sample(Reu.col$Median_age, size =length(Mau.col.nodes),
                               replace=FALSE))
}

### plot a histogram and label the mean mauritian colonisation age
## @knitr BBMMauvsReu
node.samp<-hist(mean.Reu.col,breaks=100,
                main="Histogram of mauritian vs Réunion colonisation ages",
                xlab="Mean age of resampled Réunion colonisation nodes",
                xlim=c(0,5))
abline(v=mean(node.depths$Median_age[node.depths$nodes %in% Mau.col.nodes]),
       col="red")

### test the bootstrapped non-flightloss node ages against the mean age of flight loss nodes
## @knitr BBMMauvsReuttest
t.test(mean.Reu.col, mu=mean(node.depths$Median_age[node.depths$nodes %in% Mau.col.nodes]))
BBMMauvsReu.test<-t.test(mean.Reu.col, mu=mean(node.depths$Median_age[node.depths$nodes %in% Mau.col.nodes]))

########################################################################################################
#### testing in situ node ages vs colonisation ages for Réunion
########################################################################################################

## @knitr processdata3

### identify in situ nodes for Réunion
Reu.insitu<-node.depths[node.depths$nodes %in% Reu.insitu.nodes, ]

### create a vector to fill with bootstrap values
mean.Reu.col2<-numeric(iter)

### loop across the colonisation node ages to get a vector of mean node ages where the sample size is the same as the in situ nodes
### done this way as there are more colonisation than in situ nodes
for(i in 1:iter){
  mean.Reu.col2[i]<-mean(sample(Reu.col$Median_age, size =length(Reu.insitu.nodes),
                                replace=FALSE))
}

### plot a histogram and label the mean in situ speciation age
## @knitr BBMinsituReu
node.samp<-hist(mean.Reu.col2,breaks=20,
                main="Histogram of Réunion colonisation ages",
                xlab="Mean age of resampled Réunion colonisation nodes",
                xlim=c(0,5))
abline(v=mean(node.depths$Median_age[node.depths$nodes %in% Reu.insitu.nodes]),
       col="red")

### test the bootstrapped non-flightloss node ages against the mean age of flight loss nodes
## @knitr BBMinsituReuttest
### as there are exactly the same number of Réunion colonisations as in situ speciation events, the resampling doesn't work and the t-test below isn't really appropriate
### t.test(mean.Reu.col, mu=mean(node.depths$Median_age[node.depths$nodes %in% Reu.insitu.nodes]))

### a more appropriate t-test is all insitu nodes vs all colonisation nodes.
t.test(Reu.col$Median_age, Reu.insitu$Median_age)
BBMReuinsitu.test<-t.test(Reu.col$Median_age, Reu.insitu$Median_age)

########################################################################################################
#### testing in situ node ages vs colonisation ages for Mauritius
########################################################################################################

## @knitr processdata4

### identify in situ nodes for Mauritius
Mau.insitu<-node.depths[node.depths$nodes %in% Mau.insitu.nodes, ]

### create a vector to fill with bootstrap values
mean.Mau.insitu<-numeric(iter)

### loop across the colonisation node ages to get a vector of mean node ages where the sample size is the same as the in situ nodes
for(i in 1:iter){
  mean.Mau.insitu[i]<-mean(sample(Mau.insitu$Median_age, size =length(Mau.col.nodes),
                                  replace=FALSE))
}

### plot a histogram and label the mean in situ speciation node age
## @knitr BBMinsituMau
node.samp<-hist(mean.Mau.insitu,breaks=100,
                main="Histogram of mauritian in situ speciation ages",
                xlab="Mean age of resampled mauritian in situ speciation nodes",
                xlim=c(0,5))
abline(v=mean(node.depths$Median_age[node.depths$nodes %in% Mau.col.nodes]),
       col="red")

### test the bootstrapped mauritian colonisation node ages against the mean age of in situ speciation nodes
## @knitr BBMinsituMauttest
t.test(mean.Mau.insitu, mu=mean(node.depths$Median_age[node.depths$nodes %in% Mau.col.nodes]))
BBMMauinsitu.test<-t.test(mean.Mau.insitu, mu=mean(node.depths$Median_age[node.depths$nodes %in% Mau.col.nodes]))
