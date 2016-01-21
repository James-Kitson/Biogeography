########################################################################################################
######################## Script for plotting Cratopus trees ############################
########################################################################################################

### Clear the workspace
rm(list=ls())

### open APE
library(ape)
library(phytools)
library(plyr)

### read in the trees
my.tree<-read.nexus("Data/resolved_tree.nex")

### check this has worked
is.binary.tree(my.tree)

########################################################################################################
######################## Plot the tree with node ages ############################
########################################################################################################

### Calculate the branch depths for each node (distance from tip) and the node number
node.depths<-as.numeric(branching.times(my.tree))
node.depths<-round(node.depths,7)
nodes<-seq(1, length(my.tree$node.label))
node.depths<-cbind(nodes,node.depths)
node.depths<-as.data.frame(node.depths)

### read in the node ages and heights from the vstat file
HPD<-read.csv("bipartition_ages.csv")

### add the node ages and credible intervals to the node depths
node.depths$Median_age <- round((HPD$Median[match(node.depths$node.depths,HPD$Median_height)]),2)
node.depths$Median_high <- round((HPD$CredInt_Upper[match(node.depths$node.depths,HPD$Median_height)]),2)
node.depths$Median_low <- round((HPD$CredInt_Lower[match(node.depths$node.depths,HPD$Median_height)]),2)

### identify the numbers of flight loss nodes
fl.loss<-c(57,50,39, ### FIX this

           #7,29)

### create a df of just the non-flight loss nodes
fl<-node.depths[-fl.loss,]

### set the number of bootstrap replicates
iter<-1000


########################################################################################################
#### testing age of flight loss nodes vs normal nodes
########################################################################################################

### create a vector to fill with bootstrap values
mean.fl<-numeric(iter)

### loop across the non-flightloss node ages to get a vector of mean node ages where the sample size is the same as the flight loss nodes
for(i in 1:iter){
mean.fl[i]<-mean(sample(fl$Median_age, size =length(fl.loss),
                                      replace=TRUE))
}

### plot a histogram and label the mean flightloss node age
node.samp<-hist(mean.fl,breaks=100,
     xlab="Mean nodal age",
     xlim=c(0,5))
abline(v=mean(node.depths[fl.loss,"Median_age"]),
       col="red")

### test the bootstrapped non-flightloss node ages against the mean age of flight loss nodes
t.test(mean.fl, mu=mean(node.depths[fl.loss,"Median_age"]))

########################################################################################################
#### testing age of Mauritius colonisation nodes vs Reunion colonisation nodes
########################################################################################################

### identify the numbers of Mauritius colonisation nodes
mau.col<-c(39,52,28)
### identify the numbers of Reunion colonisation nodes
reu.col<-c(41,37,34,17,16,15,10,7)

### create a df of just the reunion colonistion nodes
reu.col.age<-node.depths[reu.col,]

### create a vector to fill with bootstrap values
mean.reu.col<-numeric(iter)

### loop across the Reunion colonisation node ages to get a vector of mean node ages where the sample size is the same as the Mauritius nodes
for(i in 1:iter){
  mean.reu.col[i]<-mean(sample(reu.col.age$Median_age, size =length(mau.col),
                               replace=TRUE))
}

### plot a histogram and label the mean flightloss node age
node.samp<-hist(mean.reu.col,breaks=100,
                xlab="Mean nodal age",
                xlim=c(0,5))
abline(v=mean(node.depths[mau.col,"Median_age"]),
       col="red")

### test the bootstrapped non-flightloss node ages against the mean age of flight loss nodes
t.test(mean.reu.col, mu=mean(node.depths[mau.col,"Median_age"]))

########################################################################################################
#### testing in situ node ages vs colonisation ages for Reunion
########################################################################################################

### identify in situ nodes for Reunion
reu.insitu<-c(51,38,31,30,29,18)

### create a vector to fill with bootstrap values
mean.reu.col2<-numeric(iter)

### loop across the colonisation node ages to get a vector of mean node ages where the sample size is the same as the in situ nodes
### done this way as there are more colonisation than insitu nodes
for(i in 1:iter){
  mean.reu.col2[i]<-mean(sample(reu.col.age$Median_age, size =length(reu.insitu),
                               replace=TRUE))
}

### plot a histogram and label the mean flightloss node age
node.samp<-hist(mean.reu.col,breaks=20,
                xlab="Mean nodal age",
                xlim=c(0,5))
abline(v=mean(node.depths[reu.insitu,"Median_age"]),
       col="red")

### test the bootstrapped non-flightloss node ages against the mean age of flight loss nodes
t.test(mean.reu.col, mu=mean(node.depths[reu.insitu,"Median_age"]))

########################################################################################################
#### testing in situ node ages vs colonisation ages for Mauritius
########################################################################################################

### identify in situ nodes for Mauritius
mau.insitu<-c(36,35,33,32,4,3,2,5,6,8,9,11)

### create a df of just the reunion colonistion nodes
mau.insitu.age<-node.depths[mau.insitu,]

### create a vector to fill with bootstrap values
mean.mau.insitu<-numeric(iter)

### loop across the colonisation node ages to get a vector of mean node ages where the sample size is the same as the in situ nodes
for(i in 1:iter){
  mean.mau.insitu[i]<-mean(sample(mau.insitu.age$Median_age, size =length(mau.col),
                               replace=TRUE))
}

### plot a histogram and label the mean flightloss node age
node.samp<-hist(mean.mau.insitu,breaks=100,
                xlab="Mean nodal age",
                xlim=c(0,5))
abline(v=mean(node.depths[mau.col,"Median_age"]),
       col="red")

### test the bootstrapped non-flightloss node ages against the mean age of flight loss nodes
t.test(mean.reu.col, mu=mean(node.depths[mau.col,"Median_age"]))
