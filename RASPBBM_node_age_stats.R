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
my.trees<-read.nexus("Data/All_dating_mcorrected.nex.con.tre")

### process the trees
my.tree<-my.trees[[1]] ### raw outpout from MrBayes with polytomies
#my.tree<-multi2di(my.trees[[1]], random=TRUE) ### MrBayes output with polytomies resolved for ML character reconstruction

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
HPD<-read.csv("Data/bipartition_ages.csv")

### add the node ages and credible intervals to the node depths
node.depths$Median_age <- round((HPD$Median[match(node.depths$node.depths,HPD$Median_height)]),2)
node.depths$Median_high <- round((HPD$CredInt_Upper[match(node.depths$node.depths,HPD$Median_height)]),2)
node.depths$Median_low <- round((HPD$CredInt_Lower[match(node.depths$node.depths,HPD$Median_height)]),2)

### identify the numbers of flight loss nodes
fl.loss<-c(52,45,34,7,29)

### create a df of just the non-flight loss nodes
fl<-node.depths[-fl.loss,]

### set the number of bootstrap replicates
iter<-1000


########################################################################################################
#### testing age of flight loss nodes vs normal nodes
########################################################################################################

### create a vector to fill with bootstrap values
mean.fl<-numeric(iter)

### repeatedly sample the non-flightloss node ages to get a vector of mean node ages where the sample size is the same as the flight loss nodes
for(i in 1:iter){
  mean.fl[i]<-mean(sample(fl$Median_age, size =length(fl.loss),
                          replace=TRUE))
}

### plot a histogram and label the mean flightloss node age
png(filename="Diagrams/flightloss_vs_flight_BBM_biogeography", width=1000)
node.samp<-hist(mean.fl,breaks=100,
                main="Histogram of flightloss vs normal node ages ML bioeography",
                xlab="Mean nodal age",
                xlim=c(0,4))
abline(v=mean(node.depths[fl.loss,"Median_age"]),
       col="red")
title(sub="flight loss nodes are significantly younger, t = 41.4841, df = 999, p-value < 2.2e-16")
dev.off()
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

### plot a histogram and label the mean mauritian colonisation age
png(filename="Diagrams/Maucol_vs_Reucol_BBM_biogeography", width=1000)
node.samp<-hist(mean.reu.col,breaks=100,
                main="Histogram of mauritian vs reunion colonisation ages",
                xlab="Mean nodal age",
                xlim=c(0,5))
abline(v=mean(node.depths[mau.col,"Median_age"]),
       col="red")
title(sub="Mauritian colonisation nodes are significantly younger than reunion colonisations, t = 54.0667, df = 999, p-value < 2.2e-16")
dev.off()

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

### plot a histogram and label the mean insitu speciation age
png(filename="Diagrams/insituReu_vs_Reucol_BBM_biogeography", width=1000)
node.samp<-hist(mean.reu.col,breaks=20,
                main="Histogram of reunion insitu speciation ages",
                xlab="Mean nodal age",
                xlim=c(0,5))
abline(v=mean(node.depths[reu.insitu,"Median_age"]),
       col="red")
title(sub="reunion insitu speciation nodes are significantly younger than colonisation nodes, t = 3.2722, df = 999, p-value = 0.001104")
dev.off()

### test the bootstrapped non-flightloss node ages against the mean age of flight loss nodes
t.test(mean.reu.col, mu=mean(node.depths[reu.insitu,"Median_age"]))

########################################################################################################
#### testing in situ node ages vs colonisation ages for Mauritius - This needs to be in line with previous graph
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

### plot a histogram and label the mean insitu speciation node age
png(filename="Diagrams/insituMau_vs_Maucol_BBM_biogeography", width=1000)
node.samp<-hist(mean.mau.insitu,breaks=100,
                main="Histogram of mauritian insitu speciation ages",
                xlab="Mean nodal age",
                xlim=c(0,5))
abline(v=mean(node.depths[mau.col,"Median_age"]),
       col="red")
title(sub="Mauritian insitu speciation nodes are significantly older than colonisation nodes, t = -4.0838, df = 999, p-value = 4.786e-05")
dev.off()

### test the bootstrapped mauritian colonisation node ages against the mean age of in situ speciation nodes
t.test(mean.mau.insitu, mu=mean(node.depths[mau.col,"Median_age"]))