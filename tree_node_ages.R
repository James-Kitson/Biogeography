  ########################################################################################################
######################## Script for plotting Cratopus trees ############################
########################################################################################################

### Clear the workspace
rm(list=ls())

### set working directory
setwd("~/Documents/Projects/Cratopus project/R analyses/Biogeography")
dir.create(paste(getwd(),"/Diagrams/",sep=""))
out<-paste(getwd(),"/Diagrams/",sep="")

### open APE
library(ape)
library(phytools)
library(plyr)
library(RColorBrewer)

### read in the tree
my.trees<-read.nexus("All_dating_mcorrected.nex.con.tre")
my.tree<-my.trees[[1]]

########################################################################################################
######################## Plot the tree with node ages ############################
########################################################################################################

### Calculate the branch depths for each node (distance from tip) and the node number
node.depths<-as.numeric(branching.times(my.tree))
node.depths<-round(node.depths,7)
nodes<-seq(1, length(my.tree$node.label))
node.depths<-cbind(nodes,node.depths)
node.depths<-as.data.frame(node.depths)


########################################################################################################
####### the next bit is used when manually correcting the ages on branch lengths with a fixed clock ####
########################################################################################################

### read in the node ages and heights from the vstat file
HPD<-read.csv("bipartition_ages.csv")

### add the node ages and credible intervals to the node depths
node.depths$Median_age <- round((HPD$Median[match(node.depths$node.depths,HPD$Median_height)]),2)
node.depths$Median_high <- round((HPD$CredInt_Upper[match(node.depths$node.depths,HPD$Median_height)]),2)
node.depths$Median_low <- round((HPD$CredInt_Lower[match(node.depths$node.depths,HPD$Median_height)]),2)

my.tree$node.label<-paste("Node: ",node.depths$nodes,"\n",node.depths$Median_age,"(high_95: ",node.depths$Median_high," /low_95: ",node.depths$Median_low,")", sep="")

########################################################################################################
########################################################################################################

### read in the list of names
name<-read.csv("all_names.csv")
### read.csv turns text into factors, this gets messy later when plotting
### so make it character data
name<-data.frame(lapply(name, as.character), stringsAsFactors=FALSE)
### the next line uses match to perform the same function as vlookup in excel
my.tree$tip.label <- (name$Alt_label[match(my.tree$tip.label,name$Name)])

#### plot the tree with node ages listed on nodes
pdf(paste(out,"nodeages.pdf"),70,70)
plot(my.tree,
     show.node.label=TRUE)
dev.off()

########################################################################################################
######################## Plot the tree with support values and a nice scale ############################
########################################################################################################

### replace the codes with informative names
tree.rename<-my.trees[[1]]
### the next line uses match to perform the same function as vlookup in excel
tree.rename$tip.label <- (name$Alt_label[match(tree.rename$tip.label,name$Name)])

### COnvert node lables to numeric and round to two dp
tree.rename$node.label<-as.numeric(tree.rename$node.label)
tree.rename$node.label<-round(tree.rename$node.label,digits=2)

### check all the various labels are in the correct format
str(tree.rename$tip.label)
str(tree.rename$node.label)

#########################################################
### Plot the main combined tree with support values
pdf(file=paste(out,"nodevalues.pdf",sep=""), 30, 30)
###svg(file=paste(out,"nodevalues.svg",sep=""), 30, 30)
plot(tree.rename,
     show.node.label=FALSE,
     cex=2,
     x.lim=0.3,
     label.offset=0.001)
nodelabels(tree.rename$node.label,adj=c(1.1,1.3),frame="none",
           col=ifelse(tree.rename$node.label>0.9,"red",
                      ifelse(tree.rename$node.label>=0.5 & tree.rename$node.label<0.9,"blue","#0000ff00")),cex=1)

### make an offset for the axis as R won't draw it from the tip to the root. The offset is a negative starting point for the axis equivalent to the 
### round ing up we do at the root end of the axis i.e. if we round 4.79 Mya to 5 Mya then we need to offset by minus ~0.21Ma of distance measured in
### branch lengths. To do this we divide the root height by the root age and multiply by the difference between the oldest value on the axis and
### the oldest value on the tree.
offset<-(round_any(max(HPD$Median),0.5)-max(HPD$Median))*(max(nodeHeights(tree.rename))/max(HPD$Median))

## put on a the correct axis
axis(side=1,
     cex.axis=1,
     padj=1,
     at=seq(-offset, max(nodeHeights(tree.rename)), by=(max(nodeHeights(tree.rename))+offset)/(round_any(max(HPD$Median),0.5)/0.5)),
     labels=seq(round_any(max(HPD$Median),0.5),0,by=-0.5))

### put on the 125 kya calibration point
abline(v=max(nodeHeights(tree.rename))-((max(nodeHeights(tree.rename))/max(HPD$Median))*0.125), col="green")

dev.off()

########################################################################################################
######################## Plot the tree with node numbers ############################
########################################################################################################

### replace the codes with informative names
tree.rename2<-my.trees[[1]]
### the next line uses match to perform the same function as vlookup in excel
tree.rename2$tip.label <- (name$Alt_label[match(tree.rename2$tip.label,name$Name)])

### COnvert node lables to numeric and round to two dp
tree.rename2$node.label<-as.numeric(tree.rename2$node.label)
tree.rename2$node.label<-round(tree.rename2$node.label,digits=2)

### check all the various labels are in the correct format
str(tree.rename2$tip.label)
str(tree.rename2$node.label)

### make the node labels into node numbers
tree.rename2$node.label<-seq(1, length(tree.rename2$node.label),1)

#########################################################
### Plot the main combined tree with node numbers
pdf(file=paste(out,"nodenumbers.pdf",sep=""), 30, 30)
###svg(file=paste(out,"nodevalues.svg",sep=""), 30, 30)
plot(tree.rename2,
     show.node.label=FALSE,
     cex=2,
     x.lim=0.3,
     label.offset=0.001)
nodelabels(tree.rename2$node.label,adj=c(1.1,1.3),cex=1)

### make an offset for the axis as R won't draw it from the tip to the root. The offset is a negative starting point for the axis equivalent to the 
### round ing up we do at the root end of the axis i.e. if we round 4.79 Mya to 5 Mya then we need to offset by minus ~0.21Ma of distance measured in
### branch lengths. To do this we divide the root height by the root age and multiply by the difference between the oldest value on the axis and
### the oldest value on the tree.
offset<-(round_any(max(HPD$Median),0.5)-max(HPD$Median))*(max(nodeHeights(tree.rename2))/max(HPD$Median))

## put on a the correct axis
axis(side=1,
     cex.axis=1,
     padj=1,
     at=seq(-offset, max(nodeHeights(tree.rename2)), by=(max(nodeHeights(tree.rename2))+offset)/(round_any(max(HPD$Median),0.5)/0.5)),
     labels=seq(round_any(max(HPD$Median),0.5),0,by=-0.5))

dev.off()


########################################################################################################
######################## Do some rough phylogeography ############################
########################################################################################################

### get the tree and resolve polytomies
tree.rename3<-my.trees[[1]]
### resolve zero length branches
###tree.rename3$edge.length<-ifelse(tree.rename3$edge.length==0,0.0000000001,tree.rename3$edge.length)
### check it has worked
###tree.rename3$edge.length

### read in teh island data
islands<-read.csv("dist_ultra.csv")

### calculate the ML character reconstruction
#Cratopus_islands<-ace(islands$island,tree.rename3,type="discrete")
  
Cratopus_anc<-read.csv("RASP_anc.csv", row.names=2, header=TRUE)
Cratopus_anc<-Cratopus_anc[,2:14]

### replace all the islands that are less than 5% probability with zero
Cratopus_anc[,1:12][Cratopus_anc[,1:12] < 5] <- 0
Cratopus_anc$combined<-100-rowSums(Cratopus_anc[,1:12])

Cratopus_anc<-as.matrix(Cratopus_anc)
test<-rowSums(Cratopus_anc)
test

### COnvert node lables to numeric and round to two dp
tree.rename3$node.label<-as.numeric(tree.rename3$node.label)
tree.rename3$node.label<-round(tree.rename3$node.label,digits=2)

### set some colours for the islands
r.col<-c("#0000ff",
         "#8a2be2",
         "#deb887",
         "#ff4500",
         "#006400",
         "#1e90ff",
         "#0000cd",
         "#bdb76b",
         "#228b22",
         "#7b68ee",
         "#483d8b",
         "#e0ffff",
         "#000000")
 
### the next line uses match to perform the same function as vlookup in excel
tree.rename3$tip.label <- (name$Alt_label[match(tree.rename3$tip.label,name$Name)])
str(tree.rename3$tip.label)

### plot the tree with the character reconstruction mapped
pdf(file=paste(out,"RASP_biogeography.pdf",sep=""), 30, 30)
plot(tree.rename3, show.node.label=FALSE, label.offset=0.0, cex=2)
nodelabels(pie=Cratopus_anc, piecol=r.col, cex=0.5)
nodelabels(tree.rename3$node.label,adj=c(2,2),frame="none",
            col=ifelse(tree.rename3$node.label>0.9,"red",
                      ifelse(tree.rename3$node.label>=0.5 & tree.rename3$node.label<0.9,"blue","#0000ff00")),cex=1)
legend(x=0.005, y=25,
       legend=c("Madagascar",
                "Reunion",
                "Mauritius",
                "Rodrigues",
                "Europa",
                "Grande Glorieuse",
                "Moheli",
                "Anjouan",
                "Grande Comore",
                "Juan de Nova",
                "Aldabra",
                "Seychelles",
                "n/s"),
       fill=r.col,
       cex=1.5)

### make an offset for the axis as R won't draw it from the tip to the root. The offset is a negative starting point for the axis equivalent to the 
### round ing up we do at the root end of the axis i.e. if we round 4.79 Mya to 5 Mya then we need to offset by minus ~0.21Ma of distance measured in
### branch lengths. To do this we divide the root height by the root age and multiply by the difference between the oldest value on the axis and
### the oldest value on the tree.
offset<-(round_any(max(HPD$Median),0.5)-max(HPD$Median))*(max(nodeHeights(tree.rename3))/max(HPD$Median))

## put on a the correct axis
axis(side=1,
     cex.axis=1,
     padj=1,
     at=seq(-offset, max(nodeHeights(tree.rename3)), by=(max(nodeHeights(tree.rename3))+offset)/(round_any(max(HPD$Median),0.5)/0.5)),
     labels=seq(round_any(max(HPD$Median),0.5),0,by=-0.5))
dev.off()