########################################################################################################
######################## Script for plotting Cratopus trees ############################
########################################################################################################

## @knitr BBMtreeprocess

### Clear the workspace
rm(list=ls())

### set working directory ####
setwd("~/Documents/Projects/Cratopus project/R analyses/Biogeography")
out<-paste(getwd(),"/Diagrams/",sep="")

### open APE
library(ape)
library(phytools)
library(plyr)

### read in the tree
my.trees<-read.nexus("Data/All_dating_mcorrected.nex.con.tre")
my.tree<-my.trees[[1]]

### read in the node ages and heights from the vstat file from MrBayes - I have deleted all tip ages leaving only the root and internal nodes
HPD<-read.csv("Data/bipartition_ages.csv")

### read in the list of names
name<-read.csv("Data/all_names.csv")
### read.csv turns text into factors, this gets messy later when plotting
### so make it character data
name<-data.frame(lapply(name, as.character), stringsAsFactors=FALSE)

###Get the biogeographic reconstruvtion from RASP BBM
Cratopus_anc<-read.csv("Data/RASP_anc.csv", row.names=1, header=TRUE)

########################################################################################################
######################## Process the biogeographic output and plot the tree ############################
########################################################################################################

### replace all the islands that are less than 5% probability with zero
Cratopus_anc[,2:13][Cratopus_anc[,2:13] < 5] <- 0

### Sum all the islands probabilities again with n/s vaues zeroed
Cratopus_anc$s<-rowSums(Cratopus_anc[,c(2:14)])

## Drop total
Cratopus_anc<-Cratopus_anc[,c(2:14)]

Cratopus_anc<-as.matrix(Cratopus_anc)

### Convert node lables to numeric and round to two dp
my.tree$node.label<-as.numeric(my.tree$node.label)
my.tree$node.label<-round(my.tree$node.label,digits=2)

### the next line uses match to perform the same function as vlookup in excel and matches names from the reference file to the tree tips
my.tree$tip.label <- (name$Alt_label[match(my.tree$tip.label,name$Name)])

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

## @knitr BBMtreeplot

### plot the tree with the character reconstruction mapped
#pdf(file=paste(out,"RASpBBM_biogeography.pdf",sep=""), 30, 30)
plot(my.tree, show.node.label=FALSE, label.offset=0.0, cex=1)
nodelabels(pie=Cratopus_anc, piecol=r.col, cex=0.5)
### add bayesian support values
nodelabels(my.tree$node.label,adj=c(2,2),frame="none",
            col=ifelse(my.tree$node.label>0.9,"red",
                      ifelse(my.tree$node.label>=0.5 & my.tree$node.label<0.9,"blue","#0000ff00")),cex=1)
### add node numbers
my.tree$node.number<-seq(1, length(my.tree$node.label),1)
nodelabels(my.tree$node.number,adj=c(3.2,3.2),frame="none",col="black")
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
       cex=1)

### make an offset for the axis as R won't draw it from the tip to the root. The offset is a negative starting point for the axis equivalent to the
### round ing up we do at the root end of the axis i.e. if we round 4.79 Mya to 5 Mya then we need to offset by minus ~0.21Ma of distance measured in
### branch lengths. To do this we divide the root height by the root age and multiply by the difference between the oldest value on the axis and
### the oldest value on the tree.
offset<-(round_any(max(HPD$Median),0.5)-max(HPD$Median))*(max(nodeHeights(my.tree))/max(HPD$Median))

## put on a the correct axis
axis(side=1,
     cex.axis=1,
     padj=1,
     at=seq(-offset, max(nodeHeights(my.tree)), by=(max(nodeHeights(my.tree))+offset)/(round_any(max(HPD$Median),0.5)/0.5)),
     labels=seq(round_any(max(HPD$Median),0.5),0,by=-0.5))
#dev.off()