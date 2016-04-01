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

### read in the tree
my.tree<-read.nexus("Data/trees/Cratopus_MCC_tree.nex")
my.tree<-my.tree

### read in the node ages and heights from the vstat file from MrBayes - I have deleted all tip ages leaving only the root and internal nodes
HPD<-read.csv("Data/bipartition_ages.csv")

### read in the list of names
name<-read.csv("Data/all_names.csv")
### read.csv turns text into factors, this gets messy later when plotting
### so make it character data
name<-data.frame(lapply(name, as.character), stringsAsFactors=FALSE)

########################################################################################################
######################## Do some rough phylogeography ############################
########################################################################################################

### read in the island data
islands<-read.csv("Data/dist_ultra.csv", header=FALSE)
colnames(islands)<-c("sample","island")

### calculate the ML character reconstruction
Cratopus_anc<-ace(islands$island,my.tree,type="discrete")

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

### Convert node lables to numeric and round to two dp
my.tree$node.label<-as.numeric(my.tree$node.label)
my.tree$node.label<-round(my.tree$node.label,digits=2)

### the next line uses match to perform the same function as vlookup in excel
my.tree$tip.label <- (name$Alt_label[match(my.tree$tip.label,name$Name)])

## @ knitr MLtreeplot

### plot the tree with the character reconstruction mapped
plot(my.tree, show.node.label=FALSE, label.offset=0.0, cex=0.5)
### add biogeographic piecharts
nodelabels(pie = Cratopus_anc$lik.anc, piecol = r.col, cex = 0.75)
### make a legend
legend(x=0.01, y=10,
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
       cex=0.4)

#### make an offset for the axis as R won't draw it from the tip to the root. The offset is a negative starting point for the axis equivalent to the
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

## @ knitr dumpchunk

######################################################################################################################################
###################### Plot the tree to the output folder using a pdf dev #######################################
######################################################################################################################################
### plot the tree with the character reconstruction mapped
pdf(file="Diagrams/ML_biogeography_MCC.pdf", 30, 30)
plot(my.tree, show.node.label=FALSE, label.offset=0.0, cex=2)
### add biogeographic piecharts
nodelabels(pie = Cratopus_anc$lik.anc, piecol = r.col, cex = 0.75)
### add a legend
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

#### make an offset for the axis as R won't draw it from the tip to the root. The offset is a negative starting point for the axis equivalent to the
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

dev.off()
