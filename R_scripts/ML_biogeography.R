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
my.trees<-read.nexus("Data/All_dating_mcorrected.nex.con.tre")
my.tree<-my.trees[[1]]

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

### resolve polytomies
tree.resolve.out<-multi2di(my.tree, random=TRUE)

### write the tree to an ouptut file. The polytomies are resolved randomly so you need to write out a file and use that in future
#write.nexus(tree.resolve.out, file = "Data/resolved_tree.nex")

### read in the resolved tree
tree.resolve<-read.nexus("Data/resolved_tree.nex")

### give zero length branches an arbitrarily short length
tree.resolve$edge.length<-ifelse(tree.resolve$edge.length==0,0.0000000001,tree.resolve$edge.length)
### check it has worked
tree.resolve$edge.length

### read in the island data
islands<-read.csv("DATA/dist_ultra.csv", header=FALSE)
colnames(islands)<-c("sample","island")

### calculate the ML character reconstruction
Cratopus_anc<-ace(islands$island,tree.resolve,type="discrete")

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
tree.resolve$tip.label <- (name$Alt_label[match(tree.resolve$tip.label,name$Name)])
str(tree.resolve$tip.label)

## @ knitr MLtreeplot

### plot the tree with the character reconstruction mapped
#pdf(file="Diagrams/ML_biogeography.pdf", 30, 30)
plot(tree.resolve, show.node.label=FALSE, label.offset=0.0, cex=2)
### add biogeographic piecharts
nodelabels(pie = Cratopus_anc$lik.anc, piecol = r.col, cex = 0.75)
### add bayesian support values
nodelabels(tree.resolve$node.label,adj=c(2,2),frame="none",
           col=ifelse(tree.resolve$node.label>0.9,"red",
                      ifelse(tree.resolve$node.label>=0.5 & tree.resolve$node.label<0.9,"blue","#0000ff00")),cex=1)
### add node numbers
tree.resolve$node.number<-seq(1, length(tree.resolve$node.label),1)
nodelabels(tree.resolve$node.number,adj=c(3.2,3.2),frame="none",col="black")

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
offset<-(round_any(max(HPD$Median),0.5)-max(HPD$Median))*(max(nodeHeights(tree.resolve))/max(HPD$Median))

## put on a the correct axis
axis(side=1,
     cex.axis=1,
     padj=1,
     at=seq(-offset, max(nodeHeights(tree.resolve)), by=(max(nodeHeights(tree.resolve))+offset)/(round_any(max(HPD$Median),0.5)/0.5)),
     labels=seq(round_any(max(HPD$Median),0.5),0,by=-0.5))

#dev.off()
