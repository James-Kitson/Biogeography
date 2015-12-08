###########################################################################
############# Script for renaming and plotting trees in R #################
########## Requires MrBayes to use conformat=simple and a csv of ##########
############## the OTU codes and corresponding real names #################
###########################################################################

### Clear the workspace
rm(list=ls())

### set working directory and make objects for calling the string later
setwd("~/Documents/Projects/Crocodile project")
mod<-"/March 2015 matrix runs/Thallatosuchia_calibrated/Fixed_Thallato"
title<-"Fixed_Thallato_calibrated_relaxedclock"
work<-as.character(getwd())
out<-paste(work,mod,"/diagrams/",sep="")

### load the ape and phytools libraries
library(ape)
library(phytools)

### make a list of all files in the tree input directory
inputs<-list.files(paste(work,mod,sep=""))

### subset this to include only files that are tre files
inputs<-subset(inputs,grepl('.nex.con.tre$',inputs))

### make sure its ok
inputs

### extract a list of names from inputs for plotting file names
file.names<-substr(inputs,1,nchar(inputs)-nchar(".nex.con.tre"))
file.names

### List control command, change this from 1-n files and rerun each time
x<-1
###########################################################################
######################## Get the data #####################################
###########################################################################


### Load the tree from the working directory
Bayes.tree<-read.nexus(paste(work,mod,"/",inputs[[x]], sep=""))
### conformat=simple in MrBayes produces two trees, one with branch lengths
### and one with node probabilities, all the modification needs to happen
### on the first one so make than an object in its own right
tree<-Bayes.tree[[1]]
tree<-ladderize(tree)
tree<-root(tree,1)
### read in the list of names
name<-read.csv("~/Documents/Projects/Crocodile project/inputs/Marknames_March2015.csv")
### read.csv turns text into factors, this gets messy later when plotting
### so make it character data
name<-data.frame(lapply(name, as.character), stringsAsFactors=FALSE)
########################################################

### replace the codes with informative names
tree.rename<-tree
### the next line uses match to perform the same function as vlookup in excel
tree.rename$tip.label <- (name$name[match(tree.rename$tip.label,name$label)])

### COnvert node lables to numeric and round to two dp
tree.rename$node.label<-as.numeric(tree.rename$node.label)
tree.rename$node.label<-round(tree.rename$node.label,digits=2)

### check all the various labels are in the correct format
str(tree.rename$tip.label)
str(tree$tip.label)
str(tree.rename$node.label)

############## fancy node labels colour coded by value ###################

### make a new vector of the same length as the node label vector
fill <- character(length(tree.rename$node.label)) 

### values greater than or equal to 0.95 are called red
fill[tree.rename$node.label >= 0.95] <- "red"

### values less than 0.95 but greater than 0.75 are called blue
fill[tree.rename$node.label < 0.95 & tree.rename$node.label >= 0.75] <- "blue"

### all "NA" values (i.e. the root) and values less than 0.75 are coded as #0000ff00
### (the last two digets are the opacity resulting in no point for the root)
fill[tree.rename$node.label < 0.75] <- "#0000ff00"
fill[is.na(tree.rename$node.label)]<- "#0000ff00"

### check the vector looks right, this is then used in the nodelabels argument during plotting
fill

### it looks pretty rubbish to just have the colour dot with no edge so I need to specify that
### all the red and blue labels should have a black outline while all the clear ones have none
### I've done this the same way as above but it could be done better I'm sure.
stroke <- character(length(tree.rename$node.label)) 

### values greater than or equal to 0.95 are called red
stroke[tree.rename$node.label >= 0.95] <- "black"

### values less than 0.95 but greater than 0.75 are called blue
stroke[tree.rename$node.label < 0.95 & tree.rename$node.label >= 0.75] <- "black"

### all "NA" values (i.e. the root) and values less than 0.75 are coded as #0000ff00
### (the last two digets are the opacity resulting in no point for the root)
stroke[tree.rename$node.label < 0.75] <- "#0000ff00"
stroke[is.na(tree.rename$node.label)]<- "#0000ff00"

#########################################################
### Plot the main combined tree with support values
pdf(file=paste(out,file.names[[x]],"_",title,".pdf",sep=""), 30, 30)
###svg(file=paste(out,file.names[[x]],"_SD_02_values.svg",sep=""), 30, 30)
plot(tree.rename,
     show.node.label=FALSE,
     cex=2,
     x.lim=0.5,
     label.offset=0.001)
nodelabels(tree.rename$node.label,adj=c(1.1,1.3),frame="none",
           col=ifelse(tree.rename$node.label>0.9,"red",
                      ifelse(tree.rename$node.label>=0.75 & tree.rename$node.label<0.9,"blue","#0000ff00")),cex=1)

### make an offset for the axis as R won't draw it from the tip to the root. The offset is a negative starting point for the axis equivalent to the 
### round ing up we do at the root end of the axis i.e. if we round 227 Mya to 230 Mya then we need to offset by minus 3Ma of distance measured in
### branch lengths. To do this we divide the root height by the root age and multiply by -3.
offset<-3*(max(nodeHeights(tree.rename)/227))

### put on the titel for reference
title(main=paste(title,sep=""))

## put on a the correct axis
axis(side=1,cex.axis=1,padj=1,at=seq(-offset,max(nodeHeights(tree.rename)),by=max(nodeHeights(tree.rename))/23), labels=seq(230,0,by=-10))
abline(v=max(nodeHeights(tree.rename)))
abline(v=max(nodeHeights(tree.rename))-offset, col="blue")
abline(v=-offset, col="green")
abline(v=0,col="red")

###add.scale.bar(x=max(nodeHeights(tree.rename))/2, y=0.5)

dev.off()

max(nodeHeights(tree.rename))

### Plot the main combined tree with support colour codes
pdf(file=paste(out,file.names[[x]],"_numbers.pdf",sep=""),25,25)
plot(tree.rename,
     show.node.label=FALSE,
     cex=1.5)
nodelabels(seq(from=1,to=tree.rename$Nnode),adj=c(1,1),frame="none",col="red",cex=1.5)
title(main=paste(file.names[[x]]),
      cex.main=3)
dev.off()

