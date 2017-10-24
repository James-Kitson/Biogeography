### Scritp for making and editing the supp materialtwo
### tree for phylogeny paper
## @knitr cartoontree
### Clear the workspace
rm(list=ls())

library(devtools) ## devtools must be installed
#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
#install_github("KlausVigo/phangorn")
#install_github("liamrevell/phytools")

### load libraries
library(ape)
library(phytools)
library(plyr)
library(phangorn)

#packageVersion("phytools")
#packageVersion("phangorn")

### read in the data
sample.data<-read.csv("Data/Supplementary_names.csv", stringsAsFactors = FALSE)
tree<-read.nexus("Data/trees/Supplementary_Trees/Supp2_reps.nex")

### cut the data down to only the samples used for the analysis
sample.data.subs<-sample.data[match(tree$tip.label,sample.data$Sample),]

### reduce this data set to unique samples per calde for the backbone tree
unique.samples<-sample.data.subs[!duplicated(sample.data.subs$clade), ]
discard.tips<-sample.data.subs$Sample[!(sample.data.subs$Sample %in% unique.samples$Sample)]
unique.tree<-drop.tip(tree, tip = discard.tips)

### make the value vectors for tranforamtion table
tip.label<-unique.samples$Sample
clade.label<-unique.samples$clade
depth<-sapply(tip.label,function(x,y) 0.5*y$edge.length[which(unique.tree$edge[,2]==which(y$tip.label== x))],y=unique.tree)

### make the transofrmation table
trans<-data.frame(tip.label, clade.label, depth)
tip.counts<-count(sample.data.subs, "clade")
trans$N<-tip.counts$freq[match(trans$clade.label,tip.counts$clade)]

### make a backbone tree from the transformation table and the unique samples tree
tt<-phylo.toBackbone(unique.tree, trans)

### define the clade colours by extracting them from the dataframe with sample metadata
clade.colours<-unique.samples[,7]

## @knitr plotcartoontree
###  plot the backbone tree
#pdf("Diagrams/tree_supp2_cartoon.pdf",8,12)
plot(tt, print.clade.size=TRUE,
     col=clade.colours,
     fixed.height=TRUE)
#dev.off()

## @knitr cladetrees
root.nodes<-NULL
### find the root node of each clade in the backbone tree and plot the subtrees using a for loop
clades<-subset(tip.counts,freq>1, select = clade)
for(i in 1:nrow(clades)){
  my.clade<-subset(sample.data.subs, clade==clades[i,], select= Sample)
  root.nodes[i]<-mrca.phylo(tree, match(my.clade$Sample, tree$tip))
}

### plot all the clade trees - all of this is now in the rmarkdown as you can't call a forloop in a chunk it appears
## @knitr plotcladetrees
#for(i in 1:length(root.nodes)){
#png(paste("Diagrams/clades",i,".png", sep=""),800,1200)
#  clade.tree<-extract.clade(tree,root.nodes[i])
#  clade.tree$tip.label<-sample.data.subs$name[match(clade.tree$tip.label, sample.data.subs$Sample)]
#  plot(clade.tree,
#       type=ifelse(length(clade.tree$tip.label)>100, "fan", "phylogram"),
#       cex=ifelse(length(clade.tree$tip.label)>15, 0.5, 1))
#dev.off()
#}
