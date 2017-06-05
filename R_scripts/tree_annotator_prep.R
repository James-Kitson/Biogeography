rm(list=ls())

setwd("~/Documents/Projects/Cratopus/Biogeography/MrBayes/Tree_files/MrBayes_treesets/")

library(ape)

my_files<-list.files(pattern = "*GTRIG.nex.run.{2,3}t")

for(i in 1:length(my_files)){
trees_in<-read.nexus(file=my_files[i])
write.nexus(trees_in[2501:10001], file=paste("run",i,"_out.t", sep=''))
}

my_MCC_trees<-list.files(pattern = "*MCC")
