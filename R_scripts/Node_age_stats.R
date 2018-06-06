rm(list=ls())

colonisation_means<-NULL

for(i in 1:1000){
  node_ages<-read.csv(paste("~/Desktop/Viper_holding/DECJX/colonisation_means_unstratified_",i,".csv",sep=""))
  colonisation_means<-rbind(colonisation_means,node_ages)
}

### test whether colonisations of Mauritius are significantly younger than colonisations of Reunion
t.test(colonisation_means$Reu_Mau,colonisation_means$Mau_Reu, paired=TRUE)

### test whether colonisation of Mauritius are younger than in situ speciation on Mauritius
t.test(colonisation_means$Reu_Mau,colonisation_means$Mau_insitu, paired=TRUE)

### test whether in situ speciation events on Reunion are significantly younger than colonisations of Reunion
t.test(colonisation_means$Reu_insitu,colonisation_means$Mau_Reu, paired=TRUE)

### test whether flightloss nodes are yonger than other nodes
t.test(colonisation_means$fl_age,colonisation_means$not_fl_age, paired=TRUE)

t.test(colonisation_means$Reu_Mau_count,colonisation_means$Mau_Reu_count, paired=TRUE)

mean(colonisation_means$Reu_Mau_count)
mean(colonisation_means$Mau_Reu_count)
