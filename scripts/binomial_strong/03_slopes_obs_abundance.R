##################################################################################
## Get pop specific SNP slopes over the timeseries for "observed" data set
## Done for binomal strong snp set BF>10 in windows and BF>30 overall
## SNPs across each env varaible within the baselineset have been removed 
## Does not include 1000 permuations of random sampling for BF<0
## Also done for unique slopes across all env (histPop)
## Author Daniel Anstett
## 
## 
## Last Modified Oct 31, 2023
###################################################################################
#Functions

# Run GLM

glm_obs <-function(dfA,dfB,env_input){
  freq_env1_slope<-data.frame()
  counter<-1
  for (i in 3:dim(dfA)[2]){
    for(j in 1:12){
      chrA<-colnames(dfA)[i]
      chrB<-colnames(dfB)[i]
      
      popSNPA <- dfA %>% filter(Site==j) %>% select(Site,Year,all_of(chrA))
      popSNPB <- dfB %>% filter(Site==j) %>% select(Site,Year,all_of(chrB))
      
      colnames(popSNPA)[3]<-"snp_ID"
      colnames(popSNPB)[3]<-"snp_ID"
      snp_ab <- cbind(popSNPA$snp_ID,popSNPB$snp_ID)
      
      if(sum(!is.na(snp_ab[,1]/snp_ab[,2]))>=2 || sum(!is.na(snp_ab[,2]/snp_ab[,1]))>=2) {
        rSNP <- glm(snp_ab ~ popSNPA$Year, family = binomial)
        freq_env1_slope[counter,1]<-unique(popSNPA$Site)
        freq_env1_slope[counter,2]<-chrA
        freq_env1_slope[counter,3]<-rSNP$coefficients[2] 
        sum_snp <- summary(rSNP)
        freq_env1_slope[counter,4]<-sum_snp$coefficients[2,2]
        freq_env1_slope[counter,5]<-env_input
        freq_env1_slope[counter,6]<-"obs"
        counter<-counter+1
      } else {
        freq_env1_slope[counter,1]<-unique(popSNPA$Site)
        freq_env1_slope[counter,2]<-chrA
        freq_env1_slope[counter,3]<-NA
        freq_env1_slope[counter,4]<-NA
        freq_env1_slope[counter,5]<-env_input
        freq_env1_slope[counter,6]<-"obs"
        
        counter<-counter+1
      }
    }
  }
  colnames(freq_env1_slope)<-c("Site","snp_ID","Slope","SE","ENV","Type")
  return(freq_env1_slope)
}



###################################################################################
#Import libraries
library(tidyverse)
library(boot)

#Import timeseries frequencies
freqA_env1 <- read_csv("data/binomial_strong/freqA_env1.csv")
freqA_env2 <- read_csv("data/binomial_strong/freqA_env2.csv")
freqA_env3 <- read_csv("data/binomial_strong/freqA_env3.csv")
freqA_env4 <- read_csv("data/binomial_strong/freqA_env4.csv")
freqA_env5 <- read_csv("data/binomial_strong/freqA_env5.csv")
freqA_env6 <- read_csv("data/binomial_strong/freqA_env6.csv")
freqA_env7 <- read_csv("data/binomial_strong/freqA_env7.csv")
freqA_env8 <- read_csv("data/binomial_strong/freqA_env8.csv")
freqA_env9 <- read_csv("data/binomial_strong/freqA_env9.csv")

freqB_env1 <- read_csv("data/binomial_strong/freqB_env1.csv")
freqB_env2 <- read_csv("data/binomial_strong/freqB_env2.csv")
freqB_env3 <- read_csv("data/binomial_strong/freqB_env3.csv")
freqB_env4 <- read_csv("data/binomial_strong/freqB_env4.csv")
freqB_env5 <- read_csv("data/binomial_strong/freqB_env5.csv")
freqB_env6 <- read_csv("data/binomial_strong/freqB_env6.csv")
freqB_env7 <- read_csv("data/binomial_strong/freqB_env7.csv")
freqB_env8 <- read_csv("data/binomial_strong/freqB_env8.csv")
freqB_env9 <- read_csv("data/binomial_strong/freqB_env9.csv")


###################################################################################

## Single Loop

freq_env1_slope <- data.frame()
counter<-1
  for (i in 3:dim(freqA_env1)[2]){
    for(j in 1:12){
      chrA<-colnames(freqA_env1)[i]
      chrB<-colnames(freqB_env1)[i]
      
      popSNPA <- freqA_env1 %>% filter(Site==j) %>% select(Site,Year,all_of(chrA))
      popSNPB <- freqB_env1 %>% filter(Site==j) %>% select(Site,Year,all_of(chrB))
      
      colnames(popSNPA)[3]<-"snp_ID"
      colnames(popSNPB)[3]<-"snp_ID"
      snp_ab <- cbind(popSNPA$snp_ID,popSNPB$snp_ID)

      if(sum(!is.na(snp_ab[,1]/snp_ab[,2]))>=2 || sum(!is.na(snp_ab[,2]/snp_ab[,1]))>=2) {
      rSNP <- glm(snp_ab ~ popSNPA$Year, family = binomial)
      freq_env1_slope[counter,1]<-unique(popSNPA$Site)
      freq_env1_slope[counter,2]<-chrA
      freq_env1_slope[counter,3]<-rSNP$coefficients[2]
      sum_snp <- summary(rSNP)
      freq_env1_slope[counter,4]<-sum_snp$coefficients[2,2]
      freq_env1_slope[counter,5]<-"MAT"
      freq_env1_slope[counter,6]<-"obs"
      counter<-counter+1
      } else {
      freq_env1_slope[counter,1]<-unique(popSNPA$Site)
      freq_env1_slope[counter,2]<-chrA
      freq_env1_slope[counter,3]<-NA
      freq_env1_slope[counter,4]<-NA
      freq_env1_slope[counter,5]<-"MAT"
      freq_env1_slope[counter,6]<-"obs"
      
      counter<-counter+1
    }
      
    }
  }
  colnames(freq_env1_slope)<-c("Site","snp_ID","Slope","SE","ENV","Type")


###################################################################################
###################################################################################
#Calc glm per env
slope_env1 <- glm_obs(freqA_env1,freqB_env1,"MAT")
slope_env2 <- glm_obs(freqA_env2,freqB_env2,"MAP")
slope_env3 <- glm_obs(freqA_env3,freqB_env3,"PAS")
slope_env4 <- glm_obs(freqA_env4,freqB_env4,"EXT")
slope_env5 <- glm_obs(freqA_env5,freqB_env5,"CMD")
slope_env6 <- glm_obs(freqA_env6,freqB_env6,"Tave_wt")
slope_env7 <- glm_obs(freqA_env7,freqB_env7,"Tave_sm")
slope_env8 <- glm_obs(freqA_env8,freqB_env8,"PPT_wt")
slope_env9 <- glm_obs(freqA_env9,freqB_env9,"PPT_sm")

slope_env_all <- rbind(slope_env1,
                       slope_env2,
                       slope_env3,
                       slope_env4,
                       slope_env5,
                       slope_env6,
                       slope_env7,
                       slope_env8,
                       slope_env9)

#Get unique SNP for merged env

slope_env_all_unique <- data.frame()
for (i in 1:dim(slope_env_all)){
  slope.temp <- slope_env_all %>% filter(Site==i)
  slope.unique <- slope.temp [!duplicated(slope.temp [ , "snp_ID"]), ]
  slope_env_all_unique <- rbind(slope_env_all_unique,slope.unique)
}


###################################################################################
#Export
write_csv(slope_env1, "data/binomial_strong/slope_obs_env1.csv")
write_csv(slope_env2, "data/binomial_strong/slope_obs_env2.csv")
write_csv(slope_env3, "data/binomial_strong/slope_obs_env3.csv")
write_csv(slope_env4, "data/binomial_strong/slope_obs_env4.csv")
write_csv(slope_env5, "data/binomial_strong/slope_obs_env5.csv")
write_csv(slope_env6, "data/binomial_strong/slope_obs_env6.csv")
write_csv(slope_env7, "data/binomial_strong/slope_obs_env7.csv")
write_csv(slope_env8, "data/binomial_strong/slope_obs_env8.csv")
write_csv(slope_env9, "data/binomial_strong/slope_obs_env9.csv")
write_csv(slope_env_all, "data/binomial_strong/slope_obs_all.csv")
write_csv(slope_env_all_unique, "data/binomial_strong/slope_obs_all_unique.csv")




