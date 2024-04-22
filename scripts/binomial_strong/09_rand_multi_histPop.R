##################################################################################
## Generate stratified distribution of random non-climate associated slopes
## For unique SNP slopes across all env 
## Author Daniel Anstett
## 
## 
## Last Modified May 21, 2023
###################################################################################
#Import libraries
library(tidyverse)

#################################################################################################
#Functions for generating stratigied random distribution
#################################################################################################
#Large
#pop_snp = p1A
#freq_count = freq_count_env1_IDX
#Randomly stratified sampling is carried out for each frequency bin based on basetime starting conditions
large <- function(pop_snp,freq_count){
  rand_pop <- data.frame()
  bin_fraction<-1/length(freq_count)
  bin_step<-seq(0,1,bin_fraction)
  
  for (i in 1: length(freq_count)){
    if(i==length(freq_count)){
      p1_n<-pop_snp %>% filter(snpA>=bin_step[i] & snpA <= bin_step[i+1])
      list_env_p1_n<- sample.int(dim(p1_n)[1],freq_count[i],replace = FALSE)
      
    }else{
      p1_n<-pop_snp %>% filter(snpA>=bin_step[i] & snpA < bin_step[i+1])
      list_env_p1_n<- sample.int(dim(p1_n)[1],freq_count[i],replace = FALSE)
      
    }
   
    for (j in 1:length(list_env_p1_n)){
      rand_pop <- rbind(rand_pop ,p1_n %>% 
                          filter(chr_snp==as.character(p1_n$chr_snp[list_env_p1_n[j]])))
    }
    
  }

  return(rand_pop)
}

##Freq Count functions

#Diagnostic tool
#Used to diagnose what bin size to assign to each population. 
#Use function on data to generate bin table (freq_count_env1).
#Diagnose bin size per pop by taking into account the number and placement of zeros
freq_bins <- function(basetime){
  freq_count_calc<-data.frame()
  bin_size<-10
  for (i in 1:1){
    
#    if(i==1 || i==7){
#      bin_size<-10
#    }else if (i==3 || i==4 || i==6 || i==11){
#      bin_size<-4
#    }else if (i==8){
#      bin_size<-2
#    }else if(i==2 || i==5 || i==9 || i==10 || i==12){
#      bin_size<-5
#    }
    
    bin_fraction<-1/bin_size
    bin_step<-seq(0,1,bin_fraction)
    
    for(j in 1:bin_size){
      test_ENV <- basetime %>% 
        #filter(Site==i) %>% 
        select (-Site,-Year)
      test_ENV <- as.data.frame(test_ENV)
      test_ENV <- as.numeric(test_ENV[1,])
      if (j==bin_size){
        freq_count_calc[j,i]<-sum(test_ENV >= bin_step[j] & test_ENV <= bin_step[j+1],na.rm=T )
      }else{
        freq_count_calc[j,i]<-sum(test_ENV >= bin_step[j] & test_ENV < bin_step[j+1],na.rm=T )
      }  
    }
  }
  return(freq_count_calc)
}

#Specific bin function for pop 8 which has no values at one 0.1 increment
freq_bins_8 <- function(basetime){
  freq_count_calc<-data.frame()
  bin_size<-5
  for (i in 1:1){
    
    #    if(i==1 || i==7){
    #      bin_size<-10
    #    }else if (i==3 || i==4 || i==6 || i==11){
    #      bin_size<-4
    #    }else if (i==8){
    #      bin_size<-2
    #    }else if(i==2 || i==5 || i==9 || i==10 || i==12){
    #      bin_size<-5
    #    }
    
    bin_fraction<-1/bin_size
    bin_step<-seq(0,1,bin_fraction)
    
    for(j in 1:bin_size){
      test_ENV <- basetime %>% 
        #filter(Site==i) %>% 
        select (-Site,-Year)
      test_ENV <- as.data.frame(test_ENV)
      test_ENV <- as.numeric(test_ENV[1,])
      if (j==bin_size){
        freq_count_calc[j,i]<-sum(test_ENV >= bin_step[j] & test_ENV <= bin_step[j+1],na.rm=T )
      }else{
        freq_count_calc[j,i]<-sum(test_ENV >= bin_step[j] & test_ENV < bin_step[j+1],na.rm=T )
      }  
    }
  }
  return(freq_count_calc)
}



#Remove duplicated SNPs (dulicated across env) per population
snp_pop_unique <- function(Site_ID,env1,env2,env3,env4,env5,env6,env7,env8,env9){
  basetime_p_env1 <- basetime_env1 %>% filter(Site==Site_ID)
  basetime_p_env2 <- basetime_env2 %>% filter(Site==Site_ID)
  basetime_p_env3 <- basetime_env3 %>% filter(Site==Site_ID)
  basetime_p_env4 <- basetime_env4 %>% filter(Site==Site_ID)
  basetime_p_env5 <- basetime_env5 %>% filter(Site==Site_ID)
  basetime_p_env6 <- basetime_env6 %>% filter(Site==Site_ID)
  basetime_p_env7 <- basetime_env7 %>% filter(Site==Site_ID)
  basetime_p_env8 <- basetime_env8 %>% filter(Site==Site_ID)
  basetime_p_env9 <- basetime_env9 %>% filter(Site==Site_ID)
  
  basetime_joint_p <- cbind(basetime_p_env1,
                            basetime_p_env2,
                            basetime_p_env3,
                            basetime_p_env4,
                            basetime_p_env5,
                            basetime_p_env6,
                            basetime_p_env7,
                            basetime_p_env8,
                            basetime_p_env9)
  
  basetime_joint_p_unique <- basetime_joint_p[, !duplicated(colnames(basetime_joint_p))]
  return(basetime_joint_p_unique)
}




#################################################################################################
#Import BF>0 baseline SNPs
env1_bf0 <- read_csv("data/env1_BF0.csv")
env2_bf0 <- read_csv("data/env2_BF0.csv")
env3_bf0 <- read_csv("data/env3_BF0.csv")
env4_bf0 <- read_csv("data/env4_BF0.csv")
env5_bf0 <- read_csv("data/env5_BF0.csv")
env6_bf0 <- read_csv("data/env6_BF0.csv")
env7_bf0 <- read_csv("data/env7_BF0.csv")
env8_bf0 <- read_csv("data/env8_BF0.csv")
env9_bf0 <- read_csv("data/env9_BF0.csv")


#Import full snp table for timeseries
pop_order<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")
snp<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", header=F, sep=" ")
loci<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", header=F, sep="\t")
colnames(loci) <- c("Chromosome","SNP")
loci_united <- loci %>% unite(chr_snp,"Chromosome","SNP",sep="_")
loci_snp <-cbind(loci_united,snp) #add snp lables to rows

#Add Chromsome ID to pop order
pop_order_2 <- data.frame()
pop_order_2 [1,1] <- "chr_shp"
pop_order_2 <- rbind(pop_order_2,pop_order)

#Import pop_order table with regional and V information
pop_order_V <- read_csv("data/pop_order_V.csv")

#Filter full snp table to remove climate associated SNPs
snp_swiss <-loci_snp %>% filter (!chr_snp %in% as.character(env1_bf0$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(env2_bf0$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(env3_bf0$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(env4_bf0$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(env5_bf0$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(env6_bf0$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(env7_bf0$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(env8_bf0$chr_snp))
snp_swiss <-snp_swiss  %>% filter (!chr_snp %in% as.character(env9_bf0$chr_snp))

rm(snp)
rm(loci)
rm(loci_united)
rm(loci_snp)
rm(pop_order)
rm(env1_bf0)
rm(env2_bf0)
rm(env3_bf0)
rm(env4_bf0)
rm(env5_bf0)
rm(env6_bf0)
rm(env7_bf0)
rm(env8_bf0)
rm(env9_bf0)


#################################################################################################
# Import non-climate associated slopes with low SE
swiss_glm <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_filter.csv")

#Import obs (snp set) frequency
#already filtered to include only basetime snps with low SE
basetime_env1 <- read_csv("data/binomial_strong/env1_low.csv")
basetime_env2 <- read_csv("data/binomial_strong/env2_low.csv")
basetime_env3 <- read_csv("data/binomial_strong/env3_low.csv")
basetime_env4 <- read_csv("data/binomial_strong/env4_low.csv")
basetime_env5 <- read_csv("data/binomial_strong/env5_low.csv")
basetime_env6 <- read_csv("data/binomial_strong/env6_low.csv")
basetime_env7 <- read_csv("data/binomial_strong/env7_low.csv")
basetime_env8 <- read_csv("data/binomial_strong/env8_low.csv")
basetime_env9 <- read_csv("data/binomial_strong/env9_low.csv")


############################################################################################

#Filter out duplicated SNPs per population
#Out of function

#basetime_p_env1 <- basetime_env1 %>% filter(Site==1)
#basetime_p_env2 <- basetime_env2 %>% filter(Site==1)
#basetime_p_env3 <- basetime_env3 %>% filter(Site==1)
#basetime_p_env4 <- basetime_env4 %>% filter(Site==1)
#basetime_p_env5 <- basetime_env5 %>% filter(Site==1)
#basetime_p_env6 <- basetime_env6 %>% filter(Site==1)
#basetime_p_env7 <- basetime_env7 %>% filter(Site==1)
#basetime_p_env8 <- basetime_env8 %>% filter(Site==1)
#basetime_p_env9 <- basetime_env9 %>% filter(Site==1)

#basetime_joint_p <- cbind(basetime_p_env1,
#                           basetime_p_env2,
#                           basetime_p_env3,
#                           basetime_p_env4,
#                           basetime_p_env5,
#                           basetime_p_env6,
#                           basetime_p_env7,
#                           basetime_p_env8,
#                           basetime_p_env9)
  
#basetime_joint_p_unique <- basetime_joint_p[, !duplicated(colnames(basetime_joint_p))]
  
############################################################################################
#Apply filter duplicated function  

basetime_joint_p01 <- snp_pop_unique(1,basetime_env1,basetime_env2,basetime_env3,basetime_env4,basetime_env5,basetime_env6,basetime_env7,basetime_env8,basetime_env9)
basetime_joint_p02 <- snp_pop_unique(2,basetime_env1,basetime_env2,basetime_env3,basetime_env4,basetime_env5,basetime_env6,basetime_env7,basetime_env8,basetime_env9)
basetime_joint_p03 <- snp_pop_unique(3,basetime_env1,basetime_env2,basetime_env3,basetime_env4,basetime_env5,basetime_env6,basetime_env7,basetime_env8,basetime_env9)
basetime_joint_p04 <- snp_pop_unique(4,basetime_env1,basetime_env2,basetime_env3,basetime_env4,basetime_env5,basetime_env6,basetime_env7,basetime_env8,basetime_env9)
basetime_joint_p05 <- snp_pop_unique(5,basetime_env1,basetime_env2,basetime_env3,basetime_env4,basetime_env5,basetime_env6,basetime_env7,basetime_env8,basetime_env9)
basetime_joint_p06 <- snp_pop_unique(6,basetime_env1,basetime_env2,basetime_env3,basetime_env4,basetime_env5,basetime_env6,basetime_env7,basetime_env8,basetime_env9)
basetime_joint_p07 <- snp_pop_unique(7,basetime_env1,basetime_env2,basetime_env3,basetime_env4,basetime_env5,basetime_env6,basetime_env7,basetime_env8,basetime_env9)
basetime_joint_p08 <- snp_pop_unique(8,basetime_env1,basetime_env2,basetime_env3,basetime_env4,basetime_env5,basetime_env6,basetime_env7,basetime_env8,basetime_env9)
basetime_joint_p09 <- snp_pop_unique(9,basetime_env1,basetime_env2,basetime_env3,basetime_env4,basetime_env5,basetime_env6,basetime_env7,basetime_env8,basetime_env9)
basetime_joint_p10 <- snp_pop_unique(10,basetime_env1,basetime_env2,basetime_env3,basetime_env4,basetime_env5,basetime_env6,basetime_env7,basetime_env8,basetime_env9)
basetime_joint_p11 <- snp_pop_unique(11,basetime_env1,basetime_env2,basetime_env3,basetime_env4,basetime_env5,basetime_env6,basetime_env7,basetime_env8,basetime_env9)
basetime_joint_p12 <- snp_pop_unique(12,basetime_env1,basetime_env2,basetime_env3,basetime_env4,basetime_env5,basetime_env6,basetime_env7,basetime_env8,basetime_env9)


#Diagnose bin size per pop 
freq_count_env1_diag <- freq_bins(basetime_joint_p01)
freq_count_env2_diag <- freq_bins(basetime_joint_p02)
freq_count_env3_diag <- freq_bins(basetime_joint_p03)
freq_count_env4_diag <- freq_bins(basetime_joint_p04)
freq_count_env5_diag <- freq_bins(basetime_joint_p05)
freq_count_env6_diag <- freq_bins(basetime_joint_p06)
freq_count_env7_diag <- freq_bins(basetime_joint_p07)
freq_count_env8_diag <- freq_bins_8(basetime_joint_p08)
freq_count_env9_diag <- freq_bins(basetime_joint_p09)
freq_count_env10_diag <- freq_bins(basetime_joint_p10)
freq_count_env11_diag <- freq_bins(basetime_joint_p11)
freq_count_env12_diag <- freq_bins(basetime_joint_p12)

freq_count_diag <- cbind(freq_count_env1_diag,
                         freq_count_env2_diag,
                         freq_count_env3_diag,
                         freq_count_env4_diag,
                         freq_count_env5_diag,
                         freq_count_env6_diag,
                         freq_count_env7_diag,
                         freq_count_env8_diag,
                         freq_count_env9_diag,
                         freq_count_env10_diag,
                         freq_count_env11_diag,
                         freq_count_env12_diag)
freq_count_diag[6:10,8] <- NA
colnames(freq_count_diag) <- c("p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12")
freq_count_diag 



###################################################################################

#Diagnose bin size per pop per env variable
#Use to setup freq_bins_env function
#freq_count_env1_diag <- freq_bins(basetime_env1)
#freq_count_env1_diag

#Make frequency count table
#freq_count_env1 <- freq_bins_env(basetime_env1)
#confirm no zeros are present
#freq_count_env1

######################################################################################################
#Break up freq count into individual vectors
######################################################################################################
freq_count_env1_1 <- as.vector(freq_count_diag[,1])
freq_count_env1_2 <- as.vector(freq_count_diag[,2])
freq_count_env1_3 <- as.vector(freq_count_diag[,3])
freq_count_env1_4 <- as.vector(freq_count_diag[,4])
freq_count_env1_5 <- as.vector(freq_count_diag[,5])
freq_count_env1_6 <- as.vector(freq_count_diag[,6])
freq_count_env1_7 <- as.vector(freq_count_diag[,7])
freq_count_env1_8 <- as.vector(freq_count_diag[,8])
freq_count_env1_9 <- as.vector(freq_count_diag[,9])
freq_count_env1_10 <- as.vector(freq_count_diag[,10])
freq_count_env1_11 <- as.vector(freq_count_diag[,11])
freq_count_env1_12 <- as.vector(freq_count_diag[,12])

freq_count_env1_1<-freq_count_env1_1[!is.na(freq_count_env1_1)]
freq_count_env1_2<-freq_count_env1_2[!is.na(freq_count_env1_2)]
freq_count_env1_3<-freq_count_env1_3[!is.na(freq_count_env1_3)]
freq_count_env1_4<-freq_count_env1_4[!is.na(freq_count_env1_4)]
freq_count_env1_5<-freq_count_env1_5[!is.na(freq_count_env1_5)]
freq_count_env1_6<-freq_count_env1_6[!is.na(freq_count_env1_6)]
freq_count_env1_7<-freq_count_env1_7[!is.na(freq_count_env1_7)]
freq_count_env1_8<-freq_count_env1_8[!is.na(freq_count_env1_8)]
freq_count_env1_9<-freq_count_env1_9[!is.na(freq_count_env1_9)]
freq_count_env1_10<-freq_count_env1_10[!is.na(freq_count_env1_10)]
freq_count_env1_11<-freq_count_env1_11[!is.na(freq_count_env1_11)]
freq_count_env1_12<-freq_count_env1_12[!is.na(freq_count_env1_12)]

#######################################################################################################
#Make SNP table for neutral loci
#######################################################################################################

#Filter pop_order_V to call for 2010 regional basetime V column names
pop_V1 <- pop_order_V %>% filter(Year==2010 & Pop==1)
pop_V2 <- pop_order_V %>% filter(Year==2010 & Pop==2)
pop_V3 <- pop_order_V %>% filter(Year==2010 & Pop==3)
pop_V4 <- pop_order_V %>% filter(Year==2010 & Pop==4)
pop_V5 <- pop_order_V %>% filter(Year==2010 & Pop==5)
pop_V6 <- pop_order_V %>% filter(Year==2010 & Pop==6)

pop_V7 <- pop_order_V %>% filter(Year==2010 & Pop==7)
pop_V8 <- pop_order_V %>% filter(Year==2011 & Pop==8)
pop_V9 <- pop_order_V %>% filter(Year==2010 & Pop==9)
pop_V10 <- pop_order_V %>% filter(Year==2011 & Pop==10)
pop_V11 <- pop_order_V %>% filter(Year==2010 & Pop==11)
pop_V12 <- pop_order_V %>% filter(Year==2010 & Pop==12)

#Filter for neutral snps with low SE for each population
snp_swiss_1 <- snp_swiss %>% filter (chr_snp %in% as.character(filter(swiss_glm,Site==1)$snp_ID))
snp_swiss_2 <- snp_swiss %>% filter (chr_snp %in% as.character(filter(swiss_glm,Site==2)$snp_ID))
snp_swiss_3 <- snp_swiss %>% filter (chr_snp %in% as.character(filter(swiss_glm,Site==3)$snp_ID))
snp_swiss_4 <- snp_swiss %>% filter (chr_snp %in% as.character(filter(swiss_glm,Site==4)$snp_ID))
snp_swiss_5 <- snp_swiss %>% filter (chr_snp %in% as.character(filter(swiss_glm,Site==5)$snp_ID))
snp_swiss_6 <- snp_swiss %>% filter (chr_snp %in% as.character(filter(swiss_glm,Site==6)$snp_ID))
snp_swiss_7 <- snp_swiss %>% filter (chr_snp %in% as.character(filter(swiss_glm,Site==7)$snp_ID))
snp_swiss_8 <- snp_swiss %>% filter (chr_snp %in% as.character(filter(swiss_glm,Site==8)$snp_ID))
snp_swiss_9 <- snp_swiss %>% filter (chr_snp %in% as.character(filter(swiss_glm,Site==9)$snp_ID))
snp_swiss_10 <- snp_swiss %>% filter (chr_snp %in% as.character(filter(swiss_glm,Site==10)$snp_ID))
snp_swiss_11 <- snp_swiss %>% filter (chr_snp %in% as.character(filter(swiss_glm,Site==11)$snp_ID))
snp_swiss_12 <- snp_swiss %>% filter (chr_snp %in% as.character(filter(swiss_glm,Site==12)$snp_ID))

rm(snp_swiss)

#Filter Baseline A and B numbers for 12 basetime pops
loci_base_1 <- snp_swiss_1 %>% select(chr_snp,all_of(pop_V1$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))
loci_base_2 <- snp_swiss_2 %>% select(chr_snp,all_of(pop_V2$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))
loci_base_3 <- snp_swiss_3 %>% select(chr_snp,all_of(pop_V3$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))
loci_base_4 <- snp_swiss_4 %>% select(chr_snp,all_of(pop_V4$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))
loci_base_5 <- snp_swiss_5 %>% select(chr_snp,all_of(pop_V5$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))
loci_base_6 <- snp_swiss_6 %>% select(chr_snp,all_of(pop_V6$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))

loci_base_7 <- snp_swiss_7 %>% select(chr_snp,all_of(pop_V7$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))
loci_base_8 <- snp_swiss_8 %>% select(chr_snp,all_of(pop_V8$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))
loci_base_9 <- snp_swiss_9 %>% select(chr_snp,all_of(pop_V9$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))
loci_base_10 <- snp_swiss_10 %>% select(chr_snp,all_of(pop_V10$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))
loci_base_11 <- snp_swiss_11 %>% select(chr_snp,all_of(pop_V11$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))
loci_base_12 <- snp_swiss_12 %>% select(chr_snp,all_of(pop_V12$V_ID)) %>% mutate(region_sum= rowSums(across(where(is.numeric))))

rm(snp_swiss_1)
rm(snp_swiss_2)
rm(snp_swiss_3)
rm(snp_swiss_4)
rm(snp_swiss_5)
rm(snp_swiss_6)
rm(snp_swiss_7)
rm(snp_swiss_8)
rm(snp_swiss_9)
rm(snp_swiss_10)
rm(snp_swiss_11)
rm(snp_swiss_12)

#Calculate frequency for SNP A for 12 basetime pops
#Gives all basetime frequencies per population
p1A <- loci_base_1 %>% mutate(snpA=(loci_base_1[,2])/(region_sum)) %>% select(chr_snp,snpA)
p2A <- loci_base_2 %>% mutate(snpA=(loci_base_2[,2])/(region_sum)) %>% select(chr_snp,snpA)
p3A <- loci_base_3 %>% mutate(snpA=(loci_base_3[,2])/(region_sum)) %>% select(chr_snp,snpA)
p4A <- loci_base_4 %>% mutate(snpA=(loci_base_4[,2])/(region_sum)) %>% select(chr_snp,snpA)
p5A <- loci_base_5 %>% mutate(snpA=(loci_base_5[,2])/(region_sum)) %>% select(chr_snp,snpA)
p6A <- loci_base_6 %>% mutate(snpA=(loci_base_6[,2])/(region_sum)) %>% select(chr_snp,snpA)

p7A <- loci_base_7 %>% mutate(snpA=(loci_base_7[,2])/(region_sum)) %>% select(chr_snp,snpA)
p8A <- loci_base_8 %>% mutate(snpA=(loci_base_8[,2])/(region_sum)) %>% select(chr_snp,snpA)
p9A <- loci_base_9 %>% mutate(snpA=(loci_base_9[,2])/(region_sum)) %>% select(chr_snp,snpA)
p10A <- loci_base_10 %>% mutate(snpA=(loci_base_10[,2])/(region_sum)) %>% select(chr_snp,snpA)
p11A <- loci_base_11 %>% mutate(snpA=(loci_base_11[,2])/(region_sum)) %>% select(chr_snp,snpA)
p12A <- loci_base_12 %>% mutate(snpA=(loci_base_12[,2])/(region_sum)) %>% select(chr_snp,snpA)

rm(loci_base_1)
rm(loci_base_2)
rm(loci_base_3)
rm(loci_base_4)
rm(loci_base_5)
rm(loci_base_6)
rm(loci_base_7)
rm(loci_base_8)
rm(loci_base_9)
rm(loci_base_10)
rm(loci_base_11)
rm(loci_base_12)

#Plot Starting Frequency for all neutral sites with low SE
env_slope<- ggplot(p7A,aes(x=snpA))+
  geom_histogram(position="identity",bins=10,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Starting Frequency")+
  theme_classic()
env_slope











#Input data frame
rand_slope_out<-data.frame()

#######################################################################################################
# Stratified random sampling SNPs for a given population
for (seed_num in 1:1000){
set.seed(seed_num)

#Implement large function
rand_pop1 <- large(p1A,freq_count_env1_1)
rand_pop2 <- large(p2A,freq_count_env1_2)
rand_pop3 <- large(p3A,freq_count_env1_3)
rand_pop4 <- large(p4A,freq_count_env1_4)
rand_pop5 <- large(p5A,freq_count_env1_5)
rand_pop6 <- large(p6A,freq_count_env1_6)
rand_pop7 <- large(p7A,freq_count_env1_7)
rand_pop8 <- large(p8A,freq_count_env1_8)
rand_pop9 <- large(p9A,freq_count_env1_9)
rand_pop10 <- large(p10A,freq_count_env1_10)
rand_pop11 <- large(p11A,freq_count_env1_11)
rand_pop12 <- large(p12A,freq_count_env1_12)

#Get slopes for stratified permuted SNPs
rand_slope_1 <- filter(swiss_glm,Site==1) %>% filter (snp_ID %in% as.character(rand_pop1$chr_snp))
rand_slope_2 <- filter(swiss_glm,Site==2)  %>% filter (snp_ID %in% as.character(rand_pop2$chr_snp)) 
rand_slope_3 <- filter(swiss_glm,Site==3)  %>% filter (snp_ID %in% as.character(rand_pop3$chr_snp)) 
rand_slope_4 <- filter(swiss_glm,Site==4)  %>% filter (snp_ID %in% as.character(rand_pop4$chr_snp)) 
rand_slope_5 <- filter(swiss_glm,Site==5)  %>% filter (snp_ID %in% as.character(rand_pop5$chr_snp)) 
rand_slope_6 <- filter(swiss_glm,Site==6)  %>% filter (snp_ID %in% as.character(rand_pop6$chr_snp)) 
rand_slope_7 <- filter(swiss_glm,Site==7)  %>% filter (snp_ID %in% as.character(rand_pop7$chr_snp)) 
rand_slope_8 <- filter(swiss_glm,Site==8)  %>% filter (snp_ID %in% as.character(rand_pop8$chr_snp)) 
rand_slope_9 <- filter(swiss_glm,Site==9)  %>% filter (snp_ID %in% as.character(rand_pop9$chr_snp))
rand_slope_10 <- filter(swiss_glm,Site==10)  %>% filter (snp_ID %in% as.character(rand_pop10$chr_snp)) 
rand_slope_11 <- filter(swiss_glm,Site==11)  %>% filter (snp_ID %in% as.character(rand_pop11$chr_snp)) 
rand_slope_12 <- filter(swiss_glm,Site==12)  %>% filter (snp_ID %in% as.character(rand_pop12$chr_snp)) 


#Bind populations for each env
rand_slope_env1 <- rbind(rand_slope_1,
                      rand_slope_2,
                      rand_slope_3,
                      rand_slope_4,
                      rand_slope_5,
                      rand_slope_6,
                      rand_slope_7,
                      rand_slope_8,
                      rand_slope_9,
                      rand_slope_10,
                      rand_slope_11,
                      rand_slope_12)

rand_slope_env1 <- rand_slope_env1 %>% mutate (Seed_ID = seed_num)

#Export each joint df

rand_slope_out<-rbind(rand_slope_out,rand_slope_env1)

print(seed_num)

}

#Save large files in folder outside of github
setwd("~/Dropbox/AM_Workshop/Large_files")

write_csv(rand_slope_out, "rand_slope_histPop_strong.csv")

setwd("~/Dropbox/AM_Workshop/snp_change")















