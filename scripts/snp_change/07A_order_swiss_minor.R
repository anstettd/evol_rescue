##################################################################################
## Reorder non-climate associated SNPs so that minor allele goes first
## Minor allele is defined by the baseline dataset
## Split into 100k chucks to make managable
## Author Daniel Anstett
## 
## 
## Last Modified Dec 30, 2023
###################################################################################
#Import libraries
library(tidyverse)

###################################################################################
#Functions

#Calibrate timeseries abundance so that A is always the minor allele

minor_A <- function(df,baseline){
  swiss_minor_1 <- data.frame()
  for(i in 1:dim(df)[1]){
    
    start_swiss <- baseline[i,] %>% select(-chr_snp)
    A_col <- sum(start_swiss[,seq(1, ncol(start_swiss), by = 2)])
    B_col <- sum(start_swiss[,seq(2, ncol(start_swiss), by = 2)])   
    swiss_minor_1[i,1]<- df[i,1]
    
    for(j in 0:61){
      if(A_col >= B_col){
        #B is the minor
        swiss_minor_1[i,2+(2*j)] <- df[i,3+(2*j)]
        swiss_minor_1[i,3+(2*j)] <- df[i,2+(2*j)]
      } else{
        swiss_minor_1[i,2+(2*j)] <- df[i,2+(2*j)]
        swiss_minor_1[i,3+(2*j)] <- df[i,3+(2*j)]
      }
    }
    #print(i)
  }
  colnames(swiss_minor_1) <- colnames(df)
  return(swiss_minor_1)
}

#################################################################################################
## Manupulate entire SNP datatset for timeseries

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

#Make pop order to organize site/year headers
#pop_order_2 <- data.frame()
#pop_order_2 [1,1] <- "chr_snp"
#pop_order_2 <- rbind(pop_order_2,pop_order)
#colnames(pop_order_2)<-"chr_snp"

#Import Baseline data
#import full snp table
pop_order_base<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", 
                           header=F, sep="\t")
snp_base<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", 
                     header=F, sep=" ")
loci_base<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", 
                      header=F, sep="\t")
colnames(loci_base) <- c("Chromosome","SNP")
loci_united_base <- loci_base %>% unite(chr_snp,"Chromosome","SNP",sep="_")
loci_snp_base <-cbind(loci_united_base,snp_base) #add snp lables to rows

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

#Filter baseline by snp_swiss (reduced timeseries)
base_swiss <-loci_snp_base %>% filter (chr_snp %in% as.character(snp_swiss$chr_snp))

#Filter snp_swiss (reduced timeseries) by baseline
snp_swiss <- snp_swiss %>% filter (chr_snp %in% as.character(base_swiss$chr_snp))


tail(base_swiss)[,1]
tail(snp_swiss)[,1]

#Clear memory
rm(loci)
rm(loci_snp)
rm(snp)
rm(loci_united)
rm(env1_bf0)
rm(env2_bf0)
rm(env3_bf0)
rm(env4_bf0)
rm(env5_bf0)
rm(env6_bf0)
rm(env7_bf0)
rm(env8_bf0)
rm(env9_bf0)
rm(pop_order)



#################################################################################################

#Toy example to get minor allele per snp/pop/time

#swiss_1 <- snp_swiss[1:10,]
#base_1 <- base_swiss[1:10,]

#swiss_minor_1 <- data.frame()

#for(i in 1:10){
  
  #Get minror allele 
#  start_swiss <- base_1[i,] %>% select(-chr_snp)
#  A_col <- sum(start_swiss[,seq(1, ncol(start_swiss), by = 2)])
#  B_col <- sum(start_swiss[,seq(2, ncol(start_swiss), by = 2)])   
#  swiss_minor_1[i,1]<- swiss_1[i,1]
  
#  for(j in 0:61){
#    if(A_col >= B_col){
#      #B is the minor
#      swiss_minor_1[i,2+(2*j)] <- swiss_1[i,3+(2*j)]
#      swiss_minor_1[i,3+(2*j)] <- swiss_1[i,2+(2*j)]
#    } else{
#      swiss_minor_1[i,2+(2*j)] <- swiss_1[i,2+(2*j)]
#      swiss_minor_1[i,3+(2*j)] <- swiss_1[i,3+(2*j)]
#    }
#  }
#}

#swiss_abA <-  abA(swiss_1,pop_order_2)
#swiss_abB <-  abB(swiss_1,pop_order_2)


#################################################################################################
#Make position A the minor allele frequency


swiss_1 <- snp_swiss[1:100000,] #Split data set into 100k ch
base_1 <- base_swiss[1:100000,]
#swiss_1 <- snp_swiss[1:1000,] #Split data set into 100k ch
#base_1 <- base_swiss[1:1000,]
swiss_1_minor <- minor_A(swiss_1,base_1)
write_csv(swiss_1,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_1.csv")
rm(swiss_1)
rm(base_1)


swiss_2 <- snp_swiss[100001:200000,] 
base_2 <- base_swiss[100001:200000,] 
swiss_2 <- minor_A(swiss_2,base_2)
write_csv(swiss_2,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_2.csv")
rm(swiss_2)
rm(base_2)


swiss_3 <- snp_swiss[200001:300000,] 
base_3 <- base_swiss[200001:300000,] 
swiss_3 <- minor_A(swiss_3,base_3)
write_csv(swiss_3,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_3.csv")
rm(swiss_3)
rm(base_3)


swiss_4 <- snp_swiss[300001:400000,] 
ase_4 <- base_swiss[300001:400000,] 
swiss_4 <- minor_A(swiss_4,base_4)
write_csv(swiss_4,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_4.csv")
rm(swiss_4)
rm(base_4)


swiss_5 <- snp_swiss[400001:500000,] 
base_5 <- base_swiss[400001:500000,] 
swiss_5 <- minor_A(swiss_5,base_5)
write_csv(swiss_5,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_5.csv")
rm(swiss_5)
rm(base_5)

swiss_6 <- snp_swiss[500001:600000,] 
base_6 <- base_swiss[500001:600000,] 
swiss_6 <- minor_A(swiss_6,base_6)
write_csv(swiss_6,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_6.csv")
rm(swiss_6)
rm(base_6)


swiss_7 <- snp_swiss[600001:700000,] 
base_7 <- base_swiss[600001:700000,] 
swiss_7 <- minor_A(swiss_7,base_7)
write_csv(swiss_7,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_7.csv")
rm(swiss_7)
rm(base_7)


swiss_8 <- snp_swiss[700001:800000,] 
base_8 <- base_swiss[700001:800000,] 
swiss_8 <- minor_A(swiss_8,base_8)
write_csv(swiss_8,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_8.csv")
rm(swiss_8)
rm(base_8)

swiss_9 <- snp_swiss[800001:900000,] 
base_9 <- base_swiss[800001:900000,] 
swiss_9 <- minor_A(swiss_9,base_9)
write_csv(swiss_9,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_9.csv")
rm(swiss_9)
rm(base_9)


swiss_10 <- snp_swiss[900001:1000000,] 
base_10 <- base_swiss[900001:1000000,] 
swiss_10 <- minor_A(swiss_10,base_10)
write_csv(swiss_10,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_10.csv")
rm(swiss_10)
rm(base_10)


swiss_11 <- snp_swiss[1000001:1100000,] 
base_11 <- base_swiss[1000001:1100000,] 
swiss_11 <- minor_A(swiss_11,base_11)
write_csv(swiss_11,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_11.csv")
rm(swiss_11)
rm(base_11)


swiss_12 <- snp_swiss[1100001:1200000,] 
base_12 <- base_swiss[1100001:1200000,] 
swiss_12 <- minor_A(swiss_12,base_12)
write_csv(swiss_12,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_12.csv")
rm(swiss_12)
rm(base_12)

swiss_13 <- snp_swiss[1200001:1300000,] 
base_13 <- base_swiss[1200001:1300000,] 
swiss_13 <- minor_A(swiss_13,base_13)
write_csv(swiss_13,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_13.csv")
rm(swiss_13)
rm(base_13)

swiss_14 <- snp_swiss[1300001:1400000,] 
base_14 <- base_swiss[1300001:1400000,] 
swiss_14 <- minor_A(swiss_14,base_14)
write_csv(swiss_14,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_14.csv")
rm(swiss_14)
rm(base_14)


swiss_15 <- snp_swiss[1400001:1500000,] 
base_15 <- base_swiss[1400001:1500000,] 
swiss_15 <- minor_A(swiss_15,base_15)
write_csv(swiss_15,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_15.csv")
rm(swiss_15)
rm(base_15)

swiss_16 <- snp_swiss[1500001:1600000,] 
base_16 <- base_swiss[1500001:1600000,] 
swiss_16 <- minor_A(swiss_16,base_16)
write_csv(swiss_16,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_16.csv")
rm(swiss_16)
rm(base_16)


swiss_17 <- snp_swiss[1600001:1607008,] 
base_17 <- base_swiss[1600001:1607008,] 
swiss_17 <- minor_A(swiss_17,base_17)
write_csv(swiss_17,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_minor_17.csv")
rm(swiss_17)
rm(base_17)











