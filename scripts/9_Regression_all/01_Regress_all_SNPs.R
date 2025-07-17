##################################################################################
## Calculate GLMs for all Timeseries SNPs 
## Irrespective of climate association or baseline presence
## Split into 100k chucks to make managable
## Author Daniel Anstett
## 
## 
## Last Modified 07042025
###################################################################################
#Import libraries
library(tidyverse)

###################################################################################
##Functions

#Generate frequency matrix for prop A 
abA <- function(snp_table,pop_ID) {
  snp_prop_A<- snp_table %>% select (chr_snp)
  counter=2
  pop_num=1
  for (i in seq(1,dim(snp_table)[2]-1,2)) {
    
    snpA<-paste("V",i, sep="") #sets up string for snpA for pop num
    snpB<-paste("V",i+1, sep="") #sets up string for snpA for pop num
    P<-paste("P", pop_num, sep="") #sets up paper_ID (population ID)
    tmp<-snp_table %>% select(snpA,snpB) #graps all snps per paper_ID
    
    colnames(tmp)<-c("A", "B") #renames column headers
    
    snp_prop_A[,counter]<-as.numeric(tmp$A) #calc proportion for A, assign to output
    colnames (snp_prop_A)[counter]<-P 
    
    counter<-counter+1
    pop_num<-pop_num+1
  }
  
  colnames(snp_prop_A)<- c("chr_snp", pop_ID[,1]) #name each pop/time combination
  rownames(snp_prop_A)<- snp_prop_A$chr_snp
  snp1A_T <- as.data.frame(t(snp_prop_A)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
  colnames(snp1A_T) <- snp1A_T[1,]
  snp1A_T <- snp1A_T[-1,]
  colnames(snp1A_T)[1]<- "Site"
  colnames(snp1A_T)[2]<- "Year"
  return(snp1A_T)
}



abB <- function(snp_table,pop_ID) {
  snp_prop_A<- snp_table %>% select (chr_snp)
  counter=2
  pop_num=1
  for (i in seq(1,dim(snp_table)[2]-1,2)) {
    
    snpA<-paste("V",i, sep="") #sets up string for snpA for pop num
    snpB<-paste("V",i+1, sep="") #sets up string for snpA for pop num
    P<-paste("P", pop_num, sep="") #sets up paper_ID (population ID)
    tmp<-snp_table %>% select(snpA,snpB) #graps all snps per paper_ID
    
    colnames(tmp)<-c("A", "B") #renames column headers
    
    snp_prop_A[,counter]<-as.numeric(tmp$B) #calc proportion for A, assign to output
    colnames (snp_prop_A)[counter]<-P 
    
    counter<-counter+1
    pop_num<-pop_num+1
  }
  
  colnames(snp_prop_A)<- c("chr_snp", pop_ID[,1]) #name each pop/time combination
  rownames(snp_prop_A)<- snp_prop_A$chr_snp
  snp1A_T <- as.data.frame(t(snp_prop_A)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
  colnames(snp1A_T) <- snp1A_T[1,]
  snp1A_T <- snp1A_T[-1,]
  colnames(snp1A_T)[1]<- "Site"
  colnames(snp1A_T)[2]<- "Year"
  return(snp1A_T)
}



# Melt glm Josee
#Get Slope and SE from GLM
slope_melt <- function(dfA,dfB) {
  freq_1.melted <- reshape2::melt(dfA, id.vars = c("Site", "Year"))
  freq_2.melted <- reshape2::melt(dfB, id.vars = c("Site", "Year"))
  colnames(freq_1.melted)[3] <- "snp_ID"
  colnames(freq_1.melted)[4] <- "abA"
  colnames(freq_2.melted)[3] <- "snp_ID"
  colnames(freq_2.melted)[4] <- "abB"
  ab.melt <- cbind(freq_1.melted,freq_2.melted$abB)
  colnames(ab.melt)[5] <- "abB"
  rm(freq_1.melted)
  rm(freq_2.melted)
  
  #Extract slope and SE from binomial glm removing all instances where there is not enough data
  freq_1.slope <- group_by(ab.melt, Site, snp_ID) %>%
    arrange(Site, Year) %>% 
    summarize(Slope = ifelse(sum(!is.na(as.numeric(abA) / as.numeric(abB))) >= 2,
                             glm(cbind(as.numeric(abA), as.numeric(abB)) ~ as.numeric(Year), family = binomial)$coefficients[2],
                             ifelse(sum(!is.na(as.numeric(abB) / as.numeric(abA))) >= 2,
                                    glm(cbind(as.numeric(abA), as.numeric(abB)) ~ as.numeric(Year), family = binomial)$coefficients[2],
                                    NA)),
              SE = ifelse(sum(!is.na(as.numeric(abA) / as.numeric(abB))) >= 2,
                          summary(glm(cbind(as.numeric(abA),as.numeric(abB)) ~ as.numeric(Year), family = binomial))$coefficients[2,2],
                          ifelse(sum(!is.na(as.numeric(abB) / as.numeric(abA))) >= 2,
                                 summary(glm(cbind(as.numeric(abA),as.numeric(abB)) ~ as.numeric(Year), family = binomial))$coefficients[2,2],
                                 NA)))
  return(freq_1.slope)
}


#################################################################################################
## Import

#Import full snp table for timeseries
pop_order<-read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")
snp<-read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", header=F, sep=" ")
loci<-read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", header=F, sep="\t")
colnames(loci) <- c("Chromosome","SNP")
loci_united <- loci %>% unite(chr_snp,"Chromosome","SNP",sep="_")
snp_all <-cbind(loci_united,snp) #add snp lables to rows
dim(snp_all)

#Import population order data
pop_order<-read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")

#Make pop order to organize site/year headers
pop_order_2 <- data.frame()
#pop_order_2 [1,1] <- "chr_snp"
pop_order_2 <- rbind(pop_order_2,pop_order)
colnames(pop_order_2)<-"chr_snp"


###################################################################################
#Calculate GLMS for 2,189,349 SNPs over the drought

snp_1 <- snp_all[1:100000,]
snp_abA_1 <-  abA(snp_1,pop_order_2)
snp_abB_1 <-  abB(snp_1,pop_order_2)
snp_glm_1 <- slope_melt(snp_abA_1,snp_abB_1) #Run glm function
write_csv(snp_glm_1,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_01.csv")
rm(snp_abA_1)
rm(snp_abB_1)
rm(snp_glm_1)
rm(snp_1)

snp_2 <- snp_all[100001:200000,]
snp_abA_2 <-  abA(snp_2,pop_order_2)
snp_abB_2 <-  abB(snp_2,pop_order_2)
snp_glm_2 <- slope_melt(snp_abA_2,snp_abB_2) #Run glm function
write_csv(snp_glm_2,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_02.csv")
rm(snp_abA_2)
rm(snp_abB_2)
rm(snp_glm_2)
rm(snp_2)

snp_3 <- snp_all[200001:300000,]
snp_abA_3 <-  abA(snp_3,pop_order_2)
snp_abB_3 <-  abB(snp_3,pop_order_2)
snp_glm_3 <- slope_melt(snp_abA_3,snp_abB_3) #Run glm function
write_csv(snp_glm_3,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_03.csv")
rm(snp_abA_3)
rm(snp_abB_3)
rm(snp_glm_3)
rm(snp_3)

snp_4 <- snp_all[300001:400000,]
snp_abA_4 <-  abA(snp_4,pop_order_2)
snp_abB_4 <-  abB(snp_4,pop_order_2)
snp_glm_4 <- slope_melt(snp_abA_4,snp_abB_4) #Run glm function
write_csv(snp_glm_4,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_04.csv")
rm(snp_abA_4)
rm(snp_abB_4)
rm(snp_glm_4)
rm(snp_4)

snp_5 <- snp_all[400001:500000,]
snp_abA_5 <-  abA(snp_5,pop_order_2)
snp_abB_5 <-  abB(snp_5,pop_order_2)
snp_glm_5 <- slope_melt(snp_abA_5,snp_abB_5)
write_csv(snp_glm_5,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_05.csv")
rm(snp_abA_5)
rm(snp_abB_5)
rm(snp_glm_5)
rm(snp_5)

snp_6 <- snp_all[500001:600000,]
snp_abA_6 <-  abA(snp_6,pop_order_2)
snp_abB_6 <-  abB(snp_6,pop_order_2)
snp_glm_6 <- slope_melt(snp_abA_6,snp_abB_6)
write_csv(snp_glm_6,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_06.csv")
rm(snp_abA_6)
rm(snp_abB_6)
rm(snp_glm_6)
rm(snp_6)

snp_7 <- snp_all[600001:700000,]
snp_abA_7 <-  abA(snp_7,pop_order_2)
snp_abB_7 <-  abB(snp_7,pop_order_2)
snp_glm_7 <- slope_melt(snp_abA_7,snp_abB_7)
write_csv(snp_glm_7,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_07.csv")
rm(snp_abA_7)
rm(snp_abB_7)
rm(snp_glm_7)
rm(snp_7)

snp_8 <- snp_all[700001:800000,]
snp_abA_8 <-  abA(snp_8,pop_order_2)
snp_abB_8 <-  abB(snp_8,pop_order_2)
snp_glm_8 <- slope_melt(snp_abA_8,snp_abB_8)
write_csv(snp_glm_8,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_08.csv")
rm(snp_abA_8)
rm(snp_abB_8)
rm(snp_glm_8)
rm(snp_8)

snp_9 <- snp_all[800001:900000,]
snp_abA_9 <-  abA(snp_9,pop_order_2)
snp_abB_9 <-  abB(snp_9,pop_order_2)
snp_glm_9 <- slope_melt(snp_abA_9,snp_abB_9)
write_csv(snp_glm_9,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_09.csv")
rm(snp_abA_9)
rm(snp_abB_9)
rm(snp_glm_9)
rm(snp_9)

snp_10 <- snp_all[900001:1000000,]
snp_abA_10 <-  abA(snp_10,pop_order_2)
snp_abB_10 <-  abB(snp_10,pop_order_2)
snp_glm_10 <- slope_melt(snp_abA_10,snp_abB_10)
write_csv(snp_glm_10,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_10.csv")
rm(snp_abA_10)
rm(snp_abB_10)
rm(snp_glm_10)
rm(snp_10)

###################################################################################

snp_11 <- snp_all[1000001:1100000,]
snp_abA_11 <-  abA(snp_11,pop_order_2)
snp_abB_11 <-  abB(snp_11,pop_order_2)
snp_glm_11 <- slope_melt(snp_abA_11,snp_abB_11)
write_csv(snp_glm_11,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_11.csv")
rm(snp_abA_11)
rm(snp_abB_11)
rm(snp_glm_11)
rm(snp_11)



snp_12 <- snp_all[1100001:1200000,]
snp_abA_12 <-  abA(snp_12,pop_order_2)
snp_abB_12 <-  abB(snp_12,pop_order_2)
snp_glm_12 <- slope_melt(snp_abA_12,snp_abB_12)
write_csv(snp_glm_12,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_12.csv")
rm(snp_abA_12)
rm(snp_abB_12)
rm(snp_glm_12)
rm(snp_12)

snp_13 <- snp_all[1200001:1300000,]
snp_abA_13 <-  abA(snp_13,pop_order_2)
snp_abB_13 <-  abB(snp_13,pop_order_2)
snp_glm_13 <- slope_melt(snp_abA_13,snp_abB_13)
write_csv(snp_glm_13,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_13.csv")
rm(snp_abA_13)
rm(snp_abB_13)
rm(snp_glm_13)
rm(snp_13)

snp_14 <- snp_all[1300001:1400000,]
snp_abA_14 <-  abA(snp_14,pop_order_2)
snp_abB_14 <-  abB(snp_14,pop_order_2)
snp_glm_14 <- slope_melt(snp_abA_14,snp_abB_14)
write_csv(snp_glm_14,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_14.csv")
rm(snp_abA_14)
rm(snp_abB_14)
rm(snp_glm_14)
rm(snp_14)

snp_15 <- snp_all[1400001:1500000,]
snp_abA_15 <-  abA(snp_15,pop_order_2)
snp_abB_15 <-  abB(snp_15,pop_order_2)
snp_glm_15 <- slope_melt(snp_abA_15,snp_abB_15)
write_csv(snp_glm_15,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_15.csv")
rm(snp_abA_15)
rm(snp_abB_15)
rm(snp_glm_15)
rm(snp_15)

snp_16 <- snp_all[1500001:1600000,]
snp_abA_16 <-  abA(snp_16,pop_order_2)
snp_abB_16 <-  abB(snp_16,pop_order_2)
snp_glm_16 <- slope_melt(snp_abA_16,snp_abB_16)
write_csv(snp_glm_16,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_16.csv")
rm(snp_abA_16)
rm(snp_abB_16)
rm(snp_glm_16)
rm(snp_16)

snp_17 <- snp_all[1600001:1700000,]
snp_abA_17 <-  abA(snp_17,pop_order_2)
snp_abB_17 <-  abB(snp_17,pop_order_2)
snp_glm_17 <- slope_melt(snp_abA_17,snp_abB_17)
write_csv(snp_glm_17,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_17.csv")
rm(snp_abA_17)
rm(snp_abB_17)
rm(snp_glm_17)
rm(snp_17)

snp_18 <- snp_all[1700001:1800000,]
snp_abA_18 <-  abA(snp_18,pop_order_2)
snp_abB_18 <-  abB(snp_18,pop_order_2)
snp_glm_18 <- slope_melt(snp_abA_18,snp_abB_18)
write_csv(snp_glm_18,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_18.csv")
rm(snp_abA_18)
rm(snp_abB_18)
rm(snp_glm_18)
rm(snp_18)

snp_19 <- snp_all[1800001:1900000,]
snp_abA_19 <-  abA(snp_19,pop_order_2)
snp_abB_19 <-  abB(snp_19,pop_order_2)
snp_glm_19 <- slope_melt(snp_abA_19,snp_abB_19)
write_csv(snp_glm_19,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_19.csv")
rm(snp_abA_19)
rm(snp_abB_19)
rm(snp_glm_19)
rm(snp_19)

snp_20 <- snp_all[1900001:2000000,]
snp_abA_20 <-  abA(snp_20,pop_order_2)
snp_abB_20 <-  abB(snp_20,pop_order_2)
snp_glm_20 <- slope_melt(snp_abA_20,snp_abB_20)
write_csv(snp_glm_20,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_20.csv")
rm(snp_abA_20)
rm(snp_abB_20)
rm(snp_glm_20)
rm(snp_20)

snp_21 <- snp_all[2000001:2100000,]
snp_abA_21 <-  abA(snp_21,pop_order_2)
snp_abB_21 <-  abB(snp_21,pop_order_2)
snp_glm_21 <- slope_melt(snp_abA_21,snp_abB_21)
write_csv(snp_glm_21,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_21.csv")
rm(snp_abA_21)
rm(snp_abB_21)
rm(snp_glm_21)
rm(snp_21)

snp_22 <- snp_all[2100001:2189349,]
snp_abA_22 <-  abA(snp_22,pop_order_2)
snp_abB_22 <-  abB(snp_22,pop_order_2)
snp_glm_22 <- slope_melt(snp_abA_22,snp_abB_22)
write_csv(snp_glm_22,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_22.csv")
rm(snp_abA_22)
rm(snp_abB_22)
rm(snp_glm_22)
rm(snp_22)









