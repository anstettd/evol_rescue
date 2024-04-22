##################################################################################
## Calculate glm slope and SE for every BF>0 snp
## Computationally intensive
## Split into 100k chucks to make managable
## Author Daniel Anstett
## 
## 
## Last Modified April 26, 2023
###################################################################################
#Import libraries
library(tidyverse)

###################################################################################
#Functions

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
  
  colnames(snp_prop_A)<- pop_ID[,1] #name each pop/time combination
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
  
  colnames(snp_prop_A)<- pop_ID[,1] #name each pop/time combination
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
pop_order_2 <- data.frame()
pop_order_2 [1,1] <- "chr_shp"
pop_order_2 <- rbind(pop_order_2,pop_order)

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
#Test GLM for Swiss

#swiss_1 <- snp_swiss[1:10,] #Split data set into 100k ch
#swiss_abA <-  abA(swiss_1,pop_order_2)
#swiss_abB <-  abB(swiss_1,pop_order_2)
#swiss_propA_1 <- prop_A(swiss_1,pop_order_2) # Get frequency table

#req_1.melted <- reshape2::melt(swiss_abA, id.vars = c("Site", "Year"))
#freq_2.melted <- reshape2::melt(swiss_abB, id.vars = c("Site", "Year"))
#colnames(freq_1.melted)[3] <- "snp_ID"
#colnames(freq_1.melted)[4] <- "freq"
#colnames(freq_2.melted)[3] <- "snp_ID"
#colnames(freq_2.melted)[4] <- "freq"

#swiss_glm_1 <- slope_melt(swiss_propA_1) #Run glm function

#################################################################################################

#GLM for Swiss

swiss_1 <- snp_swiss[1:100000,] #Split data set into 100k ch
#swiss_1 <- snp_swiss[1:100,] #Split data set into 100k ch
swiss_abA_1 <-  abA(swiss_1,pop_order_2)
swiss_abB_1 <-  abB(swiss_1,pop_order_2)
swiss_glm_1 <- slope_melt(swiss_abA_1,swiss_abB_1) #Run glm function
write_csv(swiss_glm_1,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_1.csv")

rm(swiss_1)
rm(swiss_abA_1)
rm(swiss_abB_1)
rm(swiss_glm_1)


swiss_2 <- snp_swiss[100001:200000,] #Split data set into 100k ch
swiss_abA_2 <-  abA(swiss_2,pop_order_2)
swiss_abB_2 <-  abB(swiss_2,pop_order_2)
swiss_glm_2 <- slope_melt(swiss_abA_2,swiss_abB_2) #Run glm function
write_csv(swiss_glm_2,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_2.csv")


rm(swiss_2)
rm(swiss_abA_2)
rm(swiss_abB_2)
rm(swiss_glm_2)


swiss_3 <- snp_swiss[200001:300000,] #Split data set into 100k ch
swiss_abA_3 <-  abA(swiss_3,pop_order_2)
swiss_abB_3 <-  abB(swiss_3,pop_order_2)
swiss_glm_3 <- slope_melt(swiss_abA_3,swiss_abB_3) #Run glm function
write_csv(swiss_glm_3,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_3.csv")


rm(swiss_3)
rm(swiss_abA_3)
rm(swiss_abB_3)
rm(swiss_glm_3)

swiss_4 <- snp_swiss[300001:400000,] #Split data set into 100k ch
swiss_abA_4 <-  abA(swiss_4,pop_order_2)
swiss_abB_4 <-  abB(swiss_4,pop_order_2)
swiss_glm_4 <- slope_melt(swiss_abA_4,swiss_abB_4) #Run glm function
write_csv(swiss_glm_4,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_4.csv")


rm(swiss_4)
rm(swiss_abA_4)
rm(swiss_abB_4)
rm(swiss_glm_4)

swiss_5 <- snp_swiss[400001:500000,] #Split data set into 100k ch
swiss_abA_5 <-  abA(swiss_5,pop_order_2)
swiss_abB_5 <-  abB(swiss_5,pop_order_2)
swiss_glm_5 <- slope_melt(swiss_abA_5,swiss_abB_5) #Run glm function
write_csv(swiss_glm_5,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_5.csv")


rm(swiss_5)
rm(swiss_abA_5)
rm(swiss_abB_5)
rm(swiss_glm_5)

swiss_6 <- snp_swiss[500001:600000,] #Split data set into 100k ch
swiss_abA_6 <-  abA(swiss_6,pop_order_2)
swiss_abB_6 <-  abB(swiss_6,pop_order_2)
swiss_glm_6 <- slope_melt(swiss_abA_6,swiss_abB_6) #Run glm function
write_csv(swiss_glm_6,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_6.csv")


rm(swiss_6)
rm(swiss_abA_6)
rm(swiss_abB_6)
rm(swiss_glm_6)

swiss_7 <- snp_swiss[600001:700000,] #Split data set into 100k ch
swiss_abA_7 <-  abA(swiss_7,pop_order_2)
swiss_abB_7 <-  abB(swiss_7,pop_order_2)
swiss_glm_7 <- slope_melt(swiss_abA_7,swiss_abB_7) #Run glm function
write_csv(swiss_glm_7,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_7.csv")


rm(swiss_7)
rm(swiss_abA_7)
rm(swiss_abB_7)
rm(swiss_glm_7)

swiss_8 <- snp_swiss[700001:800000,] #Split data set into 100k ch
swiss_abA_8 <-  abA(swiss_8,pop_order_2)
swiss_abB_8 <-  abB(swiss_8,pop_order_2)
swiss_glm_8 <- slope_melt(swiss_abA_8,swiss_abB_8) #Run glm function
write_csv(swiss_glm_8,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_8.csv")


rm(swiss_8)
rm(swiss_abA_8)
rm(swiss_abB_8)
rm(swiss_glm_8)

swiss_9 <- snp_swiss[800001:900000,] #Split data set into 100k ch
swiss_abA_9 <-  abA(swiss_9,pop_order_2)
swiss_abB_9 <-  abB(swiss_9,pop_order_2)
swiss_glm_9 <- slope_melt(swiss_abA_9,swiss_abB_9) #Run glm function
write_csv(swiss_glm_9,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_9.csv")


rm(swiss_9)
rm(swiss_abA_9)
rm(swiss_abB_9)
rm(swiss_glm_9)

swiss_10 <- snp_swiss[900001:1000000,] #Split data set into 100k ch
swiss_abA_10 <-  abA(swiss_10,pop_order_2)
swiss_abB_10 <-  abB(swiss_10,pop_order_2)
swiss_glm_10 <- slope_melt(swiss_abA_10,swiss_abB_10) #Run glm function
write_csv(swiss_glm_10,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_10.csv")


rm(swiss_10)
rm(swiss_abA_10)
rm(swiss_abB_10)
rm(swiss_glm_10)


swiss_11 <- snp_swiss[1000001:1100000,] #Split data set into 100k ch
swiss_abA_11 <-  abA(swiss_11,pop_order_2)
swiss_abB_11 <-  abB(swiss_11,pop_order_2)
swiss_glm_11 <- slope_melt(swiss_abA_11,swiss_abB_11) #Run glm function
write_csv(swiss_glm_11,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_11.csv")


rm(swiss_11)
rm(swiss_abA_11)
rm(swiss_abB_11)
rm(swiss_glm_11)

swiss_12 <- snp_swiss[1100001:1200000,] #Split data set into 100k ch
swiss_abA_12 <-  abA(swiss_12,pop_order_2)
swiss_abB_12 <-  abB(swiss_12,pop_order_2)
swiss_glm_12 <- slope_melt(swiss_abA_12,swiss_abB_12) #Run glm function
write_csv(swiss_glm_12,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_12.csv")


rm(swiss_12)
rm(swiss_abA_12)
rm(swiss_abB_12)
rm(swiss_glm_12)

swiss_13 <- snp_swiss[1200001:1300000,] #Split data set into 100k ch
swiss_abA_13 <-  abA(swiss_13,pop_order_2)
swiss_abB_13 <-  abB(swiss_13,pop_order_2)
swiss_glm_13 <- slope_melt(swiss_abA_13,swiss_abB_13) #Run glm function
write_csv(swiss_glm_13,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_13.csv")


rm(swiss_13)
rm(swiss_abA_13)
rm(swiss_abB_13)
rm(swiss_glm_13)

swiss_14 <- snp_swiss[1300001:1400000,] #Split data set into 100k ch
swiss_abA_14 <-  abA(swiss_14,pop_order_2)
swiss_abB_14 <-  abB(swiss_14,pop_order_2)
swiss_glm_14 <- slope_melt(swiss_abA_14,swiss_abB_14) #Run glm function
write_csv(swiss_glm_14,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_14.csv")


rm(swiss_14)
rm(swiss_abA_14)
rm(swiss_abB_14)
rm(swiss_glm_14)

swiss_15 <- snp_swiss[1400001:1500000,] #Split data set into 100k ch
swiss_abA_15 <-  abA(swiss_15,pop_order_2)
swiss_abB_15 <-  abB(swiss_15,pop_order_2)
swiss_glm_15 <- slope_melt(swiss_abA_15,swiss_abB_15) #Run glm function
write_csv(swiss_glm_15,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_15.csv")


rm(swiss_15)
rm(swiss_abA_15)
rm(swiss_abB_15)
rm(swiss_glm_15)

swiss_16 <- snp_swiss[1500001:1600000,] #Split data set into 100k ch
swiss_abA_16 <-  abA(swiss_16,pop_order_2)
swiss_abB_16 <-  abB(swiss_16,pop_order_2)
swiss_glm_16 <- slope_melt(swiss_abA_16,swiss_abB_16) #Run glm function
write_csv(swiss_glm_16,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_16.csv")


rm(swiss_16)
rm(swiss_abA_16)
rm(swiss_abB_16)
rm(swiss_glm_16)

swiss_17 <- snp_swiss[1600001:1700000,] #Split data set into 100k ch
swiss_abA_17 <-  abA(swiss_17,pop_order_2)
swiss_abB_17 <-  abB(swiss_17,pop_order_2)
swiss_glm_17 <- slope_melt(swiss_abA_17,swiss_abB_17) #Run glm function
write_csv(swiss_glm_17,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_17.csv")


rm(swiss_17)
rm(swiss_abA_17)
rm(swiss_abB_17)
rm(swiss_glm_17)

swiss_18 <- snp_swiss[1700001:1800000,] #Split data set into 100k ch
swiss_abA_18 <-  abA(swiss_18,pop_order_2)
swiss_abB_18 <-  abB(swiss_18,pop_order_2)
swiss_glm_18 <- slope_melt(swiss_abA_18,swiss_abB_18) #Run glm function
write_csv(swiss_glm_18,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_18.csv")


rm(swiss_18)
rm(swiss_abA_18)
rm(swiss_abB_18)
rm(swiss_glm_18)

swiss_19 <- snp_swiss[1800001:1900000,] #Split data set into 100k ch
swiss_abA_19 <-  abA(swiss_19,pop_order_2)
swiss_abB_19 <-  abB(swiss_19,pop_order_2)
swiss_glm_19 <- slope_melt(swiss_abA_19,swiss_abB_19) #Run glm function
write_csv(swiss_glm_19,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_19.csv")


rm(swiss_19)
rm(swiss_abA_19)
rm(swiss_abB_19)
rm(swiss_glm_19)

swiss_20 <- snp_swiss[1900001:1982192,] #Split data set into 100k ch
swiss_abA_20 <-  abA(swiss_20,pop_order_2)
swiss_abB_20 <-  abB(swiss_20,pop_order_2)
swiss_glm_20 <- slope_melt(swiss_abA_20,swiss_abB_20) #Run glm function
write_csv(swiss_glm_20,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_20.csv")


rm(swiss_20)
rm(swiss_abA_20)
rm(swiss_abB_20)
rm(swiss_glm_20)


