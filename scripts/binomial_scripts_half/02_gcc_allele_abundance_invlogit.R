##################################################################################
## Finding which allele is associated with climate change
## For peak window SNPS (WZA) with BF (BayPass) > 10 and all BF>30
## Done for all alleles positively associated with climate change
## Author Daniel Anstett
## 
## 
## Last Modified Sept 18, 2023
###################################################################################
#Import libraries
library(tidyverse)
library(boot)
###################################################################################
#Functions

#Generate abundance matrix for prop A 
abA <- function(snp_table) {
  snp_prop_A<- snp_table %>% select (chr_snp)
  counter=2
  pop_num=1
  for (i in seq(1,dim(snp_table)[2]-1,2)) {
    
    snpA<-paste("V",i, sep="") #sets up string for snpA for pop num
    P<-paste("P", pop_num, sep="") #sets up paper_ID (population ID)
    tmp<-snp_table %>% select(snpA) #graps all snps per paper_ID
    
    colnames(tmp)<-"A" #renames column headers
    
    snp_prop_A[,counter]<-tmp$A #assign A
    colnames (snp_prop_A)[counter]<-P 
    
    counter<-counter+1
    pop_num<-pop_num+1
  }

  return(snp_prop_A)
}


abB <- function(snp_table) {
  snp_prop_B<- snp_table %>% select (chr_snp)
  counter=2
  pop_num=1
  for (i in seq(1,dim(snp_table)[2]-1,2)) {
    
    snpB<-paste("V",i+1, sep="") #sets up string for snpB for pop num
    P<-paste("P", pop_num, sep="") #sets up paper_ID (population ID)
    tmp<-snp_table %>% select(snpB) #graps all snps per paper_ID
    
    colnames(tmp)<-"B" #renames column headers
    
    snp_prop_B[,counter]<-tmp$B #assign B
    colnames (snp_prop_B)[counter]<-P 
    
    counter<-counter+1
    pop_num<-pop_num+1
  }
  
  return(snp_prop_B)
}



###################################################################################

#Set A as positive association with climate
FATA_p <- function(snp_base,snp_time,climate_table,env_in){
  snp_prop_A_in<-abA(snp_base) # call function
  snp_prop_B_in<-abB(snp_base) # call function
  
  snp_prop_A_in_time<-abA(snp_time)
  snp_prop_B_in_time<-abB(snp_time)
  
  freq.temp <- data.frame()
  for (i in 1:dim(snp_prop_A_in)[1]) { #for each row
    
    env_pop <- climate_table %>% select(Paper_ID,as.character(env_in)) #Select climate data
    snp_prop_A_tmp<-snp_prop_A_in %>% filter(chr_snp==snp_prop_A_in$chr_snp[i]) %>% 
      select(-chr_snp) #filter prop A data 
    snp_prop_B_tmp<-snp_prop_B_in %>% filter(chr_snp==snp_prop_B_in$chr_snp[i]) %>%
      select(-chr_snp) #filter prop B data 
    
    env_pop <- cbind(env_pop,t(snp_prop_A_tmp), t(snp_prop_B_tmp))
    colnames(env_pop)[3] <-"prop_A"
    colnames(env_pop)[4]<-"prop_B"
    
    #df_A <- env_pop  %>% select(prop_A,prop_B)
    #df_B <- env_pop  %>% select(prop_B,prop_A)
    
    df_A <- cbind(env_pop$prop_A,env_pop$prop_B)
    df_B <- cbind(env_pop$prop_B,env_pop$prop_A)
    
    lm.temp_A <- glm(df_A~env_pop[,2],family=binomial) # save glm of climate predicting prop A
    lm.temp_B <- glm(df_B~env_pop[,2],family=binomial) # save glm of climate predicting prop B
    
    #decides if A or B is positively associated 
    if (inv.logit(lm.temp_A$coefficients[2])>0.5){ 
      #print("A")
      tmp_in<-snp_prop_A_in_time %>% filter(chr_snp==snp_prop_A_in$chr_snp[i])
      freq.temp<-rbind(freq.temp,tmp_in)
      
    } else if (inv.logit(lm.temp_B$coefficients[2])>0.5){
      # print("B")
      tmp_in<-snp_prop_B_in_time %>% filter(chr_snp==snp_prop_B_in$chr_snp[i])
      freq.temp<-rbind(freq.temp,tmp_in)
      
    }
  }
  return(freq.temp)
}


FATB_p <- function(snp_base,snp_time,climate_table,env_in){
  snp_prop_A_in<-abA(snp_base) # call function
  snp_prop_B_in<-abB(snp_base) # call function
  
  snp_prop_A_in_time<-abA(snp_time)
  snp_prop_B_in_time<-abB(snp_time)
  
  freq.temp <- data.frame()
  for (i in 1:dim(snp_prop_A_in)[1]) { #for each row
    
    env_pop <- climate_table %>% select(Paper_ID,as.character(env_in)) #Select climate data
    snp_prop_A_tmp<-snp_prop_A_in %>% filter(chr_snp==snp_prop_A_in$chr_snp[i]) %>% 
      select(-chr_snp) #filter prop A data 
    snp_prop_B_tmp<-snp_prop_B_in %>% filter(chr_snp==snp_prop_B_in$chr_snp[i]) %>%
      select(-chr_snp) #filter prop B data 
    
    env_pop <- cbind(env_pop,t(snp_prop_A_tmp), t(snp_prop_B_tmp))
    colnames(env_pop)[3] <-"prop_A"
    colnames(env_pop)[4]<-"prop_B"
    
    #df_A <- env_pop  %>% select(prop_A,prop_B)
    #df_B <- env_pop  %>% select(prop_B,prop_A)
    
    df_A <- cbind(env_pop$prop_A,env_pop$prop_B)
    df_B <- cbind(env_pop$prop_B,env_pop$prop_A)
    
    lm.temp_A <- glm(df_A~env_pop[,2],family=binomial) # save glm of climate predicting prop A
    lm.temp_B <- glm(df_B~env_pop[,2],family=binomial) # save glm of climate predicting prop B
    
    #decides if A or B is positively associated 
    if (inv.logit(lm.temp_A$coefficients[2])>0.5){ 
      #print("A")
      tmp_in<-snp_prop_B_in_time %>% filter(chr_snp==snp_prop_A_in$chr_snp[i])
      freq.temp<-rbind(freq.temp,tmp_in)
      
    } else if (inv.logit(lm.temp_B$coefficients[2])>0.5){
      # print("B")
      tmp_in<-snp_prop_A_in_time %>% filter(chr_snp==snp_prop_B_in$chr_snp[i])
      freq.temp<-rbind(freq.temp,tmp_in)
      
    }
  }
  return(freq.temp)
}





FATA_n <- function(snp_base,snp_time,climate_table,env_in){
  snp_prop_A_in<-abA(snp_base) # call function
  snp_prop_B_in<-abB(snp_base) # call function
  
  snp_prop_A_in_time<-abA(snp_time)
  snp_prop_B_in_time<-abB(snp_time)
  
  freq.temp <- data.frame()
  for (i in 1:dim(snp_prop_A_in)[1]) { #for each row
    
    env_pop <- climate_table %>% select(Paper_ID,as.character(env_in)) #Select climate data
    snp_prop_A_tmp<-snp_prop_A_in %>% filter(chr_snp==snp_prop_A_in$chr_snp[i]) %>% 
      select(-chr_snp) #filter prop A data 
    snp_prop_B_tmp<-snp_prop_B_in %>% filter(chr_snp==snp_prop_B_in$chr_snp[i]) %>%
      select(-chr_snp) #filter prop B data 
    
    env_pop <- cbind(env_pop,t(snp_prop_A_tmp), t(snp_prop_B_tmp))
    colnames(env_pop)[3] <-"prop_A"
    colnames(env_pop)[4]<-"prop_B"
    
    #df_A <- env_pop  %>% select(prop_A,prop_B)
    #df_B <- env_pop  %>% select(prop_B,prop_A)
    
    df_A <- cbind(env_pop$prop_A,env_pop$prop_B)
    df_B <- cbind(env_pop$prop_B,env_pop$prop_A)
    
    lm.temp_A <- glm(df_A~env_pop[,2],family=binomial) # save glm of climate predicting prop A
    lm.temp_B <- glm(df_B~env_pop[,2],family=binomial) # save glm of climate predicting prop B
    
    #decides if A or B is positively associated 
    if (inv.logit(lm.temp_A$coefficients[2])<0.5){ 
      #print("A")
      tmp_in<-snp_prop_A_in_time %>% filter(chr_snp==snp_prop_A_in$chr_snp[i])
      freq.temp<-rbind(freq.temp,tmp_in)
      
    } else if (inv.logit(lm.temp_B$coefficients[2])<0.5){
      # print("B")
      tmp_in<-snp_prop_B_in_time %>% filter(chr_snp==snp_prop_B_in$chr_snp[i])
      freq.temp<-rbind(freq.temp,tmp_in)
      
    }
  }
  return(freq.temp)
}


FATB_n <- function(snp_base,snp_time,climate_table,env_in){
  snp_prop_A_in<-abA(snp_base) # call function
  snp_prop_B_in<-abB(snp_base) # call function
  
  snp_prop_A_in_time<-abA(snp_time)
  snp_prop_B_in_time<-abB(snp_time)
  
  freq.temp <- data.frame()
  for (i in 1:dim(snp_prop_A_in)[1]) { #for each row
    
    env_pop <- climate_table %>% select(Paper_ID,as.character(env_in)) #Select climate data
    snp_prop_A_tmp<-snp_prop_A_in %>% filter(chr_snp==snp_prop_A_in$chr_snp[i]) %>% 
      select(-chr_snp) #filter prop A data 
    snp_prop_B_tmp<-snp_prop_B_in %>% filter(chr_snp==snp_prop_B_in$chr_snp[i]) %>%
      select(-chr_snp) #filter prop B data 
    
    env_pop <- cbind(env_pop,t(snp_prop_A_tmp), t(snp_prop_B_tmp))
    colnames(env_pop)[3] <-"prop_A"
    colnames(env_pop)[4]<-"prop_B"
    
    #df_A <- env_pop  %>% select(prop_A,prop_B)
    #df_B <- env_pop  %>% select(prop_B,prop_A)
    
    df_A <- cbind(env_pop$prop_A,env_pop$prop_B)
    df_B <- cbind(env_pop$prop_B,env_pop$prop_A)
    
    lm.temp_A <- glm(df_A~env_pop[,2],family=binomial) # save glm of climate predicting prop A
    lm.temp_B <- glm(df_B~env_pop[,2],family=binomial) # save glm of climate predicting prop B
    
    #decides if A or B is positively associated 
    if (inv.logit(lm.temp_A$coefficients[2])<0.5){ 
      #print("A")
      tmp_in<-snp_prop_B_in_time %>% filter(chr_snp==snp_prop_A_in$chr_snp[i])
      freq.temp<-rbind(freq.temp,tmp_in)
      
    } else if (inv.logit(lm.temp_B$coefficients[2])<0.5){
      # print("B")
      tmp_in<-snp_prop_A_in_time %>% filter(chr_snp==snp_prop_B_in$chr_snp[i])
      freq.temp<-rbind(freq.temp,tmp_in)
      
    }
  }
  return(freq.temp)
}


########################################################################################################


#Import Baseline Climate 
climate <- read_csv("data/climate_pop.csv")

#Import pop names (site_year names)
pop_order<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")

##SNP set abundances
#Timeseries SNP abundances
snp1_time <- read_csv("data/snp_set_time_env1.csv")
snp2_time <- read_csv("data/snp_set_time_env2.csv")
snp3_time <- read_csv("data/snp_set_time_env3.csv")
snp4_time <- read_csv("data/snp_set_time_env4.csv")
snp5_time <- read_csv("data/snp_set_time_env5.csv")
snp6_time <- read_csv("data/snp_set_time_env6.csv")
snp7_time <- read_csv("data/snp_set_time_env7.csv")
snp8_time <- read_csv("data/snp_set_time_env8.csv")
snp9_time <- read_csv("data/snp_set_time_env9.csv")

#Baseline SNP abundances
env1_base <- read_csv("data/snp_set_base_env1.csv")
env2_base <- read_csv("data/snp_set_base_env2.csv")
env3_base <- read_csv("data/snp_set_base_env3.csv")
env4_base <- read_csv("data/snp_set_base_env4.csv")
env5_base <- read_csv("data/snp_set_base_env5.csv")
env6_base <- read_csv("data/snp_set_base_env6.csv")
env7_base <- read_csv("data/snp_set_base_env7.csv")
env8_base <- read_csv("data/snp_set_base_env8.csv")
env9_base <- read_csv("data/snp_set_base_env9.csv")

###################################################################################

#Ensure baseline has same SNPs as timeseries and remove any not in the timeseries 
env1_base <- env1_base %>% filter (chr_snp %in% as.character(snp1_time$chr_snp))
env2_base <- env2_base %>% filter (chr_snp %in% as.character(snp2_time$chr_snp))
env3_base <- env3_base %>% filter (chr_snp %in% as.character(snp3_time$chr_snp))
env4_base <- env4_base %>% filter (chr_snp %in% as.character(snp4_time$chr_snp))
env5_base <- env5_base %>% filter (chr_snp %in% as.character(snp5_time$chr_snp))
env6_base <- env6_base %>% filter (chr_snp %in% as.character(snp6_time$chr_snp))
env7_base <- env7_base %>% filter (chr_snp %in% as.character(snp7_time$chr_snp))
env8_base <- env8_base %>% filter (chr_snp %in% as.character(snp8_time$chr_snp))
env9_base <- env9_base %>% filter (chr_snp %in% as.character(snp9_time$chr_snp))


###################################################################################
##### Annual #####
# env 1 is MAT = Mean annual temperature (°C)
# env 2 is MAP = Mean annual precipitation (mm)
# env 3 is PAS = Precipitation as snow (mm) between August in previous year and July in current year
# env 4 is EXT = Extreme temperature over 30 years
# env 5 is CMD = Hargreaves climatic moisture deficit (mm)

##### Seasonal #####
# env 6 is Tave_wt = Winter mean temperature (°C)
# env 7 is Tave_sm = Summer mean temperature (°C)
# env 8 is PPT_wt = Winter precipitation (mm)
# env 9 is PPT_sm = Summer precipitation (mm)

## Make table with climate change associated SNP abundances for timeseries using baseline climate
freq_env1_A <- as.data.frame(FATA_p(env1_base,snp1_time,climate,"MAT"))
freq_env2_A <- as.data.frame(FATA_n(env2_base,snp2_time,climate,"MAP"))
freq_env3_A <- as.data.frame(FATA_n(env3_base,snp3_time,climate,"PAS"))
freq_env4_A <- as.data.frame(FATA_p(env4_base,snp4_time,climate,"EXT"))
freq_env5_A <- as.data.frame(FATA_p(env5_base,snp5_time,climate,"CMD"))
freq_env6_A <- as.data.frame(FATA_p(env6_base,snp6_time,climate,"Tave_wt"))
freq_env7_A <- as.data.frame(FATA_p(env7_base,snp7_time,climate,"Tave_sm"))
freq_env8_A <- as.data.frame(FATA_n(env8_base,snp8_time,climate,"PPT_wt"))
freq_env9_A <- as.data.frame(FATA_n(env9_base,snp9_time,climate,"PPT_sm"))

freq_env1_B <- as.data.frame(FATB_p(env1_base,snp1_time,climate,"MAT"))
freq_env2_B <- as.data.frame(FATB_n(env2_base,snp2_time,climate,"MAP"))
freq_env3_B <- as.data.frame(FATB_n(env3_base,snp3_time,climate,"PAS"))
freq_env4_B <- as.data.frame(FATB_p(env4_base,snp4_time,climate,"EXT"))
freq_env5_B <- as.data.frame(FATB_p(env5_base,snp5_time,climate,"CMD"))
freq_env6_B <- as.data.frame(FATB_p(env6_base,snp6_time,climate,"Tave_wt"))
freq_env7_B <- as.data.frame(FATB_p(env7_base,snp7_time,climate,"Tave_sm"))
freq_env8_B <- as.data.frame(FATB_n(env8_base,snp8_time,climate,"PPT_wt"))
freq_env9_B <- as.data.frame(FATB_n(env9_base,snp9_time,climate,"PPT_sm"))

#Make first column a row name
rownames(freq_env1_A)<- as.vector(freq_env1_A$chr_snp)
rownames(freq_env2_A)<- as.vector(freq_env2_A$chr_snp)
rownames(freq_env3_A)<- as.vector(freq_env3_A$chr_snp)
rownames(freq_env4_A)<- as.vector(freq_env4_A$chr_snp)
rownames(freq_env5_A)<- as.vector(freq_env5_A$chr_snp)
rownames(freq_env6_A)<- as.vector(freq_env6_A$chr_snp)
rownames(freq_env7_A)<- as.vector(freq_env7_A$chr_snp)
rownames(freq_env8_A)<- as.vector(freq_env8_A$chr_snp)
rownames(freq_env9_A)<- as.vector(freq_env9_A$chr_snp)

rownames(freq_env1_B)<- as.vector(freq_env1_B$chr_snp)
rownames(freq_env2_B)<- as.vector(freq_env2_B$chr_snp)
rownames(freq_env3_B)<- as.vector(freq_env3_B$chr_snp)
rownames(freq_env4_B)<- as.vector(freq_env4_B$chr_snp)
rownames(freq_env5_B)<- as.vector(freq_env5_B$chr_snp)
rownames(freq_env6_B)<- as.vector(freq_env6_B$chr_snp)
rownames(freq_env7_B)<- as.vector(freq_env7_B$chr_snp)
rownames(freq_env8_B)<- as.vector(freq_env8_B$chr_snp)
rownames(freq_env9_B)<- as.vector(freq_env9_B$chr_snp)

#Remove 1st column
freq_env1_A <- freq_env1_A %>% select(-chr_snp)
freq_env2_A <- freq_env2_A %>% select(-chr_snp)
freq_env3_A <- freq_env3_A %>% select(-chr_snp)
freq_env4_A <- freq_env4_A %>% select(-chr_snp)
freq_env5_A <- freq_env5_A %>% select(-chr_snp)
freq_env6_A <- freq_env6_A %>% select(-chr_snp)
freq_env7_A <- freq_env7_A %>% select(-chr_snp)
freq_env8_A <- freq_env8_A %>% select(-chr_snp)
freq_env9_A <- freq_env9_A %>% select(-chr_snp)

freq_env1_B <- freq_env1_B %>% select(-chr_snp)
freq_env2_B <- freq_env2_B %>% select(-chr_snp)
freq_env3_B <- freq_env3_B %>% select(-chr_snp)
freq_env4_B <- freq_env4_B %>% select(-chr_snp)
freq_env5_B <- freq_env5_B %>% select(-chr_snp)
freq_env6_B <- freq_env6_B %>% select(-chr_snp)
freq_env7_B <- freq_env7_B %>% select(-chr_snp)
freq_env8_B <- freq_env8_B %>% select(-chr_snp)
freq_env9_B <- freq_env9_B %>% select(-chr_snp)

# Add in row names to SNP frequency tables 
colnames(freq_env1_A)<- pop_order[,1] #name each pop/time combination
colnames(freq_env2_A)<- pop_order[,1] 
colnames(freq_env3_A)<- pop_order[,1] 
colnames(freq_env4_A)<- pop_order[,1] 
colnames(freq_env5_A)<- pop_order[,1] 
colnames(freq_env6_A)<- pop_order[,1] 
colnames(freq_env7_A)<- pop_order[,1] 
colnames(freq_env8_A)<- pop_order[,1] 
colnames(freq_env9_A)<- pop_order[,1] 

colnames(freq_env1_B)<- pop_order[,1] #name each pop/time combination
colnames(freq_env2_B)<- pop_order[,1] 
colnames(freq_env3_B)<- pop_order[,1] 
colnames(freq_env4_B)<- pop_order[,1] 
colnames(freq_env5_B)<- pop_order[,1] 
colnames(freq_env6_B)<- pop_order[,1] 
colnames(freq_env7_B)<- pop_order[,1] 
colnames(freq_env8_B)<- pop_order[,1] 
colnames(freq_env9_B)<- pop_order[,1] 

#Transpose and split up site_year
freq_env1_A_T_logit  <- as.data.frame(t(freq_env1_A)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env2_A_T <- as.data.frame(t(freq_env2_A)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env3_A_T <- as.data.frame(t(freq_env3_A)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env4_A_T <- as.data.frame(t(freq_env4_A)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env5_A_T <- as.data.frame(t(freq_env5_A)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env6_A_T <- as.data.frame(t(freq_env6_A)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env7_A_T <- as.data.frame(t(freq_env7_A)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env8_A_T <- as.data.frame(t(freq_env8_A)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env9_A_T <- as.data.frame(t(freq_env9_A)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))

freq_env1_B_T_logit <- as.data.frame(t(freq_env1_B)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env2_B_T <- as.data.frame(t(freq_env2_B)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env3_B_T <- as.data.frame(t(freq_env3_B)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env4_B_T <- as.data.frame(t(freq_env4_B)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env5_B_T <- as.data.frame(t(freq_env5_B)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env6_B_T <- as.data.frame(t(freq_env6_B)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env7_B_T <- as.data.frame(t(freq_env7_B)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env8_B_T <- as.data.frame(t(freq_env8_B)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))
freq_env9_B_T <- as.data.frame(t(freq_env9_B)) %>% rownames_to_column ("site_year") %>% separate(site_year, c("Site","Year"))

# Write out climate change correlated timeseries frequency table 
write_csv(freq_env1_A_T, "data/binomial_data/freqA_env1.csv")
write_csv(freq_env2_A_T, "data/binomial_data/freqA_env2.csv")
write_csv(freq_env3_A_T, "data/binomial_data/freqA_env3.csv")
write_csv(freq_env4_A_T, "data/binomial_data/freqA_env4.csv")
write_csv(freq_env5_A_T, "data/binomial_data/freqA_env5.csv")
write_csv(freq_env6_A_T, "data/binomial_data/freqA_env6.csv")
write_csv(freq_env7_A_T, "data/binomial_data/freqA_env7.csv")
write_csv(freq_env8_A_T, "data/binomial_data/freqA_env8.csv")
write_csv(freq_env9_A_T, "data/binomial_data/freqA_env9.csv")

write_csv(freq_env1_B_T, "data/binomial_data/freqB_env1.csv")
write_csv(freq_env2_B_T, "data/binomial_data/freqB_env2.csv")
write_csv(freq_env3_B_T, "data/binomial_data/freqB_env3.csv")
write_csv(freq_env4_B_T, "data/binomial_data/freqB_env4.csv")
write_csv(freq_env5_B_T, "data/binomial_data/freqB_env5.csv")
write_csv(freq_env6_B_T, "data/binomial_data/freqB_env6.csv")
write_csv(freq_env7_B_T, "data/binomial_data/freqB_env7.csv")
write_csv(freq_env8_B_T, "data/binomial_data/freqB_env8.csv")
write_csv(freq_env9_B_T, "data/binomial_data/freqB_env9.csv")








