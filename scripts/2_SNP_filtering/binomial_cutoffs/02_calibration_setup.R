##################################################################################
## Set up abundances for comparison with baseline climate data
## For peak window SNPS (WZA) with BF (BayPass) > 10 and all BF>30
## Done for all alleles positively associated with climate change
## Author Daniel Anstett
## 
## 
## Last Modified April 22, 2024
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
FATA_p <- function(snp_base,climate_table,env_in){
  snp_prop_A_in<-abA(snp_base) # call function
  snp_prop_B_in<-abB(snp_base) # call function
  
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
    
    df_A <- cbind(env_pop$prop_A,env_pop$prop_B)
    df_B <- cbind(env_pop$prop_B,env_pop$prop_A)
    
    lm.temp_A <- glm(df_A~env_pop[,2],family=binomial) # save glm of climate predicting prop A
    lm.temp_B <- glm(df_B~env_pop[,2],family=binomial) # save glm of climate predicting prop B
    
   if(inv.logit(lm.temp_A$coefficients[2])>0.5){ 
    tmp_in<-env_pop %>% select(Paper_ID, prop_A, prop_B) %>% 
      mutate (ENV=as.character(env_in), chr_snp=snp_prop_A_in$chr_snp[i], SNP_Select="A", 
              Slope=lm.temp_A$coefficients[2],
              Inv_Logit_Coeff=inv.logit(lm.temp_A$coefficients[2]),
              SE=coef(summary(lm.temp_A))[2,2])
    
    colnames(tmp_in)[2]<-"True_SNP_A"
    colnames(tmp_in)[3]<-"True_SNP_B"
    
    freq.temp<-rbind(freq.temp,tmp_in)

  } else if (inv.logit(lm.temp_B$coefficients[2])>0.5){
    tmp_in<-env_pop %>% select(Paper_ID, prop_B, prop_A) %>% 
      mutate (ENV=as.character(env_in), chr_snp=snp_prop_A_in$chr_snp[i], SNP_Select="B",
              Slope=lm.temp_B$coefficients[2],
              Inv_Logit_Coeff=inv.logit(lm.temp_B$coefficients[2]),
              SE=coef(summary(lm.temp_B))[2,2])
    
    colnames(tmp_in)[2]<-"True_SNP_A"
    colnames(tmp_in)[3]<-"True_SNP_B"
    
    freq.temp<-rbind(freq.temp,tmp_in)

  }
  }
  return(freq.temp)
}




FATA_n <- function(snp_base,climate_table,env_in){
  snp_prop_A_in<-abA(snp_base) # call function
  snp_prop_B_in<-abB(snp_base) # call function
  
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
    
    df_A <- cbind(env_pop$prop_A,env_pop$prop_B)
    df_B <- cbind(env_pop$prop_B,env_pop$prop_A)
    
    lm.temp_A <- glm(df_A~env_pop[,2],family=binomial) # save glm of climate predicting prop A
    lm.temp_B <- glm(df_B~env_pop[,2],family=binomial) # save glm of climate predicting prop B
    
    if(inv.logit(lm.temp_A$coefficients[2])<0.5){ 
      tmp_in<-env_pop %>% select(Paper_ID, prop_A, prop_B) %>% 
        mutate (ENV=as.character(env_in), chr_snp=snp_prop_A_in$chr_snp[i], SNP_Select="A", 
                Slope=lm.temp_A$coefficients[2],
                Inv_Logit_Coeff=inv.logit(lm.temp_A$coefficients[2]),
                SE=coef(summary(lm.temp_A))[2,2])
      
      colnames(tmp_in)[2]<-"True_SNP_A"
      colnames(tmp_in)[3]<-"True_SNP_B"
      
      freq.temp<-rbind(freq.temp,tmp_in)

    } else if (inv.logit(lm.temp_B$coefficients[2])<0.5){
      tmp_in<-env_pop %>% select(Paper_ID, prop_B, prop_A) %>% 
        mutate (ENV=as.character(env_in), chr_snp=snp_prop_A_in$chr_snp[i], SNP_Select="B",
                Slope=lm.temp_B$coefficients[2],
                Inv_Logit_Coeff=inv.logit(lm.temp_B$coefficients[2]),
                SE=coef(summary(lm.temp_B))[2,2])
      
      colnames(tmp_in)[2]<-"True_SNP_A"
      colnames(tmp_in)[3]<-"True_SNP_B"
      
      freq.temp<-rbind(freq.temp,tmp_in)
    }
    
  }
  return(freq.temp)
}


#Make bionomial table
binom_table <- function(df){
  abund_table<-data.frame()
  for (i in 1:dim(df)[1]){
    Binomial_A<-c(rep(1, df[i, "True_SNP_A"]), rep(0,df[i, "True_SNP_B"]))
    tmp_df<- as.data.frame(Binomial_A) %>% mutate (Paper_ID=df$Paper_ID[i], 
                                                   chr_snp=df$chr_snp[i], 
                                                   ENV=df$ENV[i])
    abund_table<-rbind(abund_table, tmp_df)
  }
  return(abund_table)
}



########################################################################################################


########################################################################################################


#Import Baseline Climate 
climate <- read_csv("data/snp_change_data/climate_pop.csv")
colnames(climate) <- c("Site_Name","Paper_ID","Latitude","Longitude", "Elevation", "MAT_raw","MAP_raw","PAS_raw",
                       "EXT_raw","CMD_raw","Tave_wt_raw","Tave_sm_raw","PPT_wt_raw","PPT_sm_raw")
climate <- climate %>% mutate(MAT=scale(MAT_raw)[,1],
                              MAP=scale(MAP_raw)[,1],
                              PAS=scale(PAS_raw)[,1],
                              EXT=scale(EXT_raw)[,1],
                              CMD=scale(CMD_raw)[,1],
                              Tave_wt=scale(Tave_wt_raw)[,1],
                              Tave_sm=scale(Tave_sm_raw)[,1],
                              PPT_wt=scale(PPT_wt_raw)[,1],
                              PPT_sm=scale(PPT_sm_raw)[,1]
                              ) %>% select(-MAT_raw,-MAP_raw,-PAS_raw,-EXT_raw,-CMD_raw,
                                           -Tave_wt_raw,-Tave_sm_raw,-PPT_wt_raw,-PPT_sm_raw)


#Import pop names (site_year names)
pop_order<-read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")

pop_order_base<-read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", 
                           header=F, sep="\t")

##SNP set abundances
#Timeseries SNP abundances
snp1_time <- read_csv("data/snp_change_data/snp_set_uncalibrated_time_env1.csv")
snp2_time <- read_csv("data/snp_change_data/snp_set_uncalibrated_time_env2.csv")
snp3_time <- read_csv("data/snp_change_data/snp_set_uncalibrated_time_env3.csv")
snp4_time <- read_csv("data/snp_change_data/snp_set_uncalibrated_time_env4.csv")
snp5_time <- read_csv("data/snp_change_data/snp_set_uncalibrated_time_env5.csv")
snp6_time <- read_csv("data/snp_change_data/snp_set_uncalibrated_time_env6.csv")
snp7_time <- read_csv("data/snp_change_data/snp_set_uncalibrated_time_env7.csv")
snp8_time <- read_csv("data/snp_change_data/snp_set_uncalibrated_time_env8.csv")
snp9_time <- read_csv("data/snp_change_data/snp_set_uncalibrated_time_env9.csv")

#Baseline SNP abundances
env1_base <- read_csv("data/snp_change_data/snp_set_uncalibrated_base_env1.csv")
env2_base <- read_csv("data/snp_change_data/snp_set_uncalibrated_base_env2.csv")
env3_base <- read_csv("data/snp_change_data/snp_set_uncalibrated_base_env3.csv")
env4_base <- read_csv("data/snp_change_data/snp_set_uncalibrated_base_env4.csv")
env5_base <- read_csv("data/snp_change_data/snp_set_uncalibrated_base_env5.csv")
env6_base <- read_csv("data/snp_change_data/snp_set_uncalibrated_base_env6.csv")
env7_base <- read_csv("data/snp_change_data/snp_set_uncalibrated_base_env7.csv")
env8_base <- read_csv("data/snp_change_data/snp_set_uncalibrated_base_env8.csv")
env9_base <- read_csv("data/snp_change_data/snp_set_uncalibrated_base_env9.csv")

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
freq_env1_A <- as.data.frame(FATA_p(env1_base,climate,"MAT"))
freq_env2_A <- as.data.frame(FATA_n(env2_base,climate,"MAP"))
freq_env3_A <- as.data.frame(FATA_n(env3_base,climate,"PAS"))
freq_env4_A <- as.data.frame(FATA_p(env4_base,climate,"EXT"))
freq_env5_A <- as.data.frame(FATA_p(env5_base,climate,"CMD"))
freq_env6_A <- as.data.frame(FATA_p(env6_base,climate,"Tave_wt"))
freq_env7_A <- as.data.frame(FATA_p(env7_base,climate,"Tave_sm"))
freq_env8_A <- as.data.frame(FATA_n(env8_base,climate,"PPT_wt"))
freq_env9_A <- as.data.frame(FATA_n(env9_base,climate,"PPT_sm"))

#Make one data frame
abund_env <- rbind(freq_env1_A,
                   freq_env2_A,
                   freq_env3_A,
                   freq_env4_A,
                   freq_env5_A,
                   freq_env6_A,
                   freq_env7_A,
                   freq_env8_A,
                   freq_env9_A)

#Get slope and SE information for unique SNPs
unique_env <- abund_env %>%
  select(ENV, chr_snp, Slope, Inv_Logit_Coeff, SE) %>%
  unique() %>% filter (SE<5)


#write_csv(unique_env,"data/baseline_SNP_slope.csv")

#Filter for high SE SNPS
abund_env <-abund_env %>% filter (SE<5)

#Filter for low Slope
abund_env_01 <- abund_env %>% filter(abs(Slope)>0.1)
abund_env_02 <- abund_env %>% filter(abs(Slope)>0.2)
abund_env_03 <- abund_env %>% filter(abs(Slope)>0.3)
abund_env_04 <- abund_env %>% filter(abs(Slope)>0.4)
abund_env_05 <- abund_env %>% filter(abs(Slope)>0.5)

##############################################################################################################
# Out of Loop
#abund_table<-data.frame()
#for (i in 1:dim(abund_env)[1]){
#  Binomial_A<-c(rep(1, abund_env[i, "True_SNP_A"]), rep(0,abund_env[i, "True_SNP_B"]))
#  tmp_df<- as.data.frame(Binomial_A) %>% mutate (Paper_ID=abund_env$Paper_ID[i], 
#                                                 chr_snp=abund_env$chr_snp[i], 
#                                                 ENV=abund_env$ENV[i])
#  abund_table<-rbind(abund_table, tmp_df)
#}
##############################################################################################################
#Make long binomial table
abund_table_01 <-binom_table(abund_env_01)
abund_table_02 <-binom_table(abund_env_02)
abund_table_03 <-binom_table(abund_env_03)
abund_table_04 <-binom_table(abund_env_04)
abund_table_05 <-binom_table(abund_env_05)



#Merge with climate data
abund_clim_01 <- left_join(abund_table_01,climate,by="Paper_ID")
abund_clim_02 <- left_join(abund_table_02,climate,by="Paper_ID")
abund_clim_03 <- left_join(abund_table_03,climate,by="Paper_ID")
abund_clim_04 <- left_join(abund_table_04,climate,by="Paper_ID")
abund_clim_05 <- left_join(abund_table_05,climate,by="Paper_ID")



#Write out file
write_csv(abund_clim_01,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/abund_table_baseline_slope_SE_std_01.csv")
write_csv(abund_clim_02,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/abund_table_baseline_slope_SE_std_02.csv")
write_csv(abund_clim_03,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/abund_table_baseline_slope_SE_std_03.csv")
write_csv(abund_clim_04,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/abund_table_baseline_slope_SE_std_04.csv")
write_csv(abund_clim_05,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/abund_table_baseline_slope_SE_std_05.csv")









