##################################################################################
## Make basetime df for observed snps
## Filter high SE from basetime observed snps
## Basetime is the starting conditions of the timeseries
## Observed snps are the climate associated snp set that was selected through WZA & BayPass
## Basetime is required to ensure random sampling during permuation is stratified by 
## starting conditions while not selecting high SE SNPs/slopes
## Author Daniel Anstett
## 
## 
## Last Modified April 26, 2023
###################################################################################
#Import libraries
library(tidyverse)

###################################################################################

#Get basetime frequencies
basetime <- function(freq_env){
  freq_env_1011 <- freq_env %>% filter(Site == 1 & Year==2010 |
                                           Site == 2 & Year==2010 |
                                           Site == 3 & Year==2010 |
                                           Site == 4 & Year==2010 |
                                           Site == 5 & Year==2010 |
                                           Site == 6 & Year==2010 |
                                           Site == 7 & Year==2010 |
                                           Site == 8 & Year==2011 |
                                           Site == 9 & Year==2010 |
                                           Site == 10 & Year==2011 |
                                           Site == 11 & Year==2010 |
                                           Site == 12 & Year==2010)
  return(freq_env_1011)
}


#Insert NA's for high SE SNP/site
insert_NA <- function(basetime_df,high_df){
  basetime_df_sorted<-as.data.frame(basetime_df[order(basetime_df$Site),])
  for (i in 1:dim(high_df)[1]){
    basetime_df_sorted[high_df$Site[i], high_df$snp_ID[i]]<-NA
  }
  return(basetime_df_sorted)
}



###################################################################################
#Import transformed timeseries frequencies from just obs (snp set)

freq_env1 <- read_csv("data/freq_env1.csv")
freq_env2 <- read_csv("data/freq_env2.csv")
freq_env3 <- read_csv("data/freq_env3.csv")
freq_env4 <- read_csv("data/freq_env4.csv")
freq_env5 <- read_csv("data/freq_env5.csv")
freq_env6 <- read_csv("data/freq_env6.csv")
freq_env7 <- read_csv("data/freq_env7.csv")
freq_env8 <- read_csv("data/freq_env8.csv")
freq_env9 <- read_csv("data/freq_env9.csv")

#Get Basetime
basetime_env1 <- basetime(freq_env1)
basetime_env2 <- basetime(freq_env2)
basetime_env3 <- basetime(freq_env3)
basetime_env4 <- basetime(freq_env4)
basetime_env5 <- basetime(freq_env5)
basetime_env6 <- basetime(freq_env6)
basetime_env7 <- basetime(freq_env7)
basetime_env8 <- basetime(freq_env8)
basetime_env9 <- basetime(freq_env9)

#Import slope and SE data
slope_env1 <- read_csv("data/binomial_data/slope_obs_env1.csv")
slope_env2 <- read_csv("data/binomial_data/slope_obs_env2.csv")
slope_env3 <- read_csv("data/binomial_data/slope_obs_env3.csv")
slope_env4 <- read_csv("data/binomial_data/slope_obs_env4.csv")
slope_env5 <- read_csv("data/binomial_data/slope_obs_env5.csv")
slope_env6 <- read_csv("data/binomial_data/slope_obs_env6.csv")
slope_env7 <- read_csv("data/binomial_data/slope_obs_env7.csv")
slope_env8 <- read_csv("data/binomial_data/slope_obs_env8.csv")
slope_env9 <- read_csv("data/binomial_data/slope_obs_env9.csv")

#Filter to remove high SE
env1_high <- slope_env1 %>% filter(SE>5)
env2_high <- slope_env2 %>% filter(SE>5)
env3_high <- slope_env3 %>% filter(SE>5)
env4_high <- slope_env4 %>% filter(SE>5)
env5_high <- slope_env5 %>% filter(SE>5)
env6_high <- slope_env6 %>% filter(SE>5)
env7_high <- slope_env7 %>% filter(SE>5)
env8_high <- slope_env8 %>% filter(SE>5)
env9_high <- slope_env9 %>% filter(SE>5)

#Insert NA
env1_low <- insert_NA(basetime_env1,env1_high)
env2_low <- insert_NA(basetime_env2,env2_high)
env3_low <- insert_NA(basetime_env3,env3_high)
env4_low <- insert_NA(basetime_env4,env4_high)
env5_low <- insert_NA(basetime_env5,env5_high)
env6_low <- insert_NA(basetime_env6,env6_high)
env7_low <- insert_NA(basetime_env7,env7_high)
env8_low <- insert_NA(basetime_env8,env8_high)
env9_low <- insert_NA(basetime_env9,env9_high)





#Export
#write_csv(env1_low, "data/env1_low.csv")
#write_csv(env2_low, "data/env2_low.csv")
#write_csv(env3_low, "data/env3_low.csv")
#write_csv(env4_low, "data/env4_low.csv")
#write_csv(env5_low, "data/env5_low.csv")
#write_csv(env6_low, "data/env6_low.csv")
#write_csv(env7_low, "data/env7_low.csv")
#write_csv(env8_low, "data/env8_low.csv")
#write_csv(env9_low, "data/env9_low.csv")


