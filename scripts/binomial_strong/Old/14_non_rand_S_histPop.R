##################################################################################
## Generate input for strength of selection graphs
## Processing for both observation and permuted CI
## Updated to count only difference between obs and upper CI rather than whole bin
## Author Daniel Anstett
## 
## 
## Last Modified May 17, 2022
###################################################################################
#Import libraries
library(tidyverse)


###################################################################################
#Functions

ns <- function(ci,env_name){
  #Filter for positive slope bin env
  env_pos <- ci %>% filter(S>0)
  env_obs <- env_pos %>% select(p1:p12) 
  env_rand <-env_pos %>% select(p1_up,p2_up,p3_up,p4_up,p5_up,p6_up,p7_up,p8_up,p9_up,p10_up,p11_up,p12_up) 
  instruction_env_P <- env_obs - env_rand 
  instruction_env_P <- cbind (instruction_env_P,env_pos$S)
  colnames(instruction_env_P)[13] <- "S"
  
  #Get non-random slopes
  #env
  slope.cumul.env <-  as.data.frame(c(1,2,3,4,5,6,7,8,9,10,11,12))
  colnames(slope.cumul.env) <- "Site"
  
  for (j in 1:12){
    instruction_env <- instruction_env_P %>% select(eval(paste("p",j, sep="")),S)
    instruction_env <- instruction_env %>% filter(instruction_env[,1]>0)
    Evol <- instruction_env[,1]*instruction_env[,2]
    slope.cumul.env[j,2] <- sum(Evol)
  }
  colnames(slope.cumul.env) <- c("Site",env_name)
  return(slope.cumul.env )
}    


ns_all <- function(ci,env_name){
  #Filter for positive slope bin env
  env_pos <- ci
  env_obs <- env_pos %>% select(p1:p12) 
  env_rand <-env_pos %>% select(p1_up,p2_up,p3_up,p4_up,p5_up,p6_up,p7_up,p8_up,p9_up,p10_up,p11_up,p12_up) 
  instruction_env_P <- env_obs - env_rand 
  instruction_env_P <- cbind (instruction_env_P,env_pos$S)
  colnames(instruction_env_P)[13] <- "S"
  
  #Get non-random slopes
  #env
  slope.cumul.env <-  as.data.frame(c(1,2,3,4,5,6,7,8,9,10,11,12))
  colnames(slope.cumul.env) <- "Site"
  
  for (j in 1:12){
    instruction_env <- instruction_env_P %>% select(eval(paste("p",j, sep="")),S)
    instruction_env <- instruction_env %>% filter(instruction_env[,1]>0)
    Evol <- instruction_env[,1]*instruction_env[,2]
    slope.cumul.env[j,2] <- sum(abs(Evol))
  }
  colnames(slope.cumul.env) <- c("Site",env_name)
  return(slope.cumul.env )
}    



###################################################################################
#Import observed slopes


#Import files that give slopes
env_obs_ci <- read_csv("data/binomial_strong/obs_ci_env_unique.csv")
env_obs_ci[is.na(env_obs_ci)] <- 0

#env_obs_ci_env <- read_csv("data/binomial_strong/obs_ci_env_ab.csv")
#env_obs_ci_env[is.na(env_obs_ci_env)] <- 0

#env_obs_ci_env2 <- env_obs_ci_env %>% filter(env=="B MAP")
#env_obs_ci_env9 <- env_obs_ci_env %>% filter(env=="I PPT_sm")  

#Run ns function to get cumulative score
cumul_pos <-ns(env_obs_ci,"cumul_pos")
cumul_all <-ns_all(env_obs_ci,"cumul_all")

#cumul_pos_env2 <-ns(env_obs_ci_env2,"cumul_pos_env2")
#cumul_all_env2 <-ns_all(env_obs_ci_env2,"cumul_all_env2")

#cumul_pos_env9 <-ns(env_obs_ci_env9,"cumul_pos_env9")
#cumul_all_env9 <-ns_all(env_obs_ci_env9,"cumul_all_env9")


#Join
cumul_grand <- left_join(cumul_pos,cumul_all,by="Site")
#cumul_grand_env2 <- left_join(cumul_pos_env2,cumul_all_env2,by="Site")
#cumul_grand_env9 <- left_join(cumul_pos_env9,cumul_all_env9,by="Site")
#cumul_grand <- cbind(cumul_grand,cumul_grand_env2[2:3],cumul_grand_env9[2:3])


##Integrate with timeseries
timeseries <- read_csv("data/offset_pop_timeseries_beagle.csv")
time_cumul <-left_join(timeseries,cumul_grand,by=c("Paper_ID"="Site"))
time_cumul$Region <- c("South",
                       "South",
                       "Center",
                       "Center",
                       "Center",
                       "Center",
                       "Center",
                       "North",
                       "North",
                       "North",
                       "North",
                       "South")

#Export
write_csv(cumul_pos,"data/binomial_strong/cumul_pos.csv")
write_csv(cumul_all,"data/binomial_strong/cumul_all.csv")
write_csv(time_cumul,"data/binomial_strong/time_cumul_beagle.csv")



