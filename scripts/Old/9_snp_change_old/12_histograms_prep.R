##################################################################################
## Generate input for strength of selection histographs
## Author Daniel Anstett
## 
## 
## Last Modified May 17, 2022
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)

###################################################################################
#Function for generating stratigied random distribution
###################################################################################

#get ranges using for loop
get_range <- function(df,site){
  env_obs <- data.frame()
  env_slope <- df %>% filter(Site==site)
  for(i in 0:49){
    env_range <- env_slope %>% filter(Slope >= (-2.5+0.1*i) & Slope< (-2.4+0.1*i))
    env_obs[i+1,1] <- dim(env_range)[1]
  }
  return(env_obs)
}



#Histogram Table 
hist_table <- function(df){
  env_obs <- as.data.frame(seq(-2.45, 2.45, by=0.1))
  env1_obs1 <- get_range(df,1)
  env1_obs2 <- get_range(df,2)
  env1_obs3 <- get_range(df,3)
  env1_obs4 <- get_range(df,4)
  env1_obs5 <- get_range(df,5)
  env1_obs6 <- get_range(df,6)
  env1_obs7 <- get_range(df,7)
  env1_obs8 <- get_range(df,8)
  env1_obs9 <- get_range(df,9)
  env1_obs10 <- get_range(df,10)
  env1_obs11 <- get_range(df,11)
  env1_obs12 <- get_range(df,12)
  obs_table <- cbind(env_obs,env1_obs1,env1_obs2,env1_obs3,env1_obs4,env1_obs5,env1_obs6,
                     env1_obs7,env1_obs8,env1_obs9,env1_obs10,env1_obs11,env1_obs12)
  colnames(obs_table) = c("S","p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12")
  
  return(obs_table)
}

#Get Confidence Intervals of permuted slope set
get_ci <- function(df){
  ci_prep <- data.frame()
  for(j in 1:12){
    for(i in 0:49){
      bars <- df %>% filter(Site==j)
      test <- bars %>% filter(Slope >= (-2.5+0.1*i) & Slope< (-2.4+0.1*i))
      count <- test %>% count(Seed_ID)
      ci_prep[1+i,j*2-1] <- as.numeric(quantile(count$n, probs = c(0, 0.95))[1])
      ci_prep[1+i,j*2] <- as.numeric(quantile(count$n, probs = c(0, 0.95))[2])
    }
  }
  colnames(ci_prep) <- c("p1_low","p1_up","p2_low","p2_up","p3_low","p3_up","p4_low","p4_up"
                         ,"p5_low","p5_up","p6_low","p6_up","p7_low","p7_up","p8_low","p8_up"
                         ,"p9_low","p9_up","p10_low","p10_up","p11_low","p11_up","p12_low","p12_up")
  
  return(ci_prep)
}


###################################################################################
#Import observed slopes

#HistPop unique SNP ascorss all env
#Updated for binomial data
obs_env_unique <- read_csv("data/snp_change_data/slope_obs_all_unique.csv") %>% filter(SE<5) 

median_pop <- obs_env_unique %>% group_by(Site) %>% summarise(median = median(Slope, na.rm = TRUE))

write_csv(median_pop, "data/snp_change_data/median_pop.csv")


###################################################################################
#Permuted Data
rand_env_unique <- read_csv("~/Dropbox/z_Documents/aLarge_files/M_gen/rand_slope_histPop_strong_50_50_no_strat.csv") %>%
  filter(Slope <= 2.5 & Slope>= -2.5)  

###############################################################################################
#Run histogram table for obs
###############################################################################################

obs_table_env_unique <- hist_table(obs_env_unique)

###############################################################################################
# Get CI
###############################################################################################
rand_table_env_unique <- get_ci(rand_env_unique)



###############################################################################################
# Merge data frames and export
obs_ci_env_unique <- cbind(obs_table_env_unique,rand_table_env_unique)


#Export
write_csv(obs_ci_env_unique, "data/snp_change_data/obs_ci_env_unique.csv")




