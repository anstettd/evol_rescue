##################################################################################
## Generate input for strength of selection graphs
## Processing for both observation and permuted CI
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
#obs_env1 <- read_csv("data/binomial_data_half/slope_obs_env1.csv") %>% filter(SE<5.5)
#obs_env2 <- read_csv("data/binomial_data_half/slope_obs_env2.csv") %>% filter(SE<5.5)
#obs_env3 <- read_csv("data/binomial_data_half/slope_obs_env3.csv") %>% filter(SE<5.5)
#obs_env4 <- read_csv("data/binomial_data_half/slope_obs_env4.csv") %>% filter(SE<5.5)
#obs_env5 <- read_csv("data/binomial_data_half/slope_obs_env5.csv") %>% filter(SE<5.5)
#obs_env6 <- read_csv("data/binomial_data_half/slope_obs_env6.csv") %>% filter(SE<5.5) %>% filter(Slope>-2.5)
#obs_env7 <- read_csv("data/binomial_data_half/slope_obs_env7.csv") %>% filter(SE<5.5)
#obs_env8 <- read_csv("data/binomial_data_half/slope_obs_env8.csv") %>% filter(SE<5.5)
#obs_env9 <- read_csv("data/binomial_data_half/slope_obs_env9.csv") %>% filter(SE<5.5)

#obs_env <- rbind(obs_env1,
#                 obs_env2,
#                 obs_env3,
#                 obs_env4,
#                 obs_env5,
#                 obs_env6,
#                 obs_env7,
#                 obs_env8,
#                 obs_env9)

#HistPop unique SNP ascorss all env
#Updated for binomial data
obs_env_unique <- read_csv("data/binomial_strong/slope_obs_all_unique.csv") %>% filter(SE<5.5) %>% filter(Slope>-2.5)

median_pop <- obs_env_unique %>% group_by(Site) %>% summarise(median = median(Slope, na.rm = TRUE))
#median_pop_env <- obs_env %>% group_by(Site, ENV) %>% summarise(median = median(Slope, na.rm = TRUE))


write_csv(median_pop, "data/binomial_strong/median_pop.csv")
#write_csv(median_pop_env, "data/binomial_strong/median_pop_env.csv")


###################################################################################
#Permuted Data
#Must removing "%>% mutate(Slope=Slope/2)"  when permutation has been re-run

setwd("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files")

#Import large files with 1000 X random slopes
#rand_env1 <- read_csv("rand_slope_env1_ab_lowSE.csv") %>% filter(Slope <= 2.5 & Slope>= -2.5)
#rand_env2 <- read_csv("rand_slope_env2_ab_lowSE.csv") %>% filter(Slope <= 2.5 & Slope>= -2.5) 
#rand_env3 <- read_csv("rand_slope_env3_ab_lowSE.csv") %>% filter(Slope <= 2.5 & Slope>= -2.5) 
#rand_env4 <- read_csv("rand_slope_env4_ab_lowSE.csv") %>% filter(Slope <= 2.5 & Slope>= -2.5) 
#rand_env5 <- read_csv("rand_slope_env5_ab_lowSE.csv") %>% filter(Slope <= 2.5 & Slope>= -2.5) 
#rand_env6 <- read_csv("rand_slope_env6_ab_lowSE.csv") %>% filter(Slope <= 2.5 & Slope>= -2.5) 
#rand_env7 <- read_csv("rand_slope_env7_ab_lowSE.csv") %>% filter(Slope <= 2.5 & Slope>= -2.5) 
#rand_env8 <- read_csv("rand_slope_env8_ab_lowSE.csv") %>% filter(Slope <= 2.5 & Slope>= -2.5) 
#rand_env9 <- read_csv("rand_slope_env9_ab_lowSE.csv") %>% filter(Slope <= 2.5 & Slope>= -2.5) 

#Updated for bionomial data
rand_env_unique <- read_csv("rand_slope_histPop_strong_minor.csv") %>% filter(Slope <= 2.5 & Slope>= -2.5)  

setwd("~/Dropbox/AM_Workshop/snp_change")

#Setup data frame


###################################################################################
#Set up histogram table obs
#Example for function hist_table
#env1
#env1_obs1 <- get_range(obs_env1,1)
#env1_obs2 <- get_range(obs_env1,2)
#env1_obs3 <- get_range(obs_env1,3)
#env1_obs4 <- get_range(obs_env1,4)
#env1_obs5 <- get_range(obs_env1,5)
#env1_obs6 <- get_range(obs_env1,6)
#env1_obs7 <- get_range(obs_env1,7)
#env1_obs8 <- get_range(obs_env1,8)
#env1_obs9 <- get_range(obs_env1,9)
#env1_obs10 <- get_range(obs_env1,10)
#env1_obs11 <- get_range(obs_env1,11)
#env1_obs12 <- get_range(obs_env1,12)
#env1_obs <- cbind(env_obs,env1_obs1,env1_obs2,env1_obs3,env1_obs4,env1_obs5,env1_obs6,
#                  env1_obs7,env1_obs8,env1_obs9,env1_obs10,env1_obs11,env1_obs12)
#colnames(env1_obs) = c("S","p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12")
###############################################################################################
#Run histogram table for obs
###############################################################################################
#obs_table_env1 <- hist_table(obs_env1)
#obs_table_env2 <- hist_table(obs_env2)
#obs_table_env3 <- hist_table(obs_env3)
#obs_table_env4 <- hist_table(obs_env4)
#obs_table_env5 <- hist_table(obs_env5)
#obs_table_env6 <- hist_table(obs_env6)
#obs_table_env7 <- hist_table(obs_env7)
#obs_table_env8 <- hist_table(obs_env8)
#obs_table_env9 <- hist_table(obs_env9)
obs_table_env_unique <- hist_table(obs_env_unique)


###############################################################################################
#Example of getting ci 

#mat_ci_prep <- data.frame()
#for(j in 1:12){
#  for(i in 0:49){
#    mat_bars <- rand_env1 %>% filter(Site==j)
#    mat_test <- mat_bars %>% filter(Slope >= (-5.0+0.2*i) & Slope< (-4.8+0.2*i))
#    mat_count <- mat_test %>% count(Seed_ID)
#    mat_ci_prep[1+i,j*2-1] <- as.numeric(quantile(mat_count$n, probs = c(0.025, 0.975))[1])
#    mat_ci_prep[1+i,j*2] <- as.numeric(quantile(mat_count$n, probs = c(0.025, 0.975))[2])
#  }
#}

#colnames(mat_ci_prep) <- c("p1_low","p1_up","p2_low","p2_up","p3_low","p3_up","p4_low","p4_up"
#                           ,"p5_low","p5_up","p6_low","p6_up","p7_low","p7_up","p8_low","p8_up"
#                           ,"p9_low","p9_up","p10_low","p10_up","p11_low","p11_up","p12_low","p12_up")

###############################################################################################
# Get CI
###############################################################################################
#Use Get CI function
#rand_table_env1 <- get_ci(rand_env1)
#rand_table_env2 <- get_ci(rand_env2)
#rand_table_env3 <- get_ci(rand_env3)
#rand_table_env4 <- get_ci(rand_env4)
#rand_table_env5 <- get_ci(rand_env5)
#rand_table_env6 <- get_ci(rand_env6)
#rand_table_env7 <- get_ci(rand_env7)
#rand_table_env8 <- get_ci(rand_env8)
#rand_table_env9 <- get_ci(rand_env9)
rand_table_env_unique <- get_ci(rand_env_unique)



###############################################################################################
# Merge data frames and export
#obs_ci_env1 <- cbind(obs_table_env1,rand_table_env1)
#obs_ci_env2 <- cbind(obs_table_env2,rand_table_env2)
#obs_ci_env3 <- cbind(obs_table_env3,rand_table_env3)
#obs_ci_env4 <- cbind(obs_table_env4,rand_table_env4)
#obs_ci_env5 <- cbind(obs_table_env5,rand_table_env5)
#obs_ci_env6 <- cbind(obs_table_env6,rand_table_env6)
#obs_ci_env7 <- cbind(obs_table_env7,rand_table_env7)
#obs_ci_env8 <- cbind(obs_table_env8,rand_table_env8)
#obs_ci_env9 <- cbind(obs_table_env9,rand_table_env9)
obs_ci_env_unique <- cbind(obs_table_env_unique,rand_table_env_unique)


#Name each dataframe as a variable
#obs_ci_env1$env <- "A MAT"
#obs_ci_env2$env <- "B MAP"
#obs_ci_env3$env <- "C PAS"
#obs_ci_env4$env <- "D EXT"
#obs_ci_env5$env <- "E CMD"
#obs_ci_env6$env <- "F Tave_wt"
#obs_ci_env7$env <- "G Tave_sm"
#obs_ci_env8$env <- "H PPT_wt"
#obs_ci_env9$env <- "I PPT_sm"

#Bind into single data frame
#obs_ci_env <- rbind(obs_ci_env1,
#                    obs_ci_env2,
#                    obs_ci_env3,
#                    obs_ci_env4,
#                    obs_ci_env5,
#                    obs_ci_env6,
#                    obs_ci_env7,
#                    obs_ci_env8,
#                    obs_ci_env9)


#Export
#write_csv(obs_ci_env, "data/binomial_strong/obs_ci_env_ab_minor.csv")
write_csv(obs_ci_env_unique, "data/binomial_strong/obs_ci_env_unique_minor.csv")




