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

#get_range
get_range <- function(df,site){
  # env Example
  env_slope <- df %>% filter(Site==site)
  
  #Make dataframe
  env_obs <- data.frame()

  #Get ranges
  env_range1 <- env_slope %>% filter(Slope >= -5.0 & Slope< -4.8)
  env_range2 <- env_slope %>% filter(Slope >= -4.8 & Slope< -4.6)
  env_range3 <- env_slope %>% filter(Slope >= -4.6 & Slope< -4.4)
  env_range4 <- env_slope %>% filter(Slope >= -4.4 & Slope< -4.2)
  env_range5 <- env_slope %>% filter(Slope >= -4.2 & Slope< -4.0)
  env_range6 <- env_slope %>% filter(Slope >= -4.0 & Slope< -3.8)
  env_range7 <- env_slope %>% filter(Slope >= -3.8 & Slope< -3.6)
  env_range8 <- env_slope %>% filter(Slope >= -3.6 & Slope< -3.4)
  env_range9 <- env_slope %>% filter(Slope >= -3.4 & Slope< -3.2)
  env_range10 <- env_slope %>% filter(Slope >= -3.2 & Slope< -3.0)
  env_range11 <- env_slope %>% filter(Slope >= -3.0 & Slope< -2.8)
  env_range12 <- env_slope %>% filter(Slope >= -2.8 & Slope< -2.6)
  
  env_range13 <- env_slope %>% filter(Slope >= -2.6 & Slope< -2.4)
  env_range14 <- env_slope %>% filter(Slope >= -2.4 & Slope < -2.2)
  env_range15 <- env_slope %>% filter(Slope >= -2.2 & Slope < -2.0)
  env_range16 <- env_slope %>% filter(Slope >= -2.0 & Slope< -1.8)
  env_range17 <- env_slope %>% filter(Slope >= -1.8 & Slope < -1.6)
  env_range18 <- env_slope %>% filter(Slope >= -1.6 & Slope < -1.4)
  env_range19 <- env_slope %>% filter(Slope >= -1.4 & Slope < -1.2)
  env_range20 <- env_slope %>% filter(Slope >= -1.2 & Slope < -1.0)
  env_range21 <- env_slope %>% filter(Slope >= -1.0 & Slope < -0.8)
  env_range22 <- env_slope %>% filter(Slope >= -0.8 & Slope < -0.6)
  env_range23 <- env_slope %>% filter(Slope >= -0.6 & Slope < -0.4)
  env_range24 <- env_slope %>% filter(Slope >= -0.4 & Slope < -0.2)
  env_range25 <- env_slope %>% filter(Slope >= -0.2 & Slope < 0.0)
  env_range26 <- env_slope %>% filter(Slope >= 0.0 & Slope < 0.2)
  env_range27 <- env_slope %>% filter(Slope >= 0.2 & Slope < 0.4)
  env_range28 <- env_slope %>% filter(Slope >= 0.3 & Slope < 0.6)
  env_range29 <- env_slope %>% filter(Slope >= 0.4 & Slope < 0.8)
  env_range30 <- env_slope %>% filter(Slope >= 0.8 & Slope < 1.0)
  env_range31 <- env_slope %>% filter(Slope >= 1.0 & Slope < 1.2)
  env_range32 <- env_slope %>% filter(Slope >= 1.2 & Slope < 1.4)
  env_range33 <- env_slope %>% filter(Slope >= 1.4 & Slope < 1.6)
  env_range34 <- env_slope %>% filter(Slope >= 1.6 & Slope < 1.8)
  env_range35 <- env_slope %>% filter(Slope >= 1.8 & Slope < 2.0)
  env_range36 <- env_slope %>% filter(Slope >= 2.0 & Slope < 2.2)
  env_range37 <- env_slope %>% filter(Slope >= 2.2 & Slope < 2.4)
  env_range38 <- env_slope %>% filter(Slope >= 2.4 & Slope < 2.6)
  
  env_range39 <- env_slope %>% filter(Slope >= 2.6 & Slope < 2.8)
  env_range40 <- env_slope %>% filter(Slope >= 2.8 & Slope < 3.0)
  env_range41 <- env_slope %>% filter(Slope >= 3.0 & Slope < 3.2)
  env_range42 <- env_slope %>% filter(Slope >= 3.2 & Slope < 3.4)
  env_range43 <- env_slope %>% filter(Slope >= 3.4 & Slope < 3.6)
  env_range44 <- env_slope %>% filter(Slope >= 3.6 & Slope < 3.8)
  env_range45 <- env_slope %>% filter(Slope >= 3.8 & Slope < 4.0)
  env_range46 <- env_slope %>% filter(Slope >= 4.0 & Slope < 4.2)
  env_range47 <- env_slope %>% filter(Slope >= 4.2 & Slope < 4.4)
  env_range48 <- env_slope %>% filter(Slope >= 4.4 & Slope < 4.6)
  env_range49 <- env_slope %>% filter(Slope >= 4.6 & Slope < 4.8)
  env_range50 <- env_slope %>% filter(Slope >= 4.8 & Slope < 5.0)
  

  
  #get dim and put into dataframe
  env_obs[1,1] <- dim(env_range1)[1]
  env_obs[2,1] <- dim(env_range2)[1]
  env_obs[3,1] <- dim(env_range3)[1]
  env_obs[4,1] <- dim(env_range4)[1]
  env_obs[5,1] <- dim(env_range5)[1]
  env_obs[6,1] <- dim(env_range6)[1]
  env_obs[7,1] <- dim(env_range7)[1]
  env_obs[8,1] <- dim(env_range8)[1]
  env_obs[9,1] <- dim(env_range9)[1]
  env_obs[10,1] <- dim(env_range10)[1]
  env_obs[11,1] <- dim(env_range11)[1]
  env_obs[12,1] <- dim(env_range12)[1]
  env_obs[13,1] <- dim(env_range13)[1]
  env_obs[14,1] <- dim(env_range14)[1]
  env_obs[15,1] <- dim(env_range15)[1]
  env_obs[16,1] <- dim(env_range16)[1]
  env_obs[17,1] <- dim(env_range17)[1]
  env_obs[18,1] <- dim(env_range18)[1]
  env_obs[19,1] <- dim(env_range19)[1]
  env_obs[20,1] <- dim(env_range20)[1]
  env_obs[21,1] <- dim(env_range21)[1]
  env_obs[22,1] <- dim(env_range22)[1]
  env_obs[23,1] <- dim(env_range23)[1]
  env_obs[24,1] <- dim(env_range24)[1]
  env_obs[25,1] <- dim(env_range25)[1]
  env_obs[26,1] <- dim(env_range26)[1]
  env_obs[27,1] <- dim(env_range27)[1]
  env_obs[28,1] <- dim(env_range28)[1]
  env_obs[29,1] <- dim(env_range29)[1]
  env_obs[30,1] <- dim(env_range30)[1]
  env_obs[31,1] <- dim(env_range31)[1]
  env_obs[32,1] <- dim(env_range32)[1]
  env_obs[33,1] <- dim(env_range33)[1]
  env_obs[34,1] <- dim(env_range34)[1]
  env_obs[35,1] <- dim(env_range35)[1]
  env_obs[36,1] <- dim(env_range36)[1]
  env_obs[37,1] <- dim(env_range37)[1]
  env_obs[38,1] <- dim(env_range38)[1]
  env_obs[39,1] <- dim(env_range39)[1]
  env_obs[40,1] <- dim(env_range40)[1]
  env_obs[41,1] <- dim(env_range41)[1]
  env_obs[42,1] <- dim(env_range42)[1]
  env_obs[43,1] <- dim(env_range43)[1]
  env_obs[44,1] <- dim(env_range44)[1]
  env_obs[45,1] <- dim(env_range45)[1]
  env_obs[46,1] <- dim(env_range46)[1]
  env_obs[47,1] <- dim(env_range47)[1]
  env_obs[48,1] <- dim(env_range48)[1]
  env_obs[49,1] <- dim(env_range49)[1]
  env_obs[50,1] <- dim(env_range50)[1]
  
  return(env_obs)
}

#Histogram Table 
hist_table <- function(df){
  env_obs <- as.data.frame(c(-4.9,-4.7,-4.5,-4.3,-4.1,-3.9,-3.7,-3.5,-3.3,-3.1,-2.9,-2.7,
                             -2.5,-2.3,-2.1,-1.9,-1.7,-1.5,-1.3,-1.1,-0.9,-0.7,-0.5,-0.3,-0.1,
                             0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,
                             3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9))
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
      test <- bars %>% filter(Slope >= (-5.0+0.2*i) & Slope< (-4.8+0.2*i))
      count <- test %>% count(Seed_ID)
      ci_prep[1+i,j*2-1] <- as.numeric(quantile(count$n, probs = c(0.025, 0.975))[1])
      ci_prep[1+i,j*2] <- as.numeric(quantile(count$n, probs = c(0.025, 0.975))[2])
    }
  }
    colnames(ci_prep) <- c("p1_low","p1_up","p2_low","p2_up","p3_low","p3_up","p4_low","p4_up"
                             ,"p5_low","p5_up","p6_low","p6_up","p7_low","p7_up","p8_low","p8_up"
                             ,"p9_low","p9_up","p10_low","p10_up","p11_low","p11_up","p12_low","p12_up")
  
  return(ci_prep)
}


###################################################################################
#Import observed slopes
obs_env1 <- read_csv("data/binomial_data/slope_obs_env1.csv") %>% filter(SE<5.5)
obs_env2 <- read_csv("data/binomial_data/slope_obs_env2.csv") %>% filter(SE<5.5)
obs_env3 <- read_csv("data/binomial_data/slope_obs_env3.csv") %>% filter(SE<5.5)
obs_env4 <- read_csv("data/binomial_data/slope_obs_env4.csv") %>% filter(SE<5.5)
obs_env5 <- read_csv("data/binomial_data/slope_obs_env5.csv") %>% filter(SE<5.5)
obs_env6 <- read_csv("data/binomial_data/slope_obs_env6.csv") %>% filter(SE<5.5) %>% filter(Slope>-5)
obs_env7 <- read_csv("data/binomial_data/slope_obs_env7.csv") %>% filter(SE<5.5)
obs_env8 <- read_csv("data/binomial_data/slope_obs_env8.csv") %>% filter(SE<5.5)
obs_env9 <- read_csv("data/binomial_data/slope_obs_env9.csv") %>% filter(SE<5.5)

obs_env <- rbind(obs_env1,
                 obs_env2,
                 obs_env3,
                 obs_env4,
                 obs_env5,
                 obs_env6,
                 obs_env7,
                 obs_env8,
                 obs_env9)

#HistPop unique SNP ascorss all env
#Updated for binomial data
obs_env_unique <- read_csv("data/binomial_data/slope_obs_all_unique.csv") %>% filter(SE<5.5) %>% filter(Slope>-5)



setwd("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files")

#Import large files with 1000 X random slopes
rand_env1 <- read_csv("rand_slope_env1_ab_lowSE.csv") %>% filter(Slope <= 5 & Slope>= -5) 
rand_env2 <- read_csv("rand_slope_env2_ab_lowSE.csv") %>% filter(Slope <= 5 & Slope>= -5) 
rand_env3 <- read_csv("rand_slope_env3_ab_lowSE.csv") %>% filter(Slope <= 5 & Slope>= -5) 
rand_env4 <- read_csv("rand_slope_env4_ab_lowSE.csv") %>% filter(Slope <= 5 & Slope>= -5) 
rand_env5 <- read_csv("rand_slope_env5_ab_lowSE.csv") %>% filter(Slope <= 5 & Slope>= -5) 
rand_env6 <- read_csv("rand_slope_env6_ab_lowSE.csv") %>% filter(Slope <= 5 & Slope>= -5) 
rand_env7 <- read_csv("rand_slope_env7_ab_lowSE.csv") %>% filter(Slope <= 5 & Slope>= -5) 
rand_env8 <- read_csv("rand_slope_env8_ab_lowSE.csv") %>% filter(Slope <= 5 & Slope>= -5) 
rand_env9 <- read_csv("rand_slope_env9_ab_lowSE.csv") %>% filter(Slope <= 5 & Slope>= -5) 

#Updated for bionomial data
rand_env_unique <- read_csv("rand_slope_histPop_lowSE_ab.csv") %>% filter(Slope <= 5 & Slope>= -5)  

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
obs_table_env1 <- hist_table(obs_env1)
obs_table_env2 <- hist_table(obs_env2)
obs_table_env3 <- hist_table(obs_env3)
obs_table_env4 <- hist_table(obs_env4)
obs_table_env5 <- hist_table(obs_env5)
obs_table_env6 <- hist_table(obs_env6)
obs_table_env7 <- hist_table(obs_env7)
obs_table_env8 <- hist_table(obs_env8)
obs_table_env9 <- hist_table(obs_env9)
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
rand_table_env1 <- get_ci(rand_env1)
rand_table_env2 <- get_ci(rand_env2)
rand_table_env3 <- get_ci(rand_env3)
rand_table_env4 <- get_ci(rand_env4)
rand_table_env5 <- get_ci(rand_env5)
rand_table_env6 <- get_ci(rand_env6)
rand_table_env7 <- get_ci(rand_env7)
rand_table_env8 <- get_ci(rand_env8)
rand_table_env9 <- get_ci(rand_env9)
rand_table_env_unique <- get_ci(rand_env_unique)



###############################################################################################
# Merge data frames and export
obs_ci_env1 <- cbind(obs_table_env1,rand_table_env1)
obs_ci_env2 <- cbind(obs_table_env2,rand_table_env2)
obs_ci_env3 <- cbind(obs_table_env3,rand_table_env3)
obs_ci_env4 <- cbind(obs_table_env4,rand_table_env4)
obs_ci_env5 <- cbind(obs_table_env5,rand_table_env5)
obs_ci_env6 <- cbind(obs_table_env6,rand_table_env6)
obs_ci_env7 <- cbind(obs_table_env7,rand_table_env7)
obs_ci_env8 <- cbind(obs_table_env8,rand_table_env8)
obs_ci_env9 <- cbind(obs_table_env9,rand_table_env9)
obs_ci_env_unique <- cbind(obs_table_env_unique,rand_table_env_unique)


#Name each dataframe as a variable
obs_ci_env1$env <- "A MAT"
obs_ci_env2$env <- "B MAP"
obs_ci_env3$env <- "C PAS"
obs_ci_env4$env <- "D EXT"
obs_ci_env5$env <- "E CMD"
obs_ci_env6$env <- "F Tave_wt"
obs_ci_env7$env <- "G Tave_sm"
obs_ci_env8$env <- "H PPT_wt"
obs_ci_env9$env <- "I PPT_sm"

#Bind into single data frame
obs_ci_env <- rbind(obs_ci_env1,
                    obs_ci_env2,
                    obs_ci_env3,
                    obs_ci_env4,
                    obs_ci_env5,
                    obs_ci_env6,
                    obs_ci_env7,
                    obs_ci_env8,
                    obs_ci_env9)


#Export
write_csv(obs_ci_env, "data/binomial_data/obs_ci_env_ab.csv")
#write_csv(obs_ci_env_unique, "data/binomial_data/obs_ci_env_unique.csv")




