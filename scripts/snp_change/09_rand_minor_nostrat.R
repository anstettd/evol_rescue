##################################################################################
## Generate random distribution of random non-climate associated slopes
## For unique SNP slopes across all env 
## Author Daniel Anstett
## 
## 
## Last Modified Dec 12, 2023
###################################################################################
#Import libraries
library(tidyverse)


###################################################################################



#################################################################################################
# Import non-climate associated slopes with low SE
swiss_glm <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_50_50.csv")

#Import obs (snp set) slopes
obs_env_unique <- read_csv("data/binomial_strong/slope_obs_all_unique.csv") %>% 
  filter(SE<5.5) %>% mutate(abs_slope = abs(Slope))


############################################################################################

rand_master <- data.frame()
for(i in 1:12){
  for (seed_num in 1:1000){
    
    set.seed(seed_num)
    rand_slope_1 <- data.frame()
    swiss_glm_1 <- swiss_glm %>% filter(Site==i)
    obs_env_unique_1 <- obs_env_unique %>% filter(Site==i)
    
    rand_slope_1 <- swiss_glm_1[sample(nrow(swiss_glm_1), dim(obs_env_unique_1)[1]), ] %>% mutate(Seed_ID=seed_num)
    rand_master <- rbind(rand_master,rand_slope_1)
  }
  print(i)
}



write_csv(rand_master, "~/Dropbox/AM_Workshop/Large_files/rand_slope_histPop_strong_50_50_no_strat.csv")










