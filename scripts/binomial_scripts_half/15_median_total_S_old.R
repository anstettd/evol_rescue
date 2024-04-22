##################################################################################
## Generate input for strength of selection graphs
## Processing for both observation and permuted CI
## Author Daniel Anstett
## 
## 
## Last Modified Aug 14, 2023
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)

###################################################################################
#HistPop unique SNP ascorss all env
#Updated for binomial data
obs_env_unique <- read_csv("data/binomial_data_half/slope_obs_all_unique.csv") %>% 
                  filter(SE<5.5) %>% mutate(abs_slope = abs(Slope))
obs_env2 <- read_csv("data/binomial_data_half/slope_obs_env2.csv") %>% filter(SE<5.5) %>%
  mutate(abs_slope = abs(Slope))
obs_env9 <- read_csv("data/binomial_data_half/slope_obs_env9.csv") %>% filter(SE<5.5) %>%
  mutate(abs_slope = abs(Slope))


#Get slope median
median_pop <- obs_env_unique %>% group_by(Site) %>% summarise(median = median(Slope, na.rm = TRUE))
median_pop_env2 <- obs_env2 %>% group_by(Site) %>% summarise(median = median(Slope, na.rm = TRUE))
median_pop_env9 <- obs_env9 %>% group_by(Site) %>% summarise(median = median(Slope, na.rm = TRUE))


#Get sum of the absolute value of all slopes
abs_pop <- obs_env_unique %>% group_by(Site) %>% summarise(abs_slope = sum(abs_slope, na.rm = TRUE))
abs_pop_env2 <- obs_env2 %>% group_by(Site) %>% summarise(abs_slope = sum(abs_slope, na.rm = TRUE))
abs_pop_env9 <- obs_env9 %>% group_by(Site) %>% summarise(abs_slope = sum(abs_slope, na.rm = TRUE))

#Get sum of all positive slopes
pos_pop <- obs_env_unique %>% filter(Slope>0) %>% group_by(Site) %>% 
  summarise(pos_slope = sum(Slope, na.rm = TRUE))
pos_pop_env2 <- obs_env2 %>% filter(Slope>0) %>% group_by(Site) %>% 
  summarise(pos_slope = sum(Slope, na.rm = TRUE))
pos_pop_env9 <- obs_env9 %>% filter(Slope>0) %>% group_by(Site) %>% 
  summarise(pos_slope = sum(Slope, na.rm = TRUE))

slope.summary <- cbind(median_pop,abs_pop[,2],pos_pop[,2],
                       median_pop_env2[,2],abs_pop_env2[,2],pos_pop_env2[,2],
                       median_pop_env9[,2],abs_pop_env9[,2],pos_pop_env9[,2])
colnames(slope.summary) <- c("Site","median","abs_slope","pos_slope",
                             "median_env2","abs_slope_env2","pos_slope_env2",
                             "median_env9","abs_slope_env9","pos_slope_env9")
  
write_csv(slope.summary, "data/binomial_data_half/slope_summary.csv")

