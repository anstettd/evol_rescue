##################################################################################
## Associate 9 climate variables
## Author Daniel Anstett
## 
## 
## Last Modified Marc 13, 2023
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)
library(Hmisc)

#Read Data
env_input <- read_csv("data/climate_data/climate_pop.csv")

#Compared the 9 climate variables below:
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

#filter
env <- env_input %>% select(Site_Name:Elevation,MAT,MAP,PAS,EXT,CMD,Tave_wt,Tave_sm,PPT_wt,PPT_sm)
env_m <- as.matrix(env[,6:14])

#Correlate
env_corr <- rcorr(env_m)

env_r <- as.data.frame(env_corr$r)
env_p <- as.data.frame(env_corr$P)

write_csv(env_r,"data/climate_data/rcorr_env_r.csv")
write_csv(env_p,"data/climate_data/rcorr_env_p.csv")


















