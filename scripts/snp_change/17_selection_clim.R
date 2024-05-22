##################################################################################
## Regress Strength of Selection vs winter precip anomaly
## Deer creek site calculated here
## Author Daniel Anstett
## 
## 
## Last Modified May 20, 2024
###################################################################################
#Library install and import
library(tidyverse)

#Import files
anoms <- read_csv("data/climate_data/climate_anomaly.csv")
slope.summary <- read_csv("data/snp_change_data/mean_median_S.csv") %>% select(Site,Median,Mean)

## CALCULATE DROUGHT ANOMALIES for Deer Creek only
deer_creek_climate <- read.csv("data/climate_data/deer_creek_climate.csv", header=T) %>% 
  select(Year,PPT_wt) 
deer_creek_8110 <- deer_creek_climate$PPT_wt[1]
deer_creek_timeseries <- deer_creek_climate %>% filter(Year>2008) %>% 
  mutate (PPT_wt.anom = log10(PPT_wt) - log10(deer_creek_climate$PPT_wt[1])
  )
deer_creek_1215 <- deer_creek_timeseries %>% filter(Year>2011) %>% filter(Year<2016)

#Join data frames
df1 <- left_join(slope.summary,anoms,by=c("Site"="Paper_ID")) %>% 
  select(Site,Median,Mean,PPT_wt_1215)

#Imput missing data
df1$PPT_wt_1215[10]<-mean(deer_creek_1215$PPT_wt.anom)

