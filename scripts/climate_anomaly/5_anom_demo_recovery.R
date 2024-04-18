##################################################################################
## Correlate pop decline with climate
## Author Daniel Anstett
## 
## 
## Last Modified Aug 14, 2023
###################################################################################
#Library install and import
library(tidyverse)
library(GGally)
library(Hmisc)


#Import data & Prepare data frame
offset_pop <- read_csv("data/genomic_data/offset_pop_beagle.csv") %>% dplyr::select(Site, Paper_ID) #just to get translation of pop names <--> numbers
offset_pop[20,1] <- "Mill Creek"
offset_pop[20,2] <- 12
demog_recovery <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv")
#cumul <- read_csv("~/Documents/Git repos/snp_change/data/binomial_data_half/time_cumul_beagle.csv")
cumul <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/snp_change/data/time_cumul_beagle.csv")
demog_recovery <- left_join(demog_recovery,offset_pop,by=c("Site"="Site")) %>% rename(Site_Name=Site)
demo_pop <- left_join(cumul,demog_recovery,by="Paper_ID") %>% filter(Paper_ID!=10) %>% filter(Paper_ID!=12)
anoms <- read_csv("data/climate_data/climate_anomaly.csv")
demo_pop <- left_join(demo_pop, anoms)


recovery.period <- demo_pop %>% 
  dplyr::select(lambda.slope.recovery, lambda.mean.recovery, 
                MAT_1619, MAP_1619, PAS_1619, CMD_1619, Tave_wt_1619, 
                Tave_sm_1619, PPT_wt_1619, PPT_sm_1619)
ggpairs(recovery.period)

demo_recovery_anom <- rcorr(as.matrix(recovery.period))



write.csv(demo_recovery_anom$r,"data/climate_data/rcorr_demo_recovery_anom_r.csv") 
write.csv(demo_recovery_anom$P,"data/climate_data/rcorr_demo_recovery_anom_p.csv") 
