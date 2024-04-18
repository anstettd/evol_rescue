##################################################################################
## Daniel Anstett
## CALCULATE DROUGHT ANOMALIES
## Use modified weather and cliamte data 
## 
## Last Modified January 22, 2020
###################################################################################

###################################################################################
#Load Libraries
library(tidyverse)

###################################################################################
#Water Year Varaibles 
#MAT, MAP, CMD
wna <- read_csv("data/climate_data/demo_climate_wateryear.csv") %>%  #Oct to Setp data
  #Rename variables and take log of MAP
  select(Site=Site, Paper_ID=Paper_ID, MAT.clim=MAT,MAP.clim=MAP,CMD.clim=CMD) %>%
  mutate(log.MAP.clim = log10(MAP.clim))
wna$Site <- as.factor(wna$Site) #Make site a factor

#Import weather data
wna2 <- read_csv("data/climate_data/demo_weather_wateryear.csv") #Import
wna2 <- wna2 %>% #Selects MAT, MAP, CMD,
  select(Site,Paper_ID,Latitude,Longitude,Year,MAT.weath=MAT,MAP.weath=MAP,CMD.weath=CMD) %>%  
  mutate(log.MAP.weath = log10(MAP.weath)) #Take log of MAP
  #separate(ID_Year1, into = c("Site", "Year"), sep = "_") #makes site/year variable
wna2$Site <- as.factor(wna2$Site) ; wna2$Year <- as.numeric(wna2$Year) #define variables

###################################################################################
# join climate and weather, calculate anomaly
wna_all <- left_join(wna2, wna, by=c("Site","Paper_ID")) %>% 
  mutate(CMD.anom = CMD.weath - CMD.clim, 
         MAT.anom = MAT.weath - MAT.clim,
         MAP.anom =  log.MAP.weath - log.MAP.clim)  #reverted to log scale
wna_final <- wna_all %>% select(Site,Paper_ID,Latitude,Longitude,Year,CMD.anom,MAT.anom,MAP.anom)

###################################################################################
#Yearly Variables
#PAS
annual_2010 <- read.csv("data/climate_data/Annual/demography_pop_20_Year_2010Y.csv", header=T) %>% mutate(Year=2010)
annual_2011 <- read.csv("data/climate_data/Annual/demography_pop_20_Year_2011Y.csv", header=T) %>% mutate(Year=2011)
annual_2012 <- read.csv("data/climate_data/Annual/demography_pop_20_Year_2012Y.csv", header=T) %>% mutate(Year=2012)
annual_2013 <- read.csv("data/climate_data/Annual/demography_pop_20_Year_2013Y.csv", header=T) %>% mutate(Year=2013)
annual_2014 <- read.csv("data/climate_data/Annual/demography_pop_20_Year_2014Y.csv", header=T) %>% mutate(Year=2014)
annual_2015 <- read.csv("data/climate_data/Annual/demography_pop_20_Year_2015Y.csv", header=T) %>% mutate(Year=2015)
annual_2016 <- read.csv("data/climate_data/Annual/demography_pop_20_Year_2016Y.csv", header=T) %>% mutate(Year=2016)
annual_2017 <- read.csv("data/climate_data/Annual/demography_pop_20_Year_2017Y.csv", header=T) %>% mutate(Year=2017)
annual_2018 <- read.csv("data/climate_data/Annual/demography_pop_20_Year_2018Y.csv", header=T) %>% mutate(Year=2018)
annual_2019 <- read.csv("data/climate_data/Annual/demography_pop_20_Year_2019Y.csv", header=T) %>% mutate(Year=2019)
annual_8110 <- read.csv("data/climate_data/Annual/demography_pop_20_Normal_1981_2010Y.csv", header=T) %>% mutate(Year=8110)

#Merge an calculate anomaly
annual_8110 <- annual_8110 %>% select(Site,PAS) %>% mutate(PAS.clim=PAS) %>% select(-PAS)
annual <- rbind(annual_2010,
                annual_2011,
                annual_2012,
                annual_2013,
                annual_2014,
                annual_2015,
                annual_2016,
                annual_2017,
                annual_2018,
                annual_2019) %>% select(Site,Paper_ID,Latitude,Longitude,Elevation,PAS,Year)
annual_final <- left_join(annual,annual_8110,by="Site") %>% 
  mutate(PAS.anom = PAS - PAS.clim) %>%
  select(-Elevation,-PAS,-PAS.clim)

###################################################################################
#Seasonal Variables
#Tave_wt, Tave_sm, PPT_wt, PPT_sm
seasonal_2010 <- read.csv("data/climate_data/Seasonal/demography_pop_20_Year_2010S.csv", header=T) %>% mutate(Year=2010)
seasonal_2011 <- read.csv("data/climate_data/Seasonal/demography_pop_20_Year_2011S.csv", header=T) %>% mutate(Year=2011)
seasonal_2012 <- read.csv("data/climate_data/Seasonal/demography_pop_20_Year_2012S.csv", header=T) %>% mutate(Year=2012)
seasonal_2013 <- read.csv("data/climate_data/Seasonal/demography_pop_20_Year_2013S.csv", header=T) %>% mutate(Year=2013)
seasonal_2014 <- read.csv("data/climate_data/Seasonal/demography_pop_20_Year_2014S.csv", header=T) %>% mutate(Year=2014)
seasonal_2015 <- read.csv("data/climate_data/Seasonal/demography_pop_20_Year_2015S.csv", header=T) %>% mutate(Year=2015)
seasonal_2016 <- read.csv("data/climate_data/Seasonal/demography_pop_20_Year_2016S.csv", header=T) %>% mutate(Year=2016)
seasonal_2017 <- read.csv("data/climate_data/Seasonal/demography_pop_20_Year_2017S.csv", header=T) %>% mutate(Year=2017)
seasonal_2018 <- read.csv("data/climate_data/Seasonal/demography_pop_20_Year_2018S.csv", header=T) %>% mutate(Year=2018)
seasonal_2019 <- read.csv("data/climate_data/Seasonal/demography_pop_20_Year_2019S.csv", header=T) %>% mutate(Year=2019)
seasonal_8110 <- read.csv("data/climate_data/Seasonal/demography_pop_20_Normal_1981_2010S.csv", header=T) %>% mutate(Year=8110)

#Merge an calculate anomaly
seasonal_8110 <- seasonal_8110 %>% select(Site,Tave_wt, Tave_sm, PPT_wt, PPT_sm)
colnames(seasonal_8110) <- c("Site","Tave_wt.clim", "Tave_sm.clim", "PPT_wt.clim", "PPT_sm.clim")
seasonal <- rbind(seasonal_2010,
                seasonal_2011,
                seasonal_2012,
                seasonal_2013,
                seasonal_2014,
                seasonal_2015,
                seasonal_2016,
                seasonal_2017,
                seasonal_2018,
                seasonal_2019) %>% select(Site,Paper_ID,Latitude,Longitude,Elevation,Tave_wt, Tave_sm, PPT_wt, PPT_sm,Year)
seasonal_final <- left_join(seasonal,seasonal_8110,by="Site") %>% 
  mutate(Tave_wt.anom = Tave_wt - Tave_wt.clim,
         Tave_sm.anom = Tave_sm - Tave_sm.clim,
         PPT_wt.anom = log10(PPT_wt) - log10(PPT_wt.clim),
         PPT_sm.anom = log10(PPT_sm) - log10(PPT_sm.clim),
         ) %>%
  select(-Elevation,-Tave_wt, -Tave_sm, -PPT_wt, -PPT_sm, -Tave_wt.clim, -Tave_sm.clim, -PPT_wt.clim, -PPT_sm.clim)

#Merge all datatests
anomaly <- left_join(wna_final,annual_final,by=c("Site","Paper_ID","Latitude","Longitude","Year"))
anomaly_final <- left_join(anomaly,seasonal_final,by=c("Site","Paper_ID","Latitude","Longitude","Year"))

#Export Anomaly
write.csv(anomaly_final ,'data/climate_data/climate_anomaly_yearly.csv') #Export file

#Filter for year
anomaly_1215_raw <- anomaly_final %>% filter(Year==2012 | Year==2013 | Year==2014 | Year==2015)
anomaly_1619_raw <- anomaly_final %>% filter(Year==2016 | Year==2017 | Year==2018 | Year==2019)

#Average Anomaly 2012 to 2015
anomaly_1215 <- anomaly_1215_raw %>% group_by(Site,Paper_ID,Latitude,Longitude) %>% 
  summarise(MAT_1215=mean(MAT.anom),
            MAP_1215=mean(MAP.anom),
            PAS_1215=mean(PAS.anom),
            CMD_1215=mean(CMD.anom),
            Tave_wt_1215=mean(Tave_wt.anom),
            Tave_sm_1215=mean(Tave_sm.anom),
            PPT_wt_1215=mean(PPT_wt.anom),
            PPT_sm_1215=mean(PPT_sm.anom)
            )

#Average Anomaly 2016 to 2019
anomaly_1619 <- anomaly_1619_raw %>% group_by(Site,Paper_ID,Latitude,Longitude) %>% 
  summarise(MAT_1619=mean(MAT.anom),
            MAP_1619=mean(MAP.anom),
            PAS_1619=mean(PAS.anom),
            CMD_1619=mean(CMD.anom),
            Tave_wt_1619=mean(Tave_wt.anom),
            Tave_sm_1619=mean(Tave_sm.anom),
            PPT_wt_1619=mean(PPT_wt.anom),
            PPT_sm_1619=mean(PPT_sm.anom)
  )

anom <- left_join(anomaly_1215,anomaly_1619,by=c("Site","Paper_ID","Latitude","Longitude"))


write.csv(anom ,'data/climate_data/climate_anomaly.csv') #Export file
###################################################################################
