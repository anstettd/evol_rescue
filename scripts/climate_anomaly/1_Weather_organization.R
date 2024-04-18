##################################################################################
## Daniel Anstett
## CALCULATE WEATHER FOR 2010-2016 
## Convert monthly temp & precip data October (of the previous year) to September Annual Weather Data
##   
## Last Modified January 22, 2020
###################################################################################

# Monthly climate data obtained from Climate NA 
# http://climatena.ca/

###################################################################################
#Load Libraries
library(tidyverse)

###################################################################################
#Import datasets and add year_actual variable
weather_2009 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_2009M.csv", header=T)
weather_2010 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_2010M.csv", header=T)
weather_2011 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_2011M.csv", header=T)
weather_2012 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_2012M.csv", header=T)
weather_2013 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_2013M.csv", header=T)
weather_2014 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_2014M.csv", header=T)
weather_2015 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_2015M.csv", header=T)
weather_2016 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_2016M.csv", header=T)
weather_2017 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_2017M.csv", header=T)
weather_2018 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_2018M.csv", header=T)
weather_2019 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_2019M.csv", header=T)

impact_all <- weather_2009 %>% select(Site,Paper_ID,Latitude,Longitude)

###################################################################################
###################################################################################
# Get weather data into Oct 1 to Sept 31 format
impact_summary <- data.frame()
for(i in 2010:2019){
  impact<- eval(parse(text=(paste("weather",i-1,sep="_")))) %>% select(Site,Paper_ID,Latitude,Longitude)
  impact<-cbind(impact,c(rep(i,20)))
  #MAT  
  impact_aT <- eval(parse(text=(paste("weather",i-1,sep="_")))) %>% select(Tave10,Tave11,Tave12)
  impact_bT <- eval(parse(text=(paste("weather",i,sep="_"))))%>% select(Tave01,Tave02,Tave03,Tave04,Tave05,
                                                                        Tave06,Tave07,Tave08,Tave09)
  impact_T <- cbind(impact_bT,impact_aT)
  MAT <- rowMeans(impact_T, na.rm = FALSE, dims = 1)
  impact <- cbind(impact,MAT)
  #MAP  
  impact_aT <- eval(parse(text=(paste("weather",i-1,sep="_")))) %>% select(PPT10,PPT11,PPT12)
  impact_bT <- eval(parse(text=(paste("weather",i,sep="_"))))%>% select(PPT01,PPT02,PPT03,PPT04,PPT05,PPT06,PPT07,PPT08,PPT09)
  impact_T <- cbind(impact_bT,impact_aT)
  MAP <- rowSums(impact_T, na.rm = FALSE, dims = 1)
  impact <- cbind(impact,MAP)
  #CMD  
  impact_aT <- eval(parse(text=(paste("weather",i-1,sep="_")))) %>% select(CMD10,CMD11,CMD12)
  impact_bT <- eval(parse(text=(paste("weather",i,sep="_"))))%>% select(CMD01,CMD02,CMD03,CMD04,CMD05,CMD06,CMD07,CMD08,CMD09)
  impact_T <- cbind(impact_bT,impact_aT)
  CMD <- rowSums(impact_T, na.rm = FALSE, dims = 1)
  impact <- cbind(impact,CMD)
  impact_summary <- rbind(impact_summary,impact)
}

colnames(impact_summary)[5] <- "Year"

write.csv(impact_summary,'data/climate_data/demo_weather_wateryear.csv') #Export file
###################################################################################
