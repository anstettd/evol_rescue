##################################################################################
## Daniel Anstett
## CALCULATE CLIMATE FOR 2010-2016 
## Convert monthly temp & precip data October (of the previous year) to September Climate Data
## Average across calcualted years 1980 to 2009 
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
weather_1979 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_1979M.csv", header=T)
weather_1980 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_1980M.csv", header=T)
weather_1981 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_1981M.csv", header=T)
weather_1982 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_1982M.csv", header=T)
weather_1983 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_1983M.csv", header=T)
weather_1984 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_1984M.csv", header=T)
weather_1985 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_1985M.csv", header=T)
weather_1986 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_1986M.csv", header=T)
weather_1987 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_1987M.csv", header=T)
weather_1988 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_1988M.csv", header=T)
weather_1989 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_1989M.csv", header=T)

weather_1990 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_1990M.csv", header=T)
weather_1991 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_1991M.csv", header=T)
weather_1992 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_1992M.csv", header=T)
weather_1993 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_1993M.csv", header=T)
weather_1994 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_1994M.csv", header=T)
weather_1995 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_1995M.csv", header=T)
weather_1996 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_1996M.csv", header=T)
weather_1997 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_1997M.csv", header=T)
weather_1998 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_1998M.csv", header=T)
weather_1999 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_1999M.csv", header=T)

weather_2000 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_2000M.csv", header=T)
weather_2001 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_2001M.csv", header=T)
weather_2002 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_2002M.csv", header=T)
weather_2003 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_2003M.csv", header=T)
weather_2004 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_2004M.csv", header=T)
weather_2005 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_2005M.csv", header=T)
weather_2006 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_2006M.csv", header=T)
weather_2007 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_2007M.csv", header=T)
weather_2008 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_2008M.csv", header=T)
weather_2009 <- read.csv("data/climate_data/Monthly/demography_pop_20_Year_2009M.csv", header=T)

###################################################################################
###################################################################################
# For loop to get Oct to Sept data from all needed years
impact_summary <- data.frame()
for(i in 1980:2009){
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
colnames(impact_summary)[5]<-"hist_year"
###################################################################################

#Write modifed yearly data
#write.csv(impact_summary,'Data/m_year.csv') #Export file

#Write 30-year average
climate_81_10 <- impact_summary %>% group_by(Site,Paper_ID,Latitude,Longitude) %>% summarise_at(c("MAT", "MAP", "CMD"), mean, na.rm=TRUE)
write.csv(climate_81_10,'data/climate_data/demo_climate_wateryear.csv') #Export file
###################################################################################
