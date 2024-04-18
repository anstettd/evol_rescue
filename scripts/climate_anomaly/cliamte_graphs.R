##################################################################################
## Daniel Anstett
## Some timeseries of weather and weather anomaly for each demography site
## 
## Last Modified January 22, 2020
###################################################################################
theme_ci <- function(){ 
  theme_classic() %+replace%    #replace elements we want to change
    theme(axis.text.x = element_text(size = 12, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
          axis.title = element_text(size = 18, face = "bold"), 
          axis.text.y = element_text(size = 12, face = "bold"))
          #strip.background = element_blank(),strip.text.x = element_text(size = 16, face = "bold"))
}
###################################################################################
#Load Libraries
library(tidyverse)
library(RColorBrewer)


###################################################################################
#Import weather data
wna2 <- read_csv("data/climate_data/demo_weather_wateryear.csv")  %>% filter(Paper_ID!=12)
wna2 <- wna2 %>% #Selects MAT, MAP, CMD,
  select(Site,Paper_ID,Latitude,Longitude,Year,MAT.weath=MAT,MAP.weath=MAP,CMD.weath=CMD) %>%  
  mutate(log.MAP.weath = log10(MAP.weath)) #Take log of MAP
#separate(ID_Year1, into = c("Site", "Year"), sep = "_") #makes site/year variable
wna2$Site <- as.factor(wna2$Site) ; wna2$Year <- as.numeric(wna2$Year) #define variables

#Import weather anomaly data
anom <- read_csv("data/climate_data/climate_anomaly_yearly.csv") %>% filter(Paper_ID!=12)

#Make Lat.Site Variable
wna2.site <- wna2 
wna2.site$Latitude <- round(wna2.site$Latitude ,digit=2) 
wna2.site <- wna2.site %>% unite(col=Lat.Site,c("Latitude","Site"),sep="_") %>% select(Lat.Site)
wna2 <- cbind(wna2,wna2.site)

anom.site <- anom 
anom.site$Latitude <- round(anom.site$Latitude ,digit=2) 
anom.site <- anom.site %>% unite(col=Lat.Site,c("Latitude","Site"),sep="_") %>% select(Lat.Site)
anom <- cbind(anom,anom.site)

###################################################################################

#Plot MAP weather
ggplot(wna2, aes(x=Year, y=log.MAP.weath))+
  geom_point()+ geom_line()+ ylab("log Mean Annual Precipitation") +
  scale_x_continuous(breaks=c(2010,2012,2014,2016,2018))+
  geom_hline(yintercept=1, linetype="dotted") +
  facet_wrap(~Lat.Site, scale="free", nrow=5) + theme_classic() #+
  #theme_ci()#+
#theme(strip.background = element_blank(), strip.text.x = element_blank(),
#     legend.title = element_blank())

#Plot MAP anomaly
ggplot(anom, aes(x=Year, y=MAP.anom))+
  geom_point()+ geom_line()+ ylab("Mean Annual Precipitation Anomaly") +
  scale_x_continuous(breaks=c(2010,2012,2014,2016,2018))+
  geom_hline(yintercept=0, linetype="dotted") +
  facet_wrap(~Lat.Site, scale="free", nrow=5) + #theme_classic() #+
theme_ci()#+
#theme(strip.background = element_blank(), strip.text.x = element_blank(),
#     legend.title = element_blank())

ggsave("graphs/Climate/MAP_time.pdf",width=12, height = 8, units = "in")


#Plot PPT_wt anomaly
ggplot(anom, aes(x=Year, y=PPT_wt.anom))+
  geom_point()+ geom_line()+ ylab("Winter Precipitation Anomaly") +
  scale_x_continuous(breaks=c(2010,2012,2014,2016,2018))+
  geom_hline(yintercept=0, linetype="dotted") +
  facet_wrap(~Lat.Site, scale="free", nrow=5) + #theme_classic() #+
  theme_ci()#+
ggsave("graphs/Climate/PPT_wt_time.pdf",width=12, height = 8, units = "in")







