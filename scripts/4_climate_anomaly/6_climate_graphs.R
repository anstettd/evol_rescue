##################################################################################
## Daniel Anstett
## Some timeseries of weather and weather anomaly for each demography site
## 
## Last Modified January 22, 2020
###################################################################################
# Remove objects and clear workspace
rm(list = ls(all=TRUE))

#Functions
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


#*******************************************************************************
### 1. Read in climate data
#*******************************************************************************

#Import demography estimates + metadata
demog_means <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv") %>% 
  dplyr::select(Site,Paper_ID,Lat.Color)


#Import weather data
wna2 <- read_csv("data/climate_data/demo_weather_wateryear.csv")  %>% filter(Paper_ID!=12)


wna2 <- wna2 %>% #Selects MAT, MAP, CMD,
  dplyr::select(Site,Paper_ID,Latitude,Longitude,Year,MAT.weath=MAT,MAP.weath=MAP,CMD.weath=CMD) %>%  
  mutate(log.MAP.weath = log10(MAP.weath)) #Take log of MAP
#separate(ID_Year1, into = c("Site", "Year"), sep = "_") #makes site/year variable
wna2$Site <- as.factor(wna2$Site) ; wna2$Year <- as.numeric(wna2$Year) #define variables

#Import weather anomaly data
anom <- read_csv("data/climate_data/climate_anomaly_yearly.csv") %>% filter(Paper_ID!=12)

#Make Lat.Site Variable
wna2.site <- wna2 
wna2.site$Latitude <- round(wna2.site$Latitude ,digit=2) 
wna2.site <- wna2.site %>% unite(col=Lat.Site,c("Latitude","Site"),sep="_") %>% dplyr::select(Lat.Site)
wna2 <- cbind(wna2,wna2.site)
anom.site <- anom 
anom.site$Latitude <- round(anom.site$Latitude ,digit=2) 
anom.site <- anom.site %>% unite(col=Lat.Site,c("Latitude","Site"),sep="_") %>% dplyr::select(Lat.Site)
anom <- cbind(anom,anom.site)
anom <- anom %>% left_join(demog_means, by=c("Site", "Paper_ID"))

#Import weather anomaly data
seasonal <- read_csv("data/climate_data/climate_seasonal.csv") %>% filter(Paper_ID!=12)
seasonal.site <- seasonal 
seasonal.site$Latitude <- round(seasonal.site$Latitude ,digit=2) 
seasonal.site <- seasonal.site %>% unite(col=Lat.Site,c("Latitude","Site"),sep="_") %>% dplyr::select(Lat.Site)
seasonal <- cbind(seasonal,seasonal.site)





#*******************************************************************************
### 2. Visualize climate variables
#*******************************************************************************

#Make populations appear in latitudinal order, filter year range
wna2_year<-wna2 %>% mutate(Site.Lat=paste(Latitude,Paper_ID,sep="_")) %>% filter(Paper_ID!=12) %>%
filter(Year<2016) 
wna2_year<-left_join(wna2_year,demog_means,by="Paper_ID")

wna2_year$Paper_ID<-as.factor(wna2_year$Paper_ID) 
wna2_year$Paper_ID<-factor(wna2_year$Paper_ID,levels=c("55","57","11","9","8","7","6","29","28","27",
                                                       "5","4","3","58","2","17","1","15","14"))


wna2_year$Paper_ID<-as.factor(wna2_year$Paper_ID) 
#wna2_year$Paper_ID<-factor(wna2_year$Paper_ID,levels=wna2_year$Paper_ID)



# plot Year vs. MAP
color_vector <- setNames(wna2_year$Lat.Color, wna2_year$Paper_ID)
MAP_plot<-ggplot() + 
  geom_line(stat="smooth",data = wna2_year, aes(y = MAP.weath, x = Year, 
                                                group=Latitude,col=factor(Paper_ID)), alpha=0.75,linewidth = 1.5,se=FALSE) + 
  xlab("Year") + 
  ylab("Mean Annual Precipitation") +
  theme_classic()+
  #scale_color_manual(values = lat_cols) +
  scale_x_continuous(breaks=c(2010,2012,2014))+
  scale_color_manual(values=color_vector)

MAP_plot<-MAP_plot + theme(legend.position = "none",
                           axis.title.x=element_text(size=24,vjust = 0, face="bold",hjust=0.5),
                           axis.text.x = element_text(size=16, face="bold", angle=0,hjust=0.5),
                           axis.text.y = element_text(size=16,face="bold"),
                           axis.title.y = element_text(size=24,vjust = 1, face="bold",hjust=0.5))#+
MAP_plot
ggsave("Graphs/Climate/MAP_test_pop.pdf",width=6, height = 8, units = "in")




# plot Year vs. Winter Precipitation Anomaly
color_vector <- setNames(anom$Lat.Color, anom$Paper_ID)
PPTwtanom_plot<-ggplot() + 
  geom_line(stat="smooth",data = anom, aes(y = PPT_wt.anom, x = Year,
            group=Latitude,col=factor(Paper_ID)), alpha=0.75,linewidth = 1.5,se=FALSE) + 
  xlab("Year") + 
  ylab("Winter Precipitation Anomaly") +
  theme_classic()+
  #scale_color_manual(values = lat_cols) +
  scale_x_continuous(breaks=c(2010,2012,2014,2016,2018))+
  scale_color_manual(values=color_vector)

PPTwtanom_plot<-PPTwtanom_plot + theme(legend.position = "none",
                           axis.title.x=element_text(size=16,vjust = 0, face="bold",hjust=0.5),
                           axis.text.x = element_text(size=14, face="bold", angle=0,hjust=0.5),
                           axis.text.y = element_text(size=14,face="bold"),
                           axis.title.y = element_text(size=16,vjust = 1, face="bold",hjust=0.5))#+
PPTwtanom_plot
ggsave("Graphs/Climate/PPTwtanom_plot.pdf",width=6, height = 8, units = "in")





