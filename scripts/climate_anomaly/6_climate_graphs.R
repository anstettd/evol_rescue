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


#*******************************************************************************
### 1. Read in climate data
#*******************************************************************************

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



#*******************************************************************************
### 2. Visualize drought anomaly
#*******************************************************************************

#Make populations appear in latitudinal order, filter year range
anom_year<-anom %>% mutate(Site.Lat=paste(Latitude,Paper_ID,sep="_")) %>% filter(Year<2016)
anom_year$Paper_ID<-as.factor(anom_year$Paper_ID) 
anom_year$Paper_ID<-factor(anom_year$Paper_ID,levels=c("55","57","11","9","8","7","6","29","28","27",
                                                       "5","4","3","58","2","17","1","15","14"))
#
Site.label<-c("43.6517_55"="55",
              "43.4137_57"="57",
              "43.37876_11"="11",
              "41.80979_9"="9",
              "41.66546_8"="8",
              "39.74298_7"="7",
              "39.39704_6"="6",
              "37.81878_29"="29",
              "37.81006_28"="28",
              "37.77722_27"="27",
              "37.53822_5"="5",
              "36.68952_4"="4",
              "36.20081_3"="3",
              "36.13756_58"="58",
              "34.28425_2"="2",
              "33.99329_17"="17",
              "32.89913_1"="1",
              "32.75206_15"="15",
              "32.65822_14"="14")


# N-S color gradient
lat_cols=colorRampPalette(brewer.pal(11,"Spectral"))
n.sites <- length(unique(anom_year$Paper_ID))
color.list <- lat_cols(n.sites)

anom_year$Paper_ID <- as.character(anom_year$Paper_ID)
#Use ggplot to generate plot with all required formating
ppt_wt_lat<-ggplot(anom_year, aes(x=Site.Lat, y=PPT_wt.anom, shape=factor(Year), col=factor(Paper_ID)))+ 
  geom_point(aes(fill=Paper_ID), size =6)+
  scale_shape_manual(values=c(48:53))+
  scale_y_continuous(name="PPT Winter Anomaly")+
  xlab("Population")+
  geom_hline(yintercept = 0, color="black", size=0.8)+
  coord_flip()+
  scale_color_manual(values= c("55"="#5E4FA2",
                               "57"="#456EB1",
                               "11"="#378EBA",
                               "9"="#54AEAC",
                               "8"="#75C8A4",
                               "7"="#9BD7A4",
                               "6"="#BEE5A0",
                               "29"="#DFF299",
                               "28"="#F1F9A9",
                               "27"="#FFFFBF",
                               "5"="#FEEDA2",
                               "4"="#FDDA86",
                               "3"="#FDBE6E",
                               "58"="#FB9F5A",
                               "2"="#F67B49",
                               "17"="#E95D47",
                               "1"="#D8434D",
                               "15"="#BC2249",
                               "14"="#9E0142"))+theme_classic()
ppt_wt_lat<-ppt_wt_lat + theme(legend.position = "none",
                               axis.title.x=element_text(size=16,vjust = 0, face="bold",hjust=0.5),
                               #axis.ticks.x = element_blank(),
                               axis.text.x = element_text(size=14, face="bold", angle=0,hjust=0.5),
                               axis.text.y = element_text(size=14,face="bold"),
                               axis.title.y = element_text(size=16,vjust = 1, face="bold",hjust=0.5))+
  scale_x_discrete(labels=Site.label)
ppt_wt_lat

#ggsave("graphs/Climate/PPT_wt_pop.pdf",width=6, height = 8, units = "in")


# plot Year vs. PPT_wt
PPT_wt_plot<-ggplot() + 
  geom_line(stat="smooth",data = anom_year, aes(y = PPT_wt.anom, x = Year, 
        group=Latitude,col=factor(Paper_ID)), alpha=0.75,linewidth = 1.5,se=FALSE) + 
  xlab("Year") + 
  ylab("Winter Precipitation Anomaly") +
  theme_classic()+
  #scale_color_manual(values = lat_cols) +
  scale_color_manual(values= c("55"="#5E4FA2",
                               "57"="#456EB1",
                               "11"="#378EBA",
                               "9"="#54AEAC",
                               "8"="#75C8A4",
                               "7"="#9BD7A4",
                               "6"="#BEE5A0",
                               "29"="#DFF299",
                               "28"="#F1F9A9",
                               "27"="#FFFFBF",
                               "5"="#FEEDA2",
                               "4"="#FDDA86",
                               "3"="#FDBE6E",
                               "58"="#FB9F5A",
                               "2"="#F67B49",
                               "17"="#E95D47",
                               "1"="#D8434D",
                               "15"="#BC2249",
                               "14"="#9E0142"))#+
  #scale_y_discrete(breaks=seq(2010,2016,2))
  #scale_x_discrete(breaks=c(2010,2012,2014,2016)) 
PPT_wt_plot<-PPT_wt_plot + theme(legend.position = "none",
                                 axis.title.x=element_text(size=16,vjust = 0, face="bold",hjust=0.5),
                                 axis.text.x = element_text(size=14, face="bold", angle=0,hjust=0.5),
                                 axis.text.y = element_text(size=14,face="bold"),
                                 axis.title.y = element_text(size=16,vjust = 1, face="bold",hjust=0.5))#+
  #scale_x_discrete(labels=Site.label)
PPT_wt_plot
ggsave("Graphs/Climate/PPT_wt_pop.pdf",width=6, height = 8, units = "in")


#*******************************************************************************
### 3. Make climate timeseries
#*******************************************************************************

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
