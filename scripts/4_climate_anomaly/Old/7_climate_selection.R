##################################################################################
#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Correlate pop decline with strength of selection
#### AUTHOR: Daniel Anstett and Amy Angert
#### DATE LAST MODIFIED: 20241016
###################################################################################

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

#Library install and import
library(tidyverse)

#Import weather data
wna2 <- read_csv("data/climate_data/demo_weather_wateryear.csv") %>%  
  mutate(MAP.log = log10(MAP))  %>% 
  filter(Paper_ID!=12) %>% 
  filter(Paper_ID!=10) %>%
  filter(Year>2011) %>% 
  filter(Year<2015) %>% 
  group_by(Paper_ID) %>%  
  summarise(mean_MAT = mean(MAT),
            mean_MAP = mean(MAP.log),
            mean_CMD = mean(CMD)) %>% ungroup()

seasonal <- read_csv("data/climate_data/climate_seasonal.csv")%>%  
  mutate(mean_PPT_wt_log = log10(PPT_wt))  %>% 
  filter(Paper_ID!=12) %>% 
  filter(Paper_ID!=10) %>%
  filter(Year>2011) %>% 
  filter(Year<2015) %>% 
  group_by(Paper_ID) %>%  
  summarise(mean_PPT_wt = mean(mean_PPT_wt_log),
            #mean_MAP = mean(MAP.log),
            #mean_CMD = mean(CMD)
            ) %>% ungroup()




#Import climate anomalies
#anoms <- read_csv("data/climate_data/climate_anomaly.csv") 

#Import selection data
slope.summary <- read_csv("data/snp_change_data/mean_median_S.csv") %>% select(Site, Median)
colnames(slope.summary) <- c("Paper_ID","Median")

#Join
anom.sel <- left_join(slope.summary,wna2,by="Paper_ID") %>% filter(Paper_ID!=10)
anom.sel <- left_join(anom.sel,seasonal,by="Paper_ID")

#For Color codes
demo_pop <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv") %>% 
  select(Paper_ID,Latitude,Lat.Color)

anom.sel <- left_join(anom.sel,demo_pop,by="Paper_ID") 


###################################################################################
#LM
lm1 <- lm(Median~mean_PPT_wt,data=anom.sel)
summary(lm1)
#Anova(lm3,type="III")

anom.sel$Lat.Color<-as.factor(anom.sel$Lat.Color)
anom.sel$Lat.Color<-factor(anom.sel$Lat.Color,levels=anom.sel$Lat.Color)

ggplot(anom.sel, aes(x=Median, y=mean_PPT_wt)) + 
  geom_smooth(method=lm,color="black",size=1.25,fill="gray71")+
  geom_point(aes(fill=anom.sel$Lat.Color), shape=21, size =6)+
  #geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.7) +
  scale_y_continuous(name="Winter Precip")+#,breaks=c(0.5,1,1.5,2,2.5))+
  scale_x_continuous(name="Median Selection")+#,breaks=c(-0.1,0,0.1,0.2))+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(name = "Latitude (°N)",labels=round(anom.sel$Latitude,1), values=as.character(anom.sel$Lat.Color)) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce height
  )
#ggsave("Graphs/Demo_selection/1_median_slope_recovery_lambda.pdf",width=8, height = 6, units = "in")

