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
library(car)

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
  select(Site, Latitude,Median,Mean,PPT_wt_1215)

#Imput missing data
df1$Latitude[10]<-42.27411
df1$PPT_wt_1215[10]<-mean(deer_creek_1215$PPT_wt.anom)

#No Mill Creek
df1.mod <- df1 %>% filter(Site!=12)


#stats
lm_mean <- lm(PPT_wt_1215~Mean,data=df1)
summary(lm_mean)
Anova(lm_mean,type="III")

lm_mod <- lm(PPT_wt_1215~Mean,data=df1.mod)
summary(lm_mod)
Anova(lm_mod,type="III")

###################################################################################
#plot(df1$PPT_wt_1215,df1$Median)

# N-S color gradient
lat_cols=c("#DC494C","#F88D51","#FDD380","#FEEB9E","#FFFFBF","#D7EF9B","#B2E0A2",
           "#88CFA4","#5FBAA8","#3F96B7","#4075B4", "#5E4FA2")


#Median slope vs. ppt_w
ggplot(df1, aes(x=PPT_wt_1215, y=Mean)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =6)+
  geom_smooth(method=lm,color="black", lty="dashed", se=FALSE)+
  scale_y_continuous(name="Strength of Selection")+ #,breaks=c(0.5,1,1.5,2,2.5))+
  scale_x_continuous(name="Winter Precipitation Anomaly")+ #,breaks=c(-0.1,0,0.1,0.2))+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=lat_cols) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce height
  )

ggsave("Graphs/Selection_climate/1_mean_ppt_wt.pdf",width=8, height = 6, units = "in")



#Median slope vs. ppt_w no Mill Creek
ggplot(df1.mod, aes(x=PPT_wt_1215, y=Mean)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =6)+
  geom_smooth(method=lm,color="black", lty="dashed", se=FALSE)+
  scale_y_continuous(name="Strength of Selection")+ #,breaks=c(0.5,1,1.5,2,2.5))+
  scale_x_continuous(name="Winter Precipitation Anomaly")+ #,breaks=c(-0.1,0,0.1,0.2))+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=lat_cols) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce height
  )

ggsave("Graphs/Selection_climate/2_mean_ppt_wt.pdf",width=8, height = 6, units = "in")




