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
library(RColorBrewer)


#Import data & Prepare data frame
offset_pop <- read_csv("data/genomic_data/offset_pop_beagle.csv") %>% dplyr::select(Site, Paper_ID) #just to get translation of pop names <--> numbers
offset_pop[20,1] <- "Mill Creek"
offset_pop[20,2] <- 12
demog_recovery <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv")
demog_recovery <- left_join(demog_recovery,offset_pop,by=c("Site"="Site")) %>% 
  rename(Site_Name=Site) %>% 
  filter(Paper_ID!=10) %>% filter(Paper_ID!=12)
anoms <- read_csv("data/climate_data/climate_anomaly.csv") %>% filter(Paper_ID!=12)

demo_pop <- left_join(demog_recovery, anoms, by="Paper_ID")

drought.period <- demo_pop %>% 
  dplyr::select(lambda.slope.decline, lambda.mean.drought, 
                MAT_1215, MAP_1215, PAS_1215, CMD_1215, Tave_wt_1215, 
                Tave_sm_1215, PPT_wt_1215, PPT_sm_1215)
ggpairs(drought.period)


demo_decline_anom <- rcorr(as.matrix(drought.period))


write.csv(demo_decline_anom$r,"data/climate_data/rcorr_demo_decline_anom_r.csv") 
write.csv(demo_decline_anom$P,"data/climate_data/rcorr_demo_decline_anom_p.csv") 


###################################################################################
#Plot decline vs significant env variables

# N-S color gradient
lat_cols=colorRampPalette(brewer.pal(10,"Spectral"))
n.sites <- length(unique(demo_pop$Paper_ID))
color.list <- lat_cols(n.sites)


#Population recovery and climate anomalies
a <- ggplot(demo_pop, aes(x=MAP_1215, y=lambda.mean.drought)) + 
  geom_point(aes(fill=as.factor(round(Latitude.x, 1))),shape=21,size =6)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Mean Lambda During Drought")+
  scale_x_continuous(name="Mean Annual Precipitation Anomaly")+
  #,breaks=c(0.04,0.045,0.05,0.055,0.06))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce hight
  )
ggsave("Graphs/Climate/1_mean_lambda_MAP.pdf",width=8, height = 6, units = "in")

b <- ggplot(demo_pop, aes(x=Tave_sm_1215, y=lambda.mean.drought)) + 
  geom_point(aes(fill=as.factor(round(Latitude.x, 1))),shape=21,size =6)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Mean Lambda During Drought")+
  scale_x_continuous(name="Summer Temperature Anomaly")+
  #,breaks=c(0.04,0.045,0.05,0.055,0.06))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce hight
  )
ggsave("Graphs/Climate/2_mean_lambda_Tave_sm.pdf",width=8, height = 6, units = "in")

c <- ggplot(demo_pop, aes(x=PPT_wt_1215, y=lambda.mean.drought)) + 
  geom_point(aes(fill=as.factor(round(Latitude.x, 1))),shape=21,size =6)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Mean Lambda During Drought")+
  scale_x_continuous(name="Winter Precipitation Anomaly")+
  #,breaks=c(0.04,0.045,0.05,0.055,0.06))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce hight
  )
ggsave("Graphs/Climate/3_mean_lambda_PPT_wt.pdf",width=8, height = 6, units = "in")

d <- ggplot(demo_pop, aes(x=PPT_sm_1215, y=lambda.mean.drought)) + 
  geom_point(aes(fill=as.factor(round(Latitude.x, 1))),shape=21,size =6)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Mean Lambda During Drought")+
  scale_x_continuous(name="Summer Precipitation Anomaly")+
  #,breaks=c(0.04,0.045,0.05,0.055,0.06))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce hight
  )

ggsave("Graphs/Climate/4_mean_lambda_PPT_sm.pdf",width=8, height = 6, units = "in")


