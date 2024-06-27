##################################################################################
#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Correlate pop decline with climate
#### AUTHOR: Daniel Anstett and Amy Angert
#### DATE LAST MODIFIED: 20240627
###################################################################################

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

#Library install and import
library(tidyverse)
library(GGally)
library(Hmisc)
library(RColorBrewer)

#Import population metadata
pop_meta <- read_csv("data/genomic_data/pop_meta_data.csv") 

#Import demography estimates
demog_means <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv") %>% select(-Region, -RegionRank)
demog_means <- left_join(demog_means,pop_meta)

#Import climate anomalies
anoms <- read_csv("data/climate_data/climate_anomaly.csv") 

#Join and filter
demo_pop <- left_join(demog_means, anoms) %>% 
  filter(Site!="Mill Creek") #remove Mill Creek because of unreliable demography data during decline period

drought.period <- demo_pop %>% 
  dplyr::select(mean.lambda.drought, 
                MAT_1215, 
                MAP_1215, 
                PAS_1215, 
                CMD_1215, 
                Tave_wt_1215, 
                Tave_sm_1215, 
                PPT_wt_1215, 
                PPT_sm_1215)
drought.period <- as.data.frame(drought.period)

decline_clim_coeff <- data.frame()
for (i in 2:9){
  lm.1 <- lm(mean.lambda.drought~drought.period[,i],data=drought.period)
  #summary(lm.1)
  decline_clim_coeff[i-1,1] <- summary(lm.1)$coefficients[2,1] #Slope
  decline_clim_coeff[i-1,2] <- summary(lm.1)$coefficients[2,4] #P-value
  decline_clim_coeff[i-1,3] <- summary(lm.1)$adj.r.squared #R2
}

write_csv(decline_clim_coeff,"data/climate_data/decline_clim_coeff.csv") 





###################################################################################
#Plot decline vs significant env variables 

# N-S color gradient
lat_cols=colorRampPalette(brewer.pal(10,"Spectral"))
n.sites <- length(unique(demo_pop$Paper_ID)) -2
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


