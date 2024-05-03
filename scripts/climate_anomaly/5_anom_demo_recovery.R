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


#Import population metadata
pop_meta <- read_csv("data/genomic_data/pop_meta_data.csv") %>% dplyr::select(Site, Paper_ID) #just to get translation of pop names <--> numbers
pop_meta[20,1] <- "Mill Creek" #Fill in missing information
pop_meta[20,2] <- 12
demog_recovery <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv")
cumul <- read_csv("data/genomic_data/pop_meta_data.csv")
demog_recovery <- left_join(demog_recovery,pop_meta,by=c("Site"="Site")) %>% rename(Site_Name=Site)
demo_pop <- left_join(cumul,demog_recovery,by="Paper_ID") %>% filter(Paper_ID!=10) %>% filter(Paper_ID!=12)
anoms <- read_csv("data/climate_data/climate_anomaly.csv")
demo_pop <- left_join(demo_pop, anoms)


recovery.period <- demo_pop %>% 
  dplyr::select(lambda.slope.recovery, lambda.mean.recovery, 
                MAT_1619, MAP_1619, PAS_1619, CMD_1619, Tave_wt_1619, 
                Tave_sm_1619, PPT_wt_1619, PPT_sm_1619)

#Make Correlations
#ggpairs(recovery.period)
demo_recovery_anom <- rcorr(as.matrix(recovery.period))

#Export
#write.csv(demo_recovery_anom$r,"data/climate_data/rcorr_demo_recovery_anom_r.csv") 
#write.csv(demo_recovery_anom$P,"data/climate_data/rcorr_demo_recovery_anom_p.csv") 


###################################################################################
#Plot decline vs significant env variables 

# N-S color gradient
lat_cols=colorRampPalette(brewer.pal(10,"Spectral"))
n.sites <- length(unique(demo_pop$Paper_ID))
color.list <- lat_cols(n.sites)


#Population recovery and climate anomalies
a <- ggplot(demo_pop, aes(x=PPT_wt_1619, y=lambda.mean.drought)) + 
  geom_point(aes(fill=as.factor(round(Latitude.x, 1))),shape=21,size =6)+
  geom_smooth(method=lm,color="black", lty="dashed", se=FALSE)+
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
a
ggsave("Graphs/Climate/mean_lambda_PPT_wt.pdf",width=8, height = 6, units = "in")






