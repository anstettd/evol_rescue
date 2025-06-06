##################################################################################
#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Correlate pop decline with climate
#### AUTHOR: Daniel Anstett and Amy Angert
#### DATE LAST MODIFIED: 20250606
###################################################################################

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

#Library install and import
library(tidyverse)
library(GGally)
library(Hmisc)
library(RColorBrewer)

#Import demography estimates (includes metadata)
demog_means <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv") %>% 
  mutate(log.decline = log(mean.lambda.drought+0.5),
         log.recovery = log(mean.lambda.recovery+0.5))

#Import climate anomalies
anoms <- read_csv("data/climate_data/climate_anomaly.csv") 

#Join and filter
demo_pop <- left_join(demog_means, anoms, by=c("Site", "Latitude", "Paper_ID")) #necessary to specify this list because Mill Creek longitude is not the same between files, so anom values turn to NA if default joining by all shared columns

drought.period <- demo_pop %>% 
  dplyr::select(log.decline, 
                MAT_1214, 
                MAP_1214, 
                PAS_1214, 
                CMD_1214, 
                Tave_wt_1214, 
                Tave_sm_1214, 
                PPT_wt_1214, 
                PPT_sm_1214)
drought.period <- as.data.frame(drought.period)

decline_clim_coeff <- data.frame()
for (i in 2:9){
  decline_clim_coeff[i-1,1] <- names(drought.period[i])
  lm.1 <- lm(scale(log.decline)~scale(drought.period[,i]),data=drought.period)
  #summary(lm.1)
  decline_clim_coeff[i-1,2] <- summary(lm.1)$coefficients[2,1] #Slope
  decline_clim_coeff[i-1,3] <- summary(lm.1)$coefficients[2,4] #P-value
  decline_clim_coeff[i-1,4] <- summary(lm.1)$adj.r.squared #R2
}

names(decline_clim_coeff) <- c("Var", "Slope", "P", "AdjR2")
write_csv(decline_clim_coeff,"data/climate_data/decline_clim_coeff.csv") 





###################################################################################
#Plot decline vs significant climate anomalies 

#Set color as factor
demo_pop$Lat.Color<-as.factor(demo_pop$Lat.Color)
demo_pop$Lat.Color<-factor(demo_pop$Lat.Color,levels=demo_pop$Lat.Color)


#Winter Precipitation Anomaly
a <- ggplot(demo_pop, aes(x=PPT_wt_1214, y=log.decline)) + 
  geom_smooth(method=lm,color="black",size=1.25,fill="gray71")+
  geom_point(aes(fill=Lat.Color), shape=21, size=6)+
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.7) +
  scale_y_continuous(name="Mean Pop. Growth During Drought")+
  scale_x_continuous(name="Winter Precipitation Anomaly (2012-2014)")+
  #,breaks=c(0.04,0.045,0.05,0.055,0.06))+
  scale_fill_manual(labels=round(demo_pop$Latitude,1), values=as.character(demo_pop$Lat.Color)) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=22, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce height
  )
a
ggsave("Graphs/Climate/1_drought_lambda_PPT_wt.pdf",width=8, height = 6, units = "in")
