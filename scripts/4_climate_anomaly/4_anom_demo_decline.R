##################################################################################
#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Correlate pop decline with climate
#### AUTHOR: Daniel Anstett and Amy Angert
#### DATE LAST MODIFIED: 20251211
###################################################################################

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

#Library install and import
library(tidyverse)
library(GGally)
library(Hmisc)
library(RColorBrewer)
library(MASS)
library(sfsmisc)
library(plotly)
#Import demography estimates (includes metadata)
demog_means <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv") 

#Import climate anomalies
anoms <- read_csv("data/climate_data/climate_anomaly.csv") 

#Join and filter
demo_pop <- left_join(demog_means, anoms, by=c("Site", "Latitude", "Paper_ID")) #necessary to specify this list because Mill Creek longitude is not the same between files, so anom values turn to NA if default joining by all shared columns

drought.period <- demo_pop %>% 
  dplyr::select(mean.r.drought, 
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
  lm.1 <- lm(scale(mean.r.drought)~scale(drought.period[,i]),data=drought.period)
  rlm.1 <- rlm(scale(mean.r.drought)~scale(drought.period[,i]),data=drought.period, maxit=500)
  rlm.f <- f.robftest(rlm.1)
  #summary(lm.1)
  #summary(rlm.1)
  decline_clim_coeff[i-1,2] <- summary(lm.1)$coefficients[2,1] #Slope
  decline_clim_coeff[i-1,3] <- summary(lm.1)$coefficients[2,4] #P-value
  decline_clim_coeff[i-1,4] <- summary(lm.1)$adj.r.squared #R2
  decline_clim_coeff[i-1,5] <- rlm.1$coefficients[[2]] #robust slope
  decline_clim_coeff[i-1,6] <- rlm.f$p.value #robust P-value
}

names(decline_clim_coeff) <- c("Var", "OLS_Slope", "OLS_P", "OLS_AdjR2", "RR_Slope", "RR_P")
write_csv(decline_clim_coeff,"data/climate_data/decline_clim_coeff.csv") 





###################################################################################
#Plot demographic decline vs significant climate anomalies 

#Set color as factor
demo_pop$Lat.Color<-as.factor(demo_pop$Lat.Color)
demo_pop$Lat.Color<-factor(demo_pop$Lat.Color,levels=demo_pop$Lat.Color)


# Demographic decline ~ Winter Precipitation Anomaly (FIGURE 1E)
a <- ggplot(demo_pop, aes(x=PPT_wt_1214, y=mean.r.drought)) + 
  geom_smooth(method=lm,color="black",linetype="dashed",size=1.25,fill="gray75")+
  geom_smooth(method=MASS::rlm,color="black",size=1.25,fill="gray50")+
  geom_point(aes(fill=Lat.Color), shape=21, size=6)+
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.7) +
  scale_y_continuous(name="Mean Pop. Growth During Drought")+
  scale_x_continuous(name="Winter Precipitation Anomaly (2012-2014)")+
  #,breaks=c(0.04,0.045,0.05,0.055,0.06))+
  #scale_fill_manual(labels=round(demo_pop$Latitude,1), 
  #                 values=as.character(demo_pop$Lat.Color)) +
  scale_fill_manual(
    name = "Latitude (°N)",
    values = as.character(demo_pop$Lat.Color),  # colors in order
    labels = demo_pop$Paper_ID                   # labels replacing latitude
  )+
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_classic() + theme(
    axis.text.x = element_text(size=24, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=22, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce height
  )
a
ggsave("Graphs/Climate/1_drought_meanr_PPTwt_anom_1E.pdf",width=8, height = 7, units = "in")



#Demographic decline ~ Summer Temperature Anomaly (not shown in manuscript)
b <- ggplot(demo_pop, aes(x=Tave_sm_1214, y=mean.r.drought)) + 
  geom_smooth(method=lm,color="black",size=1.25,linetype="dashed", fill="gray75")+
  geom_smooth(method=MASS::rlm,color="black",size=1.25, fill="gray50")+
  geom_point(aes(fill=Lat.Color), shape=21, size=6)+
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.7) +
  scale_y_continuous(name="Mean Pop. Growth During Drought")+
  scale_x_continuous(name="Summer Temperature Anomaly (2012-2014)")+
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
b
ggsave("Graphs/Climate/1_drought_meanr_Tavesm_anom.pdf",width=8, height = 6, units = "in")

#Demographic decline ~ Climatic Moisture Deficit Anomaly (not shown in manuscript)
c <- ggplot(demo_pop, aes(x=CMD_1214, y=mean.r.drought)) + 
  geom_smooth(method=lm,color="black",size=1.25,linetype="dashed",fill="gray75")+
  geom_smooth(method=MASS::rlm,color="black",size=1.25,fill="gray50")+
  geom_point(aes(fill=Lat.Color), shape=21, size=6)+
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.7) +
  scale_y_continuous(name="Mean Pop. Growth During Drought")+
  scale_x_continuous(name="Moisture Deficit Anomaly (2012-2014)")+
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
c
ggsave("Graphs/Climate/1_drought_meanr_CMD_anom.pdf",width=8, height = 6, units = "in")
