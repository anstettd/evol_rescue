##################################################################################
#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Correlate population recovery with climate
#### AUTHOR: Daniel Anstett and Amy Angert
#### DATE LAST MODIFIED: 20240628
###################################################################################

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

# Libraries
library(tidyverse)
library(GGally)
library(Hmisc)
library(RColorBrewer)

#Import demography estimates + metadata
demog_means <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv") 

#Import climate anomalies
anoms <- read_csv("data/climate_data/climate_anomaly.csv") 

#Join and filter
demo_pop <- left_join(demog_means, anoms, by=c("Site", "Latitude", "Paper_ID")) #necessary to specify this list because Mill Creek longitude is not the same between files, so anom values turn to NA if default joining by all shared columns  
  
recovery.period <- demo_pop %>% 
  dplyr::select(mean.lambda.recovery, 
                MAT_1517, 
                MAP_1517, 
                PAS_1517, 
                CMD_1517, 
                Tave_wt_1517, 
                Tave_sm_1517, 
                PPT_wt_1517, 
                PPT_sm_1517)

recovery.period <- as.data.frame(recovery.period)

recovery_clim_coeff <- data.frame()
for (i in 2:9){
  recovery_clim_coeff[i-1,1] <- names(recovery.period[i])
  lm.1 <- lm(mean.lambda.recovery~recovery.period[,i],data=recovery.period)
  #summary(lm.1)
  recovery_clim_coeff[i-1,2] <- summary(lm.1)$coefficients[2,1] #Slope
  recovery_clim_coeff[i-1,3] <- summary(lm.1)$coefficients[2,4] #P-value
  recovery_clim_coeff[i-1,4] <- summary(lm.1)$adj.r.squared #R2
}
names(recovery_clim_coeff) <- c("Var", "Slope", "P", "AdjR2")

write_csv(recovery_clim_coeff,"data/climate_data/recovery_clim_coeff_alldemo.csv") 




###################################################################################
#Population recovery and climate anomalies (no outliers culled)

a <- ggplot(demo_pop, aes(x=PPT_wt_1517, y=mean.lambda.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =6)+
  geom_smooth(method=lm,color="black", lty="dashed", se=FALSE)+
  scale_y_continuous(name="Mean Lambda After Drought")+#, limits=c(0,2.5), breaks=seq(0,2.5,0.5))+
  scale_x_continuous(name="Winter Precipitation Anomaly (2015-2017)")+
  #,breaks=c(0.04,0.045,0.05,0.055,0.06))+
  scale_fill_manual(values=demo_pop$Lat.Color) +
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
a
ggsave("Graphs/Climate/recovery_lambda_PPT_wt.pdf",width=8, height = 6, units = "in")

b <- ggplot(demo_pop, aes(x=CMD_1517, y=mean.lambda.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =6)+
  geom_smooth(method=lm,color="black", lty="dashed", se=FALSE)+
  scale_y_continuous(name="Mean Lambda After Drought")+#, limits=c(0,2.5), breaks=seq(0,2.5,0.5))+
  scale_x_continuous(name="Climatic Moisure Deficit Anomaly (2015-2017)")+
  #,breaks=c(0.04,0.045,0.05,0.055,0.06))+
  scale_fill_manual(values=demo_pop$Lat.Color) +
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
b
ggsave("Graphs/Climate/recovery_lambda_CMD.pdf",width=8, height = 6, units = "in")



#We should probably subject these to the same outlier checks as for pi and slope analyses?

#Cook's distance check for influential outliers
mod.check <- lm(mean.lambda.recovery~PPT_wt_1517,data=demo_pop)
cooksd <- cooks.distance(mod.check)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels

#Cull outlier populations
demo_pop_cull1 <- demo_pop %>% 
  filter(Paper_ID!=12) %>% #remove influential outlier Mill Creek
  filter(Paper_ID!=27) #remove influential outlier Buck Meadows

#Cook's distance check for further influential outliers
mod.check_cull1 <- lm(mean.lambda.recovery~PPT_wt_1517,data=demo_pop_cull1)
cooksd_cull1 <- cooks.distance(mod.check_cull1)
plot(cooksd_cull1, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd_cull1, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd_cull1)+1, y=cooksd_cull1, labels=ifelse(cooksd_cull1>4*mean(cooksd_cull1, na.rm=T),names(cooksd_cull1),""), col="red")  # add labels

#Cull further outlier populations
demo_pop_cull2 <- demo_pop_cull1 %>% 
  filter(Paper_ID!=15) %>% #remove influential outlier Kitchen Creek
  filter(Paper_ID!=4) #remove influential outlier Redwood

#Cook's distance check for further influential outliers
mod.check_cull2 <- lm(mean.lambda.recovery~PPT_wt_1517,data=demo_pop_cull2)
cooksd_cull2 <- cooks.distance(mod.check_cull2)
plot(cooksd_cull2, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd_cull2, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd_cull2)+1, y=cooksd_cull2, labels=ifelse(cooksd_cull2>4*mean(cooksd_cull2, na.rm=T),names(cooksd_cull2),""), col="red")  # add labels

#Cull further outlier populations
demo_pop_cull3 <- demo_pop_cull2 %>% 
  filter(Paper_ID!=14) %>% #remove influential outlier Hauser
  filter(Paper_ID!=1) #remove influential outlier Sweetwater

#Cook's distance check for further influential outliers
mod.check_cull3 <- lm(mean.lambda.recovery~PPT_wt_1517,data=demo_pop_cull3)
cooksd_cull3 <- cooks.distance(mod.check_cull3)
plot(cooksd_cull3, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd_cull3, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd_cull3)+1, y=cooksd_cull3, labels=ifelse(cooksd_cull3>4*mean(cooksd_cull3, na.rm=T),names(cooksd_cull3),""), col="red")  # add labels

#Cull further outlier populations
demo_pop_cull4 <- demo_pop_cull3 %>% 
  filter(Paper_ID!=58)  #remove influential outlier SFMF Tule
 
#Cook's distance check for further influential outliers
mod.check_cull4 <- lm(mean.lambda.recovery~PPT_wt_1517,data=demo_pop_cull4)
cooksd_cull4 <- cooks.distance(mod.check_cull4)
plot(cooksd_cull4, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd_cull4, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd_cull4)+1, y=cooksd_cull4, labels=ifelse(cooksd_cull4>4*mean(cooksd_cull4, na.rm=T),names(cooksd_cull4),""), col="red")  # add labels

#Cull further outlier populations
demo_pop_cull5 <- demo_pop_cull4 %>% 
  filter(Paper_ID!=17)  #remove influential outlier Whitewater

#Cook's distance check for further influential outliers
mod.check_cull5 <- lm(mean.lambda.recovery~PPT_wt_1517,data=demo_pop_cull5)
cooksd_cull5 <- cooks.distance(mod.check_cull5)
plot(cooksd_cull5, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd_cull5, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd_cull5)+1, y=cooksd_cull5, labels=ifelse(cooksd_cull5>4*mean(cooksd_cull5, na.rm=T),names(cooksd_cull5),""), col="red")  # add labels
#no more!


recovery.period.cull <- demo_pop_cull5 %>% 
  dplyr::select(mean.lambda.recovery, 
                MAT_1517, 
                MAP_1517, 
                PAS_1517, 
                CMD_1517, 
                Tave_wt_1517, 
                Tave_sm_1517, 
                PPT_wt_1517, 
                PPT_sm_1517)

recovery.period.cull <- as.data.frame(recovery.period.cull)

recovery_clim_coeff_cull <- data.frame()
for (i in 2:9){
  recovery_clim_coeff_cull[i-1,1] <- names(recovery.period.cull[i])
  lm.1 <- lm(mean.lambda.recovery~recovery.period.cull[,i],data=recovery.period.cull)
  #summary(lm.1)
  recovery_clim_coeff_cull[i-1,2] <- summary(lm.1)$coefficients[2,1] #Slope
  recovery_clim_coeff_cull[i-1,3] <- summary(lm.1)$coefficients[2,4] #P-value
  recovery_clim_coeff_cull[i-1,4] <- summary(lm.1)$adj.r.squared #R2
}
names(recovery_clim_coeff_cull) <- c("Var", "Slope", "P", "AdjR2")

write_csv(recovery_clim_coeff_cull,"data/climate_data/recovery_clim_coeff_culleddemo.csv") 


