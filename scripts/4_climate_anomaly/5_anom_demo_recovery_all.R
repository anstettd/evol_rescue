##################################################################################
#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Correlate population recovery with climate
#### AUTHOR: Daniel Anstett and Amy Angert
#### DATE LAST MODIFIED: 20250607
###################################################################################

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

# Libraries
library(tidyverse)
library(GGally)
library(Hmisc)
library(RColorBrewer)
library(MASS)
library(sfsmisc)

#Import demography estimates + metadata
demog_means <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv")

#Import climate anomalies
anoms <- read_csv("data/climate_data/climate_anomaly.csv") 

#Join and filter
demo_pop <- left_join(demog_means, anoms, by=c("Site", "Latitude", "Paper_ID")) #necessary to specify this list because Mill Creek longitude is not the same between files, so anom values turn to NA if default joining by all shared columns  
  
recovery.period <- demo_pop %>% 
  dplyr::select(mean.r.recovery, 
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
  lm.1 <- lm(scale(mean.r.recovery)~scale(recovery.period[,i]),data=recovery.period)
  rlm.1 <- rlm(scale(mean.r.recovery)~scale(recovery.period[,i]),data=recovery.period, maxit=500)
  rlm.f <- f.robftest(rlm.1)
  #summary(lm.1)
  #summary(rlm.1)
  recovery_clim_coeff[i-1,2] <- summary(lm.1)$coefficients[2,1] #OLS slope
  recovery_clim_coeff[i-1,3] <- summary(lm.1)$coefficients[2,4] #OLS P-value
  recovery_clim_coeff[i-1,4] <- summary(lm.1)$adj.r.squared #OLS R2
  recovery_clim_coeff[i-1,5] <- rlm.1$coefficients[[2]] #robust slope
  recovery_clim_coeff[i-1,6] <- rlm.f$p.value #robust P-value
}
names(recovery_clim_coeff) <- c("Var", "OLS_Slope", "OLS_P", "OLS_AdjR2", "RR_Slope", "RR_P")

write_csv(recovery_clim_coeff,"data/climate_data/recovery_clim_coeff_alldemo.csv") 




###################################################################################
#Population recovery and climate anomalies 

demo_pop$Lat.Color<-as.factor(demo_pop$Lat.Color)
demo_pop$Lat.Color<-factor(demo_pop$Lat.Color,levels=demo_pop$Lat.Color)

a <- ggplot(demo_pop, aes(x=MAP_1517, y=mean.r.recovery)) +
  geom_smooth(method=lm,color="black",size=1.25, linetype="dashed", fill="gray71")+
  geom_smooth(method=MASS::rlm,color="black",size=1.25,fill="gray40")+
  geom_point(aes(fill=demo_pop$Lat.Color), shape=21, size =6)+
  geom_hline(yintercept = 1, linetype = "dotted", color = "black", size = 0.7) +
  scale_y_continuous(name="Mean Pop. Growth after Drought")+#, limits=c(0,2.5), breaks=seq(0,2.5,0.5))+
  scale_x_continuous(name="Annual Precipitation Anomaly (2015-2017)")+
  #,breaks=c(0.04,0.045,0.05,0.055,0.06))+
  #scale_fill_manual(name = "Latitude (°N)",labels=round(demo_pop$Latitude,1), values=as.character(demo_pop$Lat.Color)) +
  scale_fill_manual(values=as.character(demo_pop$Lat.Color), 
                    labels = unique(demo_pop$Paper_ID) ) +
#  scale_color_manual(values=demo_pop$Lat.Color, 
#                     labels = unique(demo_pop$Paper_ID )) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=22, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce height
  )+
  guides(color = guide_legend(reverse = TRUE, override.aes = list(linetype = 0)),
         fill  = guide_legend(reverse = TRUE))
a
ggsave("Graphs/Climate/2_recovery_r_MAP.pdf",width=8, height = 6, units = "in")




