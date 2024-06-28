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
demo_pop <- left_join(demog_means, anoms)  
  
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
#Population recovery and climate anomalies

a <- ggplot(demo_pop, aes(x=PPT_wt_1517, y=mean.lambda.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =6)+
  geom_smooth(method=lm,color="black", lty="dashed", se=FALSE)+
  scale_y_continuous(name="Mean Lambda After Drought", limits=c(0,2.5), breaks=seq(0,2.5,0.5))+
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






