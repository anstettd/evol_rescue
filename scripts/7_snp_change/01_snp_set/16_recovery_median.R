##################################################################################
#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Test whether lambda recovery is predicted by selection slopes
#### AUTHOR: Daniel Anstett and Amy Angert
#### DATE LAST MODIFIED: 20250608
###################################################################################

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

#Library install and import
library(tidyverse) 
library(car)
library(ggrepel)
library(RColorBrewer)
library(MASS)
library(sfsmisc)

#Import demography data + metadata
demo_pop <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv")

#Import selection data
slope.summary <- read_csv("data/snp_change_2/mean_median_S_all.csv") %>% dplyr::select(Site, Median, Mean)

#Join frames
slope_pop <- left_join(demo_pop, slope.summary, by=c("Paper_ID"="Site"))
slope_pop <- slope_pop %>% dplyr::select(Latitude,Site,Paper_ID,Lat.Color,mean.r.recovery,Median)
slope_pop <- drop_na(slope_pop)

#Visualize scatter plot
ggplot(slope_pop, aes(x=Median, y=mean.r.recovery)) + geom_point()

mod_S <- lm(mean.r.recovery~Median,data=slope_pop)
summary(mod_S) # P = 0.011, R2 = 0.52
Anova(mod_S,type="III")
qqnorm(resid(mod_S))
qqline(resid(mod_S))

#Robust regression instead
rob.mod_S <- rlm(slope_pop$mean.r.recovery~Median,data=slope_pop)
summary(rob.mod_S) 
f.robftest(rob.mod_S, var="Median") #p-value = 0.01238


######################################################################################### Graphs of lambda recovery vs median selection slope 

slope_pop$Lat.Color<-as.factor(slope_pop$Lat.Color)
slope_pop$Lat.Color<-factor(slope_pop$Lat.Color,levels=slope_pop$Lat.Color)

ggplot(slope_pop, aes(x=Median, y=mean.r.recovery)) + 
  geom_smooth(method=lm,color="black", size=1.8, linetype="dashed", fill="gray75")+
  geom_smooth(method=MASS::rlm, color="black", size=1.8, fill="grey50") +
  #geom_smooth(data=slope_pop_cull1, aes(x=Median, y=log(mean.lambda.recovery+0.5)), method=lm, color="black", size=1.8, linetype="longdash", fill="grey35") +
  geom_point(aes(fill=Lat.Color), shape=21, size =6)+
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.7) +
  scale_y_continuous(name="Mean Pop. Growth after Drought")+#,breaks=c(0.5,1,1.5,2,2.5))+
  scale_x_continuous(name="Median Response to Selection")+#,breaks=c(-0.1,0,0.1,0.2))+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(name = "Latitude (°N)",labels=round(slope_pop$Latitude,1),
                    values=as.character(slope_pop$Lat.Color)) +
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
ggsave("Graphs/Demography_2/01_median_slope_logr.pdf",width=8, height = 6, units = "in")





