##################################################################################
#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Test whether lambda recovery is predicted by selection slopes
#### AUTHOR: Daniel Anstett and Amy Angert
#### DATE LAST MODIFIED: 20250310
###################################################################################

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

#Library install and import
library(tidyverse) 
library(car)
library(ggrepel)
library(RColorBrewer)

#Import demography data + metadata
demo_pop <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv")

#Import selection data
slope.summary <- read_csv("data/snp_change_data/mean_median_S.csv") %>% select(Site, Median, Mean)

#Join frames
slope_pop <- left_join(demo_pop, slope.summary, by=c("Paper_ID"="Site"))

#Visualize scatter plot
ggplot(slope_pop, aes(x=Median, y=mean.lambda.recovery)) + geom_point()
# looks like an extreme outlier - needs test

#Cook's distance check for influential outliers
mod.check <- lm(mean.lambda.recovery~Median,data=slope_pop)
cooksd <- cooks.distance(mod.check)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels

#Cull outlier populations
slope_pop_cull1 <- slope_pop %>% 
  filter(Paper_ID!=4) %>% #remove influential outlier Redwood Creek
  droplevels()

#Cook's distance check for further influential outliers
mod.check <- lm(mean.lambda.recovery~Median,data=slope_pop_cull1)
cooksd_cull1 <- cooks.distance(mod.check)
plot(cooksd_cull1, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd_cull1, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd_cull1)+1, y=cooksd_cull1, labels=ifelse(cooksd_cull1>4*mean(cooksd_cull1, na.rm=T),names(cooksd_cull1),""), col="red")  # add labels

#Cull outlier population
slope_pop_cull2 <- slope_pop_cull1 %>% 
  filter(Paper_ID!=1) %>% #remove influential outlier Sweetwater R
  droplevels()

#Cook's distance check for further influential outliers
mod.check <- lm(mean.lambda.recovery~Median,data=slope_pop_cull2)
cooksd_cull2 <- cooks.distance(mod.check)
plot(cooksd_cull2, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd_cull2, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd_cull2)+1, y=cooksd_cull2, labels=ifelse(cooksd_cull2>4*mean(cooksd_cull2, na.rm=T),namescooksd_cull2,""), col="red")  # add labels


########################################################################################
#Lambda recovery vs median selection slope 


#### No removal #### 
slope_pop_graph <- drop_na(slope_pop)
lm3 <- lm(mean.lambda.recovery~Median,data=slope_pop_graph)
summary(lm3)
Anova(lm3,type="III")

slope_pop_graph$Lat.Color<-as.factor(slope_pop_graph$Lat.Color)
slope_pop_graph$Lat.Color<-factor(slope_pop_graph$Lat.Color,levels=slope_pop_graph$Lat.Color)

ggplot(slope_pop_graph, aes(x=Median, y=mean.lambda.recovery)) + 
  geom_smooth(method=lm,color="black",size=1.25,fill="gray71")+
  geom_point(aes(fill=slope_pop_graph$Lat.Color), shape=21, size =6)+
  geom_hline(yintercept = 1, linetype = "dotted", color = "black", size = 0.7) +
  scale_y_continuous(name="Mean Lambda after Drought")+#,breaks=c(0.5,1,1.5,2,2.5))+
  scale_x_continuous(name="Median Strength of Selection")+#,breaks=c(-0.1,0,0.1,0.2))+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(name = "Latitude (°N)",labels=round(slope_pop_graph$Latitude,1), values=as.character(slope_pop_graph$Lat.Color)) +
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
ggsave("Graphs/Demo_selection/1_median_slope_recovery_lambda.pdf",width=8, height = 6, units = "in")


#### Outlier removal #### 
slope_pop_graph_cull2 <- drop_na(slope_pop_cull2)
lm1 <- lm(mean.lambda.recovery~Median,data=slope_pop_graph_cull2)
summary(lm1)
Anova(lm1,type="III")

slope_pop_graph_cull2$Lat.Color<-as.factor(slope_pop_graph_cull2$Lat.Color)
slope_pop_graph_cull2$Lat.Color<-factor(slope_pop_graph_cull2$Lat.Color,levels=slope_pop_graph_cull2$Lat.Color)

ggplot(slope_pop_graph, aes(x=Median, y=mean.lambda.recovery)) + 
  geom_smooth(method=lm,color="black",size=1.8,fill="gray75")+
  geom_point(aes(fill=slope_pop_graph$Lat.Color), shape=21, size =6)+
  geom_hline(yintercept = 1, linetype = "dotted", color = "black", size = 0.7) +
  geom_smooth(data=slope_pop_graph_cull2, aes(x=Median, y=mean.lambda.recovery), method=lm, color="black", linetype="dashed", size=1.8, fill="gray50") +
  scale_y_continuous(name="Mean Lambda after Drought")+#,breaks=c(0.5,1,1.5,2,2.5))+
  scale_x_continuous(name="Median Strength of Selection")+#,breaks=c(-0.1,0,0.1,0.2))+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(name = "Latitude (°N)",labels=round(slope_pop_graph$Latitude,1),
                    values=as.character(slope_pop_graph$Lat.Color)) +
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
ggsave("Graphs/Demo_selection/2_median_slope_recovery_lambda.pdf",width=8, height = 6, units = "in")


## Overlay slopes with and without outlier removal

ggplot(slope_pop_graph_cull2, aes(x=Median, y=mean.lambda.recovery)) + 
  geom_smooth(method=lm,color="black",size=1.8,fill="gray75")+
  geom_point(aes(fill=slope_pop_graph_cull2$Lat.Color), shape=21, size =6)+
  geom_hline(yintercept = 1, linetype = "dotted", color = "black", size = 0.7) +
  geom_abline(slope=coef(lm3)[2], intercept=coef(lm3)[1]) + 
  scale_y_continuous(name="Mean Lambda after Drought")+#,breaks=c(0.5,1,1.5,2,2.5))+
  scale_x_continuous(name="Median Strength of Selection")+#,breaks=c(-0.1,0,0.1,0.2))+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(name = "Latitude (°N)",labels=round(slope_pop_graph_cull2$Latitude,1),
                    values=as.character(slope_pop_graph_cull2$Lat.Color)) +
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
ggsave("Graphs/Demo_selection/3_median_slope_recovery_lambda.pdf",width=8, height = 6, units = "in")


