##################################################################################
## Plot lambda mean against slope summaries
## Author Daniel Anstett
## 
## 
## Last Modified Aug 14, 2023
###################################################################################
#Library install and import
library(tidyverse) 
library(car)
library(ggrepel)
library(RColorBrewer)

#Import data & Prepare data frame
pop_meta <- read_csv("data/genomic_data/pop_meta_data.csv") %>% dplyr::select(Site, Paper_ID) #just to get translation of pop names <--> numbers
pop_meta[20,1] <- "Mill Creek"
pop_meta[20,2] <- 12
demog_recovery <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv")
demog_recovery <- left_join(demog_recovery,pop_meta,by=c("Site"="Site")) %>% rename(Site_Name=Site)

#Import Mean Median Selection data
slope.summary <- read_csv("data/snp_change_data/mean_median_S.csv") %>%  select(Site,Median,Mean)

#Remove unneeded sites
demo_pop <- left_join(slope.summary,demog_recovery,by=c("Site"="Paper_ID")) %>% filter(Site!=10) %>% 
  filter(Site!=12) 


#stats
lm1 <- lm(lambda.mean.recovery~Median,data=demo_pop)
summary(lm1)
Anova(lm1,type="III")

lm2 <- lm(lambda.mean.recovery~Mean,data=demo_pop)
summary(lm2)
Anova(lm2,type="III")

lm3 <- lm(Median~Latitude,data=demo_pop)
summary(lm3)
Anova(lm3,type="III")


###########################################################################################################
# N-S color gradient
#lat_cols=colorRampPalette(brewer.pal(10,"Spectral"))
#n.sites <- length(unique(demo_pop$Site))
#color.list <- lat_cols(n.sites)


###########################################################################################################

# N-S color gradient
lat_cols=c("#DC494C","#F88D51","#FDD380","#FEEB9E","#FFFFBF","#D7EF9B","#B2E0A2","#88CFA4","#5FBAA8","#3F96B7")



#Median slope vs. lambda.slope.recovery
ggplot(demo_pop, aes(x=Median, y=lambda.mean.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =6)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Mean Lambda after Drought",breaks=c(0.5,1,1.5,2,2.5))+
  scale_x_continuous(name="Median Slope",breaks=c(-0.1,0,0.1,0.2))+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=lat_cols) +
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
ggsave("Graphs/Selection_demo/1_median_slope_recovery_lambda.pdf",width=8, height = 6, units = "in")

#Median slope vs. lambda.slope.recovery
ggplot(demo_pop, aes(x=Mean, y=lambda.mean.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =6)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Mean Lambda after Drought",breaks=c(0.5,1,1.5,2,2.5))+
  scale_x_continuous(name="Mean Slope",breaks=c(-0.1,0,0.1,0.2))+
  scale_fill_manual(values=lat_cols) +
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
ggsave("Graphs/Selection_demo/2_mean_slope_recovery_lambda.pdf",width=8, height = 6, units = "in")


#Median slope vs. lambda.slope.recovery
ggplot(demo_pop, aes(x=Latitude, y=Median)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =6)+
  geom_smooth(method=lm,color="black", lty="dashed", se=FALSE)+
  scale_y_continuous(name="Median",breaks=c(-0.1,0,0.1,0.2),limits = c(-0.1,0.2))+
  scale_x_continuous(name="Latitude")+
  scale_fill_manual(values=lat_cols) +
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

ggsave("Graphs/Selection_demo/3_median_latitude.pdf",width=8, height = 6, units = "in")






