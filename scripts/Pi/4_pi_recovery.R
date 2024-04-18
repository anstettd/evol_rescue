##################################################################################
## Plot pi against lambda.slope & lambda mean
## Author Daniel Anstett
## 
## 
## Last Modified Aug 3, 2023
###################################################################################

#Library install and import
library(tidyverse) 
library(car)
library(ggrepel)
library(RColorBrewer)
library(cowplot)

#Import data
offset_pop <- read_csv("data/genomic_data/offset_pop_beagle.csv") %>% dplyr::select(Site, Paper_ID) #just to get translation of pop names <--> numbers
offset_pop[20,1] <- "Mill Creek"
offset_pop[20,2] <- 12
demog_recovery <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv")
pi_raw <- read_csv("data/genomic_data/raw_pi.csv")

demog_recovery <- left_join(demog_recovery,offset_pop,by=c("Site"="Site")) %>% rename(Site_Name=Site)
pi_pop <- left_join(demog_recovery,pi_raw,by=c("Paper_ID"="Site")) %>% filter(Paper_ID!=12) 


#stats
lm1 <- lm(lambda.slope.recovery~pi_snp_set,data=pi_pop)
lm2 <- lm(lambda.slope.recovery~pi_all_snps,data=pi_pop)
lm3 <- lm(lambda.mean.recovery~pi_snp_set,data=pi_pop)
lm4 <- lm(lambda.mean.recovery~pi_all_snps,data=pi_pop)
lm5 <- lm(pi_snp_set~Latitude,data=pi_pop)
lm5b <- lm(pi_snp_set~poly(Latitude,2),data=pi_pop)
lm6 <- lm(pi_all_snps~Latitude,data=pi_pop)


summary(lm1)
Anova(lm1,type="III")

summary(lm2)
Anova(lm2,type="III")

summary(lm3)
Anova(lm3,type="III")

summary(lm4)
Anova(lm4,type="III")

summary(lm5)
Anova(lm5,type="III")

summary(lm5b)
Anova(lm5b,type="III")

summary(lm6)
Anova(lm6,type="III")









###########################################################################################################

# N-S color gradient
lat_cols=colorRampPalette(brewer.pal(11,"Spectral"))
n.sites <- length(unique(pi_pop$Site_Name)) -2
color.list <- lat_cols(n.sites)


#Lambda Slope

#pi snp set
ggplot(pi_pop, aes(x=pi_snp_set, y=lambda.slope.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Rate of Increase in Lambda")+
  scale_x_continuous(name="Pi (Climate SNP)")+
                     #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce height
  )
ggsave("Graphs/lambda_pi/Recovery/1_pi_lambda_recovery_snpset.pdf",width=8, height = 6, units = "in")


#global pi
ggplot(pi_pop, aes(x=pi_all_snps, y=lambda.slope.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Rate of Increase in Lambda")+
  scale_x_continuous(name="Pi (Genome-Wide)")+
                     #,breaks=c(0.04,0.045,0.05,0.055,0.06))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce hight
  )
ggsave("Graphs/lambda_pi/Recovery/2_pi_lambda_recovery_global.pdf",width=8, height = 6, units = "in")


#Mean Lambda

#pi snp set
ggplot(pi_pop, aes(x=pi_snp_set, y=lambda.mean.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =6)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Mean Lambda after Drought",breaks=seq(0,5,1))+
  scale_x_continuous(name="Pi (Climate SNP)")+
                    # ,breaks=c(0.025,0.03,0.035,0.04,0.045))+
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
ggsave("Graphs/lambda_pi/Recovery/3_pi_mean_lambda_recovery_snpset.pdf",width=8, height = 6, units = "in")


#global pi
ggplot(pi_pop, aes(x=pi_all_snps, y=lambda.mean.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =6)+
  geom_smooth(method=lm,color="black", lty="dashed", se=FALSE)+
  scale_y_continuous(name="Mean Lambda after Drought")+
  scale_x_continuous(name="Pi (Genome-Wide)")+
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

ggsave("Graphs/lambda_pi/Recovery/4_pi_mean_lambda_recovery_global.pdf",width=8, height = 6, units = "in")


#global pi vs Latitude
ggplot(pi_pop, aes(x=Latitude, y=pi_all_snps)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Genome-wide Pi")+
  scale_x_continuous(name="Latitude")+
  #,breaks=c(0.04,0.045,0.05,0.055,0.06))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce hight
  )
ggsave("Graphs/lambda_pi/Recovery/5_pi_global_latitude.pdf",width=8, height = 6, units = "in")

#climate pi vs Latitude
ggplot(pi_pop, aes(x=Latitude, y=pi_snp_set)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Climate-associated Pi")+
  scale_x_continuous(name="Latitude")+
  #,breaks=c(0.04,0.045,0.05,0.055,0.06))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce hight
  )
ggsave("Graphs/lambda_pi/Recovery/6_pi_climate_latitude.pdf",width=8, height = 6, units = "in")







