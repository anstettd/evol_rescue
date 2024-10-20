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

#Import data
offset_pop <- read_csv("data/genomic_data/offset_pop_beagle.csv") %>% select(Site, Paper_ID) #just to get translation of pop names <--> numbers
offset_pop[20,1] <- "Mill Creek"
offset_pop[20,2] <- 12
demog_decline <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv")
pi_raw <- read_csv("data/genomic_data/raw_pi.csv")

demog_decline <- left_join(demog_decline,offset_pop,by=c("Site"="Site")) %>% rename(Site_Name=Site)
pi_pop <- left_join(demog_decline,pi_raw,by=c("Paper_ID"="Site")) %>% filter(Paper_ID!=12) 


#stats
lm1 <- lm(lambda.slope.decline~pi_snp_set,data=pi_pop)
lm2 <- lm(lambda.slope.decline~pi_all_snps,data=pi_pop)
lm3 <- lm(lambda.mean.drought~pi_snp_set,data=pi_pop)
lm4 <- lm(lambda.mean.drought~pi_all_snps,data=pi_pop)


summary(lm1)
Anova(lm1,type="III")

summary(lm2)
Anova(lm2,type="III")

summary(lm3)
Anova(lm3,type="III")

summary(lm4)
Anova(lm4,type="III")











###########################################################################################################

# N-S color gradient
lat_cols=colorRampPalette(brewer.pal(11,"Spectral"))
n.sites <- length(unique(pi_pop$Site_Name))
color.list <- lat_cols(n.sites)


#Lambda Slope

#pi snp set
ggplot(pi_pop, aes(x=pi_snp_set, y=lambda.slope.decline)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Rate of Decline in Lambda")+
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
ggsave("Graphs/lambda_pi/Decline/1_pi_lambda_decline_snpset.pdf",width=8, height = 6, units = "in")


#global pi
ggplot(pi_pop, aes(x=pi_all_snps, y=lambda.slope.decline)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Rate of Decline in Lambda")+
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
ggsave("Graphs/lambda_pi/Decline/2_pi_lambda_decline_global.pdf",width=8, height = 6, units = "in")


#Mean Lambda

#pi snp set
ggplot(pi_pop, aes(x=pi_snp_set, y=lambda.mean.drought)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Mean Lambda during Drought")+
  scale_x_continuous(name="Pi (Climate SNP)")+
                    # ,breaks=c(0.025,0.03,0.035,0.04,0.045))+
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
ggsave("Graphs/lambda_pi/Decline/3_pi_mean_lambda_drought_snpset.pdf",width=8, height = 6, units = "in")


#global pi
ggplot(pi_pop, aes(x=pi_all_snps, y=lambda.mean.drought)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Mean Lambda during Drought")+
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
ggsave("Graphs/lambda_pi/Decline/4_pi_mean_lambda_drought_global.pdf",width=8, height = 6, units = "in")
