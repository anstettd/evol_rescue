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

#Import population metadata
pop_meta <- read_csv("data/genomic_data/pop_meta_data.csv") %>% dplyr::select(Site, Paper_ID) #just to get translation of pop names <--> numbers
pop_meta[20,1] <- "Mill Creek" #Fill in missing information
pop_meta[20,2] <- 12

#Import Demography Data
demog_recovery <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv")

#Import PI for demography popoulations
pi_raw <- read_csv("data/genomic_data/raw_pi.csv")

#Joing Data Sets
demog_recovery <- left_join(demog_recovery,pop_meta,by=c("Site"="Site")) %>% rename(Site_Name=Site)
pi_pop <- left_join(demog_recovery,pi_raw,by=c("Paper_ID"="Site")) %>% 
  filter(Paper_ID!=12) #%>% #remove site with unreliable demography data
  #filter(Paper_ID!=27) 


#Run linear models
lm3 <- lm(lambda.mean.recovery~pi_snp_set,data=pi_pop)
lm4 <- lm(lambda.mean.recovery~pi_all_snps,data=pi_pop)

#Get summary and Anova for each model
summary(lm3)
Anova(lm3,type="III")
summary(lm4)
Anova(lm4,type="III")


###########################################################################################################

# N-S color gradient
lat_cols=colorRampPalette(brewer.pal(11,"Spectral"))
n.sites <- length(unique(pi_pop$Site_Name)) -2
color.list <- lat_cols(n.sites)


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
ggsave("Graphs/Demography_pi/1_pi_demography_snpset.pdf",width=8, height = 6, units = "in")


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

ggsave("Graphs/Demography_pi/2_pi_demography_global.pdf",width=8, height = 6, units = "in")






