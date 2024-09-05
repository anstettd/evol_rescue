##################################################################################
## Pi Latitude Plots
## Author Daniel Anstett
## 
## 
## Last Modified Jan 20, 2023
###################################################################################

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

#Library install and import
library(tidyverse) 
library(car)
library(ggrepel)
library(RColorBrewer)
library(cowplot)

###################################################################################
#Import data & Prepare data frame
#Baseline & Timeseries
all_pop <- read_csv("data/genomic_data/Baseline_Timeseries_pops_final2.csv")%>% filter(Paper_ID<56) %>% filter(Paper_ID!=12)

all_pop <- all_pop %>% mutate()

#Pi
pi_df <- read_csv("data/genomic_data/baseline_pi.csv")
#mean(pi_df$pi_all_snps)
pi_all_pop <-left_join(all_pop,pi_df,by=c("Paper_ID"="Site"))
###################################################################################

#stats
lm1 <- lm(pi_snp_set~Lat,data=pi_all_pop)
summary(lm1)
Anova(lm1,type="III")

lm2 <- lm(pi_snp_set~poly(Lat,2),data=pi_all_pop)
summary(lm2)
Anova(lm2,type="III")

lm3 <- lm(pi_all_snps~Lat,data=pi_all_pop)
summary(lm3)
Anova(lm3,type="III")

lm4 <- lm(pi_snp_set~pi_all_snps,data=pi_all_pop)
summary(lm4)
Anova(lm4,type="III")


###########################################################################################################
#Make Latitude-pi graphs

# N-S color gradient

#color.list <- c("Baseline"="grey40","Timeseries_Demography"="red","Timeseries"="cyan2","Demography"="orange")



#Pi (Climate-Associated) vs. Latitude
ggplot(pi_all_pop, aes(x=Lat, y=pi_snp_set,shape=Type)) + 
  #geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =6)+
  geom_point(aes(x=Lat,y=pi_snp_set),size =4)+
  #stat_smooth(method=lm,color="black",formula = y ~ x+I(x^2))+
  scale_y_continuous(name="Pi (Climate Associated)")+
  scale_x_continuous(name="Latitude")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_shape_manual(values = c(1,17,18,16))+
  #scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),#legend.position="none",
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce height
  )
ggsave("Graphs/Pi_latitude/1_lat_pi_snp_set.pdf",width=11, height = 5.5, units = "in")


#Pi (Genome-wide) vs. Latitude
ggplot(pi_all_pop, aes(x=Lat, y=pi_all_snps,shape=Type)) + 
  #geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =6)+
  geom_point(aes(x=Lat,y=pi_all_snps),size =4)+
  #stat_smooth(method =lm,color="black")+
  scale_y_continuous(name="Pi (Genome-Wide)")+
  scale_x_continuous(name="Latitude")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_shape_manual(values = c(1,17,18,16))+
  #scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(), #legend.position="none",
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce height
  )
ggsave("Graphs/Pi_latitude/2_lat_pi_global.pdf",width=11, height = 5.5, units = "in")


#Adaptive vs neutral Pi
ggplot(pi_all_pop, aes(x=pi_all_snps, y=pi_snp_set,shape=Type)) + 
  #geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =6)+
  geom_point(aes(x=pi_all_snps,y=pi_snp_set),size =4)+
  #stat_smooth(method =lm,color="black")+
  scale_y_continuous(name="Pi (Climate Associated)", breaks=c(0, 0.1, 0.2, 0.3, 0.4)) +
  scale_x_continuous(name="Pi (Genome-Wide)", breaks=c(0, 0.1, 0.2, 0.3)) +
  scale_shape_manual(values = c(1,17,18,16))+
  #scale_fill_manual(values=color.list) +
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
ggsave("Graphs/Pi_latitude/3_pi_vs_pi.pdf",width=11, height = 5.5, units = "in")




#Pi  
ggplot(pi_all_pop, aes(x=Lat, y=pi_snp_set,fill=Type)) + 
  #geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =6)+
  geom_point(aes(x=Lat,y=pi_snp_set),shape=21,size =4)+
  #stat_smooth(method=lm,color="black",formula = y ~ x+I(x^2))+
  scale_y_continuous(name="Pi (Climate Associated)")+
  scale_x_continuous(name="Latitude")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
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
ggsave("Graphs/Pi_latitude/legend.pdf",width=6, height = 5.5, units = "in")




####################### Old with colors ####################### 

#Pi  
ggplot(pi_all_pop, aes(x=Lat, y=pi_snp_set,fill=Type)) + 
  #geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =6)+
  geom_point(aes(x=Lat,y=pi_snp_set),shape=21,size =4)+
  #stat_smooth(method=lm,color="black",formula = y ~ x+I(x^2))+
  scale_y_continuous(name="Pi (Climate Associated)")+
  scale_x_continuous(name="Latitude")+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),legend.position="none",
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce height
  )
#ggsave("Graphs/Pi_latitude/1_lat_pi_snp_set.pdf",width=6, height = 5.5, units = "in")
