##################################################################################
## Make SNP slopes histograms for obs data only
## Author Daniel Anstett
## 
## For all 9 env
## Last Modified April 25, 2023
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)

#Full Data
env_all <- read_csv("data/binomial_data_half/slope_obs_all_unique_base_slope.csv")
env_low <- env_all %>% filter(SE<5) %>% filter(SE_base<5)
median_pop <- env_low %>% group_by(Site) %>% summarise(median = median(Slope, na.rm = TRUE))
#env_low$Site <- as.factor(env_low$Site)
#env_low$Site <- factor(env_low$Site, levels = c("1","12","2","3","4","5","6","7","8","9","10","11"))


#0th to 10th percentile

#Remove < 0.1 slope
env_01 <- env_low %>% filter(Slope_base>=0.1)
median_pop_01 <- env_01 %>% group_by(Site) %>% summarise(median = median(Slope, na.rm = TRUE))

#Remove < 0.2 slope
env_02 <- env_low %>% filter(Slope_base>=0.2)
median_pop_02 <- env_02 %>% group_by(Site) %>% summarise(median = median(Slope, na.rm = TRUE))

#Remove < 0.3 slope
env_03<- env_low %>% filter(Slope_base>=0.3)
median_pop_03 <- env_03 %>% group_by(Site) %>% summarise(median = median(Slope, na.rm = TRUE))



#Export Medians
median_table <- cbind(median_pop,median_pop_01[,2],median_pop_02[,2],median_pop_03[,2])
colnames(median_table) <- c("Site","Median","Median_0.1","Median_0.2","Median_0.3")
write_csv(median_table,"graphs/histograms/baseline_slope/median_table_rm_baseslope.csv")

###################################################################################

#All baseline slope percentiles
env_slope<- ggplot(env_low,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",fill="lightblue1",alpha=0.5)+
  geom_vline(xintercept=0) +
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection",limits=c(-2.2,2))+
  theme_classic()
env_slope <- env_slope  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
env_slope <- env_slope + facet_wrap(.~Site) +   
  geom_vline(data = median_pop, aes(xintercept = median), linetype="dashed",color="firebrick2")
env_slope
ggsave("graphs/histograms/baseline_slope/1_all.pdf",width=10, height = 7.5, units = "in")



#less_0.1
env_slope<- ggplot(env_01,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",fill="lightblue1",alpha=0.5)+
  geom_vline(xintercept=0) +
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection",limits=c(-2.2,2))+
  theme_classic()
env_slope <- env_slope  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
env_slope <- env_slope + facet_wrap(.~Site) +   
  geom_vline(data = median_pop_01, aes(xintercept = median), linetype="dashed",color="firebrick2")
env_slope
ggsave("graphs/histograms/baseline_slope/2_less_0.1.pdf",width=10, height = 7.5, units = "in")



#less_0.2
env_slope<- ggplot(env_02,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",fill="lightblue1",alpha=0.5)+
  geom_vline(xintercept=0) +
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection",limits=c(-2.2,2))+
  theme_classic()
env_slope <- env_slope  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
env_slope <- env_slope + facet_wrap(.~Site) +   
  geom_vline(data = median_pop_02, aes(xintercept = median), linetype="dashed",color="firebrick2")
env_slope
ggsave("graphs/histograms/baseline_slope/3_most_less_0.2.pdf",width=10, height = 7.5, units = "in")


#less_0.3
env_slope<- ggplot(env_03,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",fill="lightblue1",alpha=0.5)+
  geom_vline(xintercept=0) +
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection",limits=c(-2.2,2))+
  theme_classic()
env_slope <- env_slope  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
env_slope <- env_slope + facet_wrap(.~Site) +   
  geom_vline(data = median_pop_03, aes(xintercept = median), linetype="dashed",color="firebrick2")
env_slope
ggsave("graphs/histograms/baseline_slope/4_best_less_0.3.pdf",width=10, height = 7.5, units = "in")




