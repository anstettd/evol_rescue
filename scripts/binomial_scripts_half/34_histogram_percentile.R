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
env_0_10 <- env_low %>% filter(percentile<=0.1)
median_pop_0_10 <- env_0_10 %>% group_by(Site) %>% summarise(median = median(Slope, na.rm = TRUE))

#10th to 100th percentile
env_10_100 <- env_low %>% filter(percentile>=0.1)
median_pop_10_100 <- env_10_100 %>% group_by(Site) %>% summarise(median = median(Slope, na.rm = TRUE))

#90th to 100th percentile
env_90_100 <- env_low %>% filter(percentile>=0.9)
median_pop_90_100 <- env_90_100 %>% group_by(Site) %>% summarise(median = median(Slope, na.rm = TRUE))

#Export Medians
median_table <- cbind(median_pop,median_pop_0_10[,2],median_pop_10_100[,2],median_pop_90_100[,2])
colnames(median_table) <- c("Site","Median","Median_0to10","Median_10to100","Median_90to100")
write_csv(median_table,"graphs/histograms/baseline_slope/median_table.csv")

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



#0th to 10th percentile
env_slope<- ggplot(env_0_10,aes(x=Slope))+
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
  geom_vline(data = median_pop_0_10, aes(xintercept = median), linetype="dashed",color="firebrick2")
env_slope
ggsave("graphs/histograms/baseline_slope/2worst_0_10th.pdf",width=10, height = 7.5, units = "in")



#10th to 100th percentile
env_slope<- ggplot(env_10_100,aes(x=Slope))+
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
  geom_vline(data = median_pop_10_100, aes(xintercept = median), linetype="dashed",color="firebrick2")
env_slope
ggsave("graphs/histograms/baseline_slope/3_most_10_100th.pdf",width=10, height = 7.5, units = "in")


#90th to 100th percentile
env_slope<- ggplot(env_90_100,aes(x=Slope))+
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
  geom_vline(data = median_pop_90_100, aes(xintercept = median), linetype="dashed",color="firebrick2")
env_slope
ggsave("graphs/histograms/baseline_slope/4_best_10_90th.pdf",width=10, height = 7.5, units = "in")




