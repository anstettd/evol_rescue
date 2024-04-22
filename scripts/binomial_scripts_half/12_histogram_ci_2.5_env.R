##################################################################################
## Generate input for strength of selection graphs
## Processing for both observation and permuted CI
## Author Daniel Anstett
## 
## Modified from https://jkzorz.github.io/2020/05/17/Error-bars.html
## Last Modified May 17, 2022
###################################################################################
#Function

theme_ci <- function(){ 
  theme_classic() %+replace%    #replace elements we want to change
    theme(axis.text.x = element_text(size = 14, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
          axis.title = element_text(size = 18, face = "bold"), 
          axis.text.y = element_text(size = 16, face = "bold"),
          strip.background = element_blank(),strip.text.x = element_text(size = 16, face = "bold"))
}



###################################################################################
#Import libraries
library(tidyverse)

#Import files
env_obs_ci <- read_csv("data/binomial_data_half/obs_ci_env_ab.csv")  %>% filter(S <= 1.25 & S>= -1.25) 

#Import medians
median_pop_env <- read_csv("data/binomial_data_half/median_pop_env.csv")


###################################################################################
## Histogram graphs with CI
###################################################################################

#make input dataframe tidy

tidy_env_obs_ci<- env_obs_ci %>% pivot_longer(starts_with("p"), names_to = "Site", values_to = "Count") %>%
  separate(Site, sep = "_", into = c("Site", "CI")) 

tidy_env_obs_ci$CI<- tidy_env_obs_ci$CI %>% replace_na("Bin_Count")
tidy_env_obs_ci$Site<-sub("p","",tidy_env_obs_ci$Site)

tidy_env_obs_ci<- tidy_env_obs_ci %>% pivot_wider(names_from = CI,  values_from = Count)

tidy_env_obs_ci$Site <- as.factor(tidy_env_obs_ci$Site)
tidy_env_obs_ci$Site <- factor(tidy_env_obs_ci$Site, levels = c(1,12,2,3,4,5,6,7,8,9,10,11))



############
ENV_Plot<- tidy_env_obs_ci %>% filter(env=="A MAT")
median_pop_env_plot<-median_pop_env %>% filter(ENV=="MAT") %>% select(-ENV)


env1_hist <- ggplot(ENV_Plot,aes(x=S,y=Bin_Count,ymin=low,ymax=up))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection (MAT)", y = "Number of SNPs") +
  scale_y_continuous(limits=c(0,50))+ 
  theme_ci() + facet_wrap(.~Site)+
  geom_vline(data = median_pop_env_plot, aes(xintercept = median), linetype="dashed")

env1_hist 
ggsave("graphs/histograms/env/slope_ci_env1.pdf",width=12, height = 8, units = "in")

############
ENV_Plot<- tidy_env_obs_ci %>% filter(env=="B MAP")
median_pop_env_plot<-median_pop_env %>% filter(ENV=="MAP") %>% select(-ENV)


env2_hist <- ggplot(ENV_Plot,aes(x=S,y=Bin_Count,ymin=low,ymax=up))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection (MAP)", y = "Number of SNPs") +
  scale_y_continuous(limits=c(0,50))+ 
  theme_ci() + facet_wrap(.~Site)+
  geom_vline(data = median_pop_env_plot, aes(xintercept = median), linetype="dashed")

env2_hist 
ggsave("graphs/histograms/env/slope_ci_env2.pdf",width=12, height = 8, units = "in")

############
ENV_Plot<- tidy_env_obs_ci %>% filter(env=="C PAS")
median_pop_env_plot<-median_pop_env %>% filter(ENV=="PAS") %>% select(-ENV)


env3_hist <- ggplot(ENV_Plot,aes(x=S,y=Bin_Count,ymin=low,ymax=up))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection (PAS)", y = "Number of SNPs") +
  scale_y_continuous(limits=c(0,50))+ 
  theme_ci() + facet_wrap(.~Site)+
  geom_vline(data = median_pop_env_plot, aes(xintercept = median), linetype="dashed")

env3_hist 
ggsave("graphs/histograms/env/slope_ci_env3.pdf",width=12, height = 8, units = "in")

############
ENV_Plot<- tidy_env_obs_ci %>% filter(env=="D EXT")
median_pop_env_plot<-median_pop_env %>% filter(ENV=="EXT") %>% select(-ENV)


env4_hist <- ggplot(ENV_Plot,aes(x=S,y=Bin_Count,ymin=low,ymax=up))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection (EXT)", y = "Number of SNPs") +
  scale_y_continuous(limits=c(0,50))+ 
  theme_ci() + facet_wrap(.~Site)+
  geom_vline(data = median_pop_env_plot, aes(xintercept = median), linetype="dashed")

env4_hist 
ggsave("graphs/histograms/env/slope_ci_env4.pdf",width=12, height = 8, units = "in")

############
ENV_Plot<- tidy_env_obs_ci %>% filter(env=="E CMD")
median_pop_env_plot<-median_pop_env %>% filter(ENV=="CMD") %>% select(-ENV)


env5_hist <- ggplot(ENV_Plot,aes(x=S,y=Bin_Count,ymin=low,ymax=up))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection (CMD)", y = "Number of SNPs") +
  scale_y_continuous(limits=c(0,50))+ 
  theme_ci() + facet_wrap(.~Site)+
  geom_vline(data = median_pop_env_plot, aes(xintercept = median), linetype="dashed")

env5_hist 
ggsave("graphs/histograms/env/slope_ci_env5.pdf",width=12, height = 8, units = "in")

############
ENV_Plot<- tidy_env_obs_ci %>% filter(env=="F Tave_wt")
median_pop_env_plot<-median_pop_env %>% filter(ENV=="Tave_wt") %>% select(-ENV)


env6_hist <- ggplot(ENV_Plot,aes(x=S,y=Bin_Count,ymin=low,ymax=up))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection (Tave_wt)", y = "Number of SNPs") +
  scale_y_continuous(limits=c(0,50))+ 
  theme_ci() + facet_wrap(.~Site)+
  geom_vline(data = median_pop_env_plot, aes(xintercept = median), linetype="dashed")

env6_hist 
ggsave("graphs/histograms/env/slope_ci_env6.pdf",width=12, height = 8, units = "in")

############
ENV_Plot<- tidy_env_obs_ci %>% filter(env=="G Tave_sm")
median_pop_env_plot<-median_pop_env %>% filter(ENV=="Tave_sm") %>% select(-ENV)


env7_hist <- ggplot(ENV_Plot,aes(x=S,y=Bin_Count,ymin=low,ymax=up))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection (Tave_sm)", y = "Number of SNPs") +
  scale_y_continuous(limits=c(0,50))+ 
  theme_ci() + facet_wrap(.~Site)+
  geom_vline(data = median_pop_env_plot, aes(xintercept = median), linetype="dashed")

env7_hist 
ggsave("graphs/histograms/env/slope_ci_env7.pdf",width=12, height = 8, units = "in")

############
ENV_Plot<- tidy_env_obs_ci %>% filter(env=="H PPT_wt")
median_pop_env_plot<-median_pop_env %>% filter(ENV=="PPT_wt") %>% select(-ENV)


env8_hist <- ggplot(ENV_Plot,aes(x=S,y=Bin_Count,ymin=low,ymax=up))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection (PPT_wt)", y = "Number of SNPs") +
  scale_y_continuous(limits=c(0,50))+ 
  theme_ci() + facet_wrap(.~Site)+
  geom_vline(data = median_pop_env_plot, aes(xintercept = median), linetype="dashed")

env8_hist 
ggsave("graphs/histograms/env/slope_ci_env8.pdf",width=12, height = 8, units = "in")

############
ENV_Plot<- tidy_env_obs_ci %>% filter(env=="I PPT_sm")
median_pop_env_plot<-median_pop_env %>% filter(ENV=="PPT_sm") %>% select(-ENV)


env9_hist <- ggplot(ENV_Plot,aes(x=S,y=Bin_Count,ymin=low,ymax=up))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection (PPT_sm)", y = "Number of SNPs") +
  scale_y_continuous(limits=c(0,50))+ 
  theme_ci() + facet_wrap(.~Site)+
  geom_vline(data = median_pop_env_plot, aes(xintercept = median), linetype="dashed")

env9_hist 
ggsave("graphs/histograms/env/slope_ci_env9.pdf",width=12, height = 8, units = "in")


