##################################################################################
## Make strenth of selection histograms
## Author Daniel Anstett
## 
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
env_obs_ci_unique <- read_csv("data/snp_change_data/obs_ci_env_unique.csv")

#Import Medians
median_pop <- read_csv("data/snp_change_data/median_pop.csv")


#Isolate each pop and lable
env_p1 <- env_obs_ci_unique %>% select(S,p1,p1_low,p1_up) %>% mutate(Site=1,pop_lable="A Site 1")
env_p12 <- env_obs_ci_unique %>% select(S,p12,p1_low,p12_up) %>% mutate(Site=12,pop_lable="B Site 12")
env_p2 <- env_obs_ci_unique %>% select(S,p2,p1_low,p2_up) %>% mutate(Site=2,pop_lable="C Site 2")
env_p3 <- env_obs_ci_unique %>% select(S,p3,p1_low,p3_up) %>% mutate(Site=3,pop_lable="D Site 3")
env_p4 <- env_obs_ci_unique %>% select(S,p4,p1_low,p4_up) %>% mutate(Site=4,pop_lable="E Site 4")
env_p5 <- env_obs_ci_unique %>% select(S,p5,p1_low,p5_up) %>% mutate(Site=5,pop_lable="F Site 5")
env_p6 <- env_obs_ci_unique %>% select(S,p6,p1_low,p6_up) %>% mutate(Site=6,pop_lable="G Site 6")
env_p7 <- env_obs_ci_unique %>% select(S,p7,p1_low,p7_up) %>% mutate(Site=7,pop_lable="H Site 7")
env_p8 <- env_obs_ci_unique %>% select(S,p8,p1_low,p8_up) %>% mutate(Site=8,pop_lable="I Site 8")
env_p9 <- env_obs_ci_unique %>% select(S,p9,p1_low,p9_up) %>% mutate(Site=9,pop_lable="J Site 9")
env_p10 <- env_obs_ci_unique %>% select(S,p10,p1_low,p10_up) %>% mutate(Site=10,pop_lable="K Site 10")
env_p11 <- env_obs_ci_unique %>% select(S,p11,p1_low,p11_up) %>% mutate(Site=11,pop_lable="L Site 11")

#Rename columns
colnames(env_p1) = c("S","obs","low","high","Site","pop_lable")
colnames(env_p2) = c("S","obs","low","high","Site","pop_lable")
colnames(env_p3) = c("S","obs","low","high","Site","pop_lable")
colnames(env_p4) = c("S","obs","low","high","Site","pop_lable")
colnames(env_p5) = c("S","obs","low","high","Site","pop_lable")
colnames(env_p6) = c("S","obs","low","high","Site","pop_lable")
colnames(env_p7) = c("S","obs","low","high","Site","pop_lable")
colnames(env_p8) = c("S","obs","low","high","Site","pop_lable")
colnames(env_p9) = c("S","obs","low","high","Site","pop_lable")
colnames(env_p10) = c("S","obs","low","high","Site","pop_lable")
colnames(env_p11) = c("S","obs","low","high","Site","pop_lable")
colnames(env_p12) = c("S","obs","low","high","Site","pop_lable")

#Bind long dataset
env_histPop <- rbind(env_p1,
                     env_p12,
                     env_p2,
                     env_p3,
                     env_p4,
                     env_p5,
                     env_p6,
                     env_p7,
                     env_p8,
                     env_p9,
                     env_p10,
                     env_p11)

env_histPop_25 <- env_histPop %>% filter(S <= 1.25 & S>= -1.25) 


site_unique <- env_histPop %>% select(Site,pop_lable) 
site_unique <- unique(site_unique)

median_pop <- left_join(median_pop,site_unique, by="Site")

env_histPop_1 <- env_histPop_25 %>% filter(Site==3 | Site==11)
env_histPop_2 <- env_histPop_25 %>% filter(Site==2 | Site==3 | Site==11)



###################################################################################
## Histogram graphs with CI
###################################################################################

# -2.5 to 2.5
histPop <- ggplot(env_histPop ,aes(x=S,y=obs,ymin=low,ymax=high))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  #geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection", y = "Number of SNPs") +
  #scale_y_continuous(limits=c(0,40))+ 
  theme_ci() + facet_wrap(.~pop_lable) +
  geom_vline(data = median_pop, aes(xintercept = median), linetype="dashed")
histPop 
ggsave("Graphs/Selection/01_selection_2.5.pdf", histPop, width=12, height = 8, units = "in")


# -1.25 to 1.25
histPop <- ggplot(env_histPop_25 ,aes(x=S,y=obs,ymin=low,ymax=high))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  #geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection", y = "Number of SNPs") +
  scale_y_continuous(limits=c(0,125),breaks=seq(0,125,by=25))+ 
  theme_ci() + facet_wrap(.~pop_lable)+
  geom_vline(data = median_pop, aes(xintercept = median), linetype="dashed")

histPop
ggsave("Graphs/Selection/02_selection_1.25.pdf", histPop, width=12, height = 8, units = "in")



###################################################################################
## Individual Graphs
###################################################################################

#Sites 3, 11
histPop1 <- ggplot(env_histPop_1 ,aes(x=S,y=obs,ymin=low,ymax=high))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  #geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection", y = "Number of SNPs") +
  scale_y_continuous(breaks=seq(0,100,by=25))+ 
  theme_ci() + facet_wrap(.~pop_lable, ncol = 4)+ theme(strip.text.x = element_text(size=0))+
  geom_vline(data = median_pop_filter_1, aes(xintercept = median), linetype="dashed")

histPop1
#Export 
ggsave("Graphs/Selection/03_s_pop_3_11.pdf",width=10, height = 4, units = "in")


#Sites 2, 3, 11

histPop1 <- ggplot(env_histPop_2 ,aes(x=S,y=obs,ymin=low,ymax=high))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  #geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection", y = "Number of SNPs") +
  scale_y_continuous(breaks=seq(0,100,by=25))+ 
  theme_ci() + facet_wrap(.~pop_lable, ncol = 4) + theme(strip.text.x = element_text(size=0)) +
  geom_vline(data = median_pop_filter_1, aes(xintercept = median), linetype="dashed")

histPop1
#Export 
ggsave("Graphs/Selection/04_s_pop_1_3_11.pdf",width=11, height = 3.5, units = "in")




