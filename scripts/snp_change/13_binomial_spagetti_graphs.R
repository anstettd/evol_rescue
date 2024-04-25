################################################################################################################
## Regression plots of BF>20 SNPs for random and observed
## Author Julia Anstett
## 
## 
## Last Modified Aug 16, 2023
################################################################################################################
#Import libraries
library(tidyverse)
library(reshape2)

################################################################################################################
#Functions

#Define Theme

theme_spaghetti <- function(){ 
  
  theme_classic() %+replace%    #replace elements we want to change
    theme(legend.position = "none",
          axis.text.x = element_text(size = 14, face = "bold", angle = 45,hjust = 1, vjust = 1), 
          axis.title = element_text(size =16, face = "bold"), 
          axis.text.y = element_text(size = 14, face = "bold"),legend.title = element_blank(),
          legend.text = element_text(size=12,face="bold"),
          strip.background = element_blank(), 
          strip.text.x=element_text(size=14,face="bold",hjust=0,vjust=-1.2))
}

theme( 
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0, face="bold"),
  axis.title.y = element_text(color="black", size=20, vjust = 1.6,face="bold",hjust=0.5 ,angle=90)
)

############

#Import timeseries frequencies
freqA_env1 <- read_csv("data/binomial_strong/freqA_env1.csv")
freqA_env2 <- read_csv("data/binomial_strong/freqA_env2.csv")
freqA_env3 <- read_csv("data/binomial_strong/freqA_env3.csv")
freqA_env4 <- read_csv("data/binomial_strong/freqA_env4.csv")
freqA_env5 <- read_csv("data/binomial_strong/freqA_env5.csv")
freqA_env6 <- read_csv("data/binomial_strong/freqA_env6.csv")
freqA_env7 <- read_csv("data/binomial_strong/freqA_env7.csv")
freqA_env8 <- read_csv("data/binomial_strong/freqA_env8.csv")
freqA_env9 <- read_csv("data/binomial_strong/freqA_env9.csv")

freqB_env1 <- read_csv("data/binomial_strong/freqB_env1.csv")
freqB_env2 <- read_csv("data/binomial_strong/freqB_env2.csv")
freqB_env3 <- read_csv("data/binomial_strong/freqB_env3.csv")
freqB_env4 <- read_csv("data/binomial_strong/freqB_env4.csv")
freqB_env5 <- read_csv("data/binomial_strong/freqB_env5.csv")
freqB_env6 <- read_csv("data/binomial_strong/freqB_env6.csv")
freqB_env7 <- read_csv("data/binomial_strong/freqB_env7.csv")
freqB_env8 <- read_csv("data/binomial_strong/freqB_env8.csv")
freqB_env9 <- read_csv("data/binomial_strong/freqB_env9.csv")



freq_1A.melted <- reshape2::melt(freqA_env1, id.vars = c("Site", "Year"))
colnames(freq_1A.melted)<-c("Site", "Year", "SNP_ID", "Count_A")
freq_1B.melted <- reshape2::melt(freqB_env1, id.vars = c("Site", "Year"))
colnames(freq_1B.melted)<-c("Site", "Year", "SNP_ID", "Count_B")
freq_table_1<-left_join(freq_1A.melted, freq_1B.melted, by=c("Site", "SNP_ID", "Year")) %>% mutate (env="env1")


freq_2A.melted <- reshape2::melt(freqA_env2, id.vars = c("Site", "Year"))
colnames(freq_2A.melted)<-c("Site", "Year", "SNP_ID", "Count_A")
freq_2B.melted <- reshape2::melt(freqB_env2, id.vars = c("Site", "Year"))
colnames(freq_2B.melted)<-c("Site", "Year", "SNP_ID", "Count_B")
freq_table_2<-left_join(freq_2A.melted, freq_2B.melted, by=c("Site", "SNP_ID", "Year")) %>% mutate (env="env2")


freq_3A.melted <- reshape2::melt(freqA_env3, id.vars = c("Site", "Year"))
colnames(freq_3A.melted)<-c("Site", "Year", "SNP_ID", "Count_A")
freq_3B.melted <- reshape2::melt(freqB_env3, id.vars = c("Site", "Year"))
colnames(freq_3B.melted)<-c("Site", "Year", "SNP_ID", "Count_B")
freq_table_3<-left_join(freq_3A.melted, freq_3B.melted, by=c("Site", "SNP_ID", "Year")) %>% mutate (env="env3")



freq_4A.melted <- reshape2::melt(freqA_env4, id.vars = c("Site", "Year"))
colnames(freq_4A.melted)<-c("Site", "Year", "SNP_ID", "Count_A")
freq_4B.melted <- reshape2::melt(freqB_env4, id.vars = c("Site", "Year"))
colnames(freq_4B.melted)<-c("Site", "Year", "SNP_ID", "Count_B")
freq_table_4<-left_join(freq_4A.melted, freq_4B.melted, by=c("Site", "SNP_ID", "Year")) %>% mutate (env="env4")


freq_5A.melted <- reshape2::melt(freqA_env5, id.vars = c("Site", "Year"))
colnames(freq_5A.melted)<-c("Site", "Year", "SNP_ID", "Count_A")
freq_5B.melted <- reshape2::melt(freqB_env5, id.vars = c("Site", "Year"))
colnames(freq_5B.melted)<-c("Site", "Year", "SNP_ID", "Count_B")
freq_table_5<-left_join(freq_5A.melted, freq_5B.melted, by=c("Site", "SNP_ID", "Year")) %>% mutate (env="env5")



freq_6A.melted <- reshape2::melt(freqA_env6, id.vars = c("Site", "Year"))
colnames(freq_6A.melted)<-c("Site", "Year", "SNP_ID", "Count_A")
freq_6B.melted <- reshape2::melt(freqB_env6, id.vars = c("Site", "Year"))
colnames(freq_6B.melted)<-c("Site", "Year", "SNP_ID", "Count_B")
freq_table_6<-left_join(freq_6A.melted, freq_6B.melted, by=c("Site", "SNP_ID", "Year")) %>% mutate (env="env6")


freq_7A.melted <- reshape2::melt(freqA_env7, id.vars = c("Site", "Year"))
colnames(freq_7A.melted)<-c("Site", "Year", "SNP_ID", "Count_A")
freq_7B.melted <- reshape2::melt(freqB_env7, id.vars = c("Site", "Year"))
colnames(freq_7B.melted)<-c("Site", "Year", "SNP_ID", "Count_B")
freq_table_7<-left_join(freq_7A.melted, freq_7B.melted, by=c("Site", "SNP_ID", "Year")) %>% mutate (env="env7")


freq_8A.melted <- reshape2::melt(freqA_env8, id.vars = c("Site", "Year"))
colnames(freq_8A.melted)<-c("Site", "Year", "SNP_ID", "Count_A")
freq_8B.melted <- reshape2::melt(freqB_env8, id.vars = c("Site", "Year"))
colnames(freq_8B.melted)<-c("Site", "Year", "SNP_ID", "Count_B")
freq_table_8<-left_join(freq_8A.melted, freq_8B.melted, by=c("Site", "SNP_ID", "Year")) %>% mutate (env="env8")


freq_9A.melted <- reshape2::melt(freqA_env9, id.vars = c("Site", "Year"))
colnames(freq_9A.melted)<-c("Site", "Year", "SNP_ID", "Count_A")
freq_9B.melted <- reshape2::melt(freqB_env9, id.vars = c("Site", "Year"))
colnames(freq_9B.melted)<-c("Site", "Year", "SNP_ID", "Count_B")
freq_table_9<-left_join(freq_9A.melted, freq_9B.melted, by=c("Site", "SNP_ID", "Year")) %>% mutate (env="env9")




freq_table_all<-rbind(freq_table_1, 
                      freq_table_2, 
                      freq_table_3,  
                      freq_table_4,  
                      freq_table_5,  
                      freq_table_6, 
                      freq_table_7,  
                      freq_table_8,  
                      freq_table_9) %>% unite(col = "Filter_ID", c("Site", "Year","SNP_ID"), sep = "/") %>%
  distinct(Filter_ID, .keep_all=T) %>% separate(Filter_ID, c("Site", "Year","SNP_ID"), sep = "/") %>% select(-env)


#SE data
env_all <- read_csv("data/binomial_strong/slope_obs_all_unique.csv")
colnames(env_all) <- c("Site","SNP_ID","Slope","SE","ENV","Type")
#env_low <- env_all %>% filter(SE<5.5) 
freq_table_all$Site <- as.numeric(freq_table_all$Site)
freq_table_all <- left_join(freq_table_all,env_all,by=c("Site","SNP_ID")) %>% 
  filter(SE<5.5) %>% 
  select(-Slope, -ENV, -SE, -Type)


freq_table_Final<-data.frame()

for (i in 1:dim(freq_table_all)[1]){
  
  Binomial_A<-c(rep(1, freq_table_all[i, "Count_A"]), rep(0,freq_table_all[i, "Count_B"]))

  tmp_df<- as.data.frame(Binomial_A ) %>% mutate (Site=freq_table_all$Site[i], 
                                                           SNP_ID=freq_table_all$SNP_ID[i], 
                                                           Year=freq_table_all$Year[i])
  
  freq_table_Final<-rbind(freq_table_Final, tmp_df)
  
  print(i)
}


freq_table_Final$Site <- as.factor(freq_table_Final$Site)
freq_table_Final$Site <- factor(freq_table_Final$Site, levels = c(1,12,2,3,4,5,6,7,8,9,10,11))




#All Sites
ggplot(data=freq_table_Final ,aes(Year,Binomial_A,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~Site) + theme_spaghetti()
ggsave("graphs/spaghetii/Strong/spaghetii_obs.pdf",width=8, height = 7, units = "in")










#Combination sites
freq_env_ver1 <- freq_table_Final %>% filter(Site==2 | Site==3 | Site==11)
ggplot(data=freq_env_ver1 ,aes(Year,Binomial_A,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~Site) + theme_spaghetti() +  
  theme(strip.text.x = element_text(size=0),axis.title = element_text(size =18, face = "bold"),
        axis.title.x = element_blank(), axis.text.x = element_text(size = 16, face = "bold", angle = 45,hjust = 1, vjust = 1))
ggsave("graphs/spaghetii/1_3_11spaghetii_obs.pdf",width=14, height = 5, units = "in")

freq_env_ver2 <- freq_table_Final %>% filter(Site==1 | Site==3 | Site==6 | Site==11)
ggplot(data=freq_env_ver2 ,aes(Year,Binomial_A,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~Site, ncol = 4) + theme_spaghetti()
#ggsave("graphs/spaghetii/1_3_6_11spaghetii_obs.pdf",width=11, height = 4, units = "in")


#Selected Sites
freq_env_1 <- freq_table_Final %>% filter(Site==1)
ggplot(data=freq_env_1 ,aes(Year,Binomial_A,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="Year") + theme_spaghetti()
#ggsave("graphs/spaghetii/01_spaghetii_obs.pdf",width=4, height = 4, units = "in")


freq_env_3 <- freq_table_Final %>% filter(Site==3)
ggplot(data=freq_env_3 ,aes(Year,Binomial_A,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="Year") + theme_spaghetti()
#ggsave("graphs/spaghetii/03_spaghetii_obs.pdf",width=4, height = 4, units = "in")


freq_env_6 <- freq_table_Final %>% filter(Site==6)
ggplot(data=freq_env_6 ,aes(Year,Binomial_A,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="Year") + theme_spaghetti()
#ggsave("graphs/spaghetii/06_spaghetii_obs.pdf",width=4, height = 4, units = "in")


freq_env_11 <- freq_table_Final %>% filter(Site==11)
ggplot(data=freq_env_11 ,aes(Year,Binomial_A,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="Year") + theme_spaghetti()
#ggsave("graphs/spaghetii/11_spaghetii_obs.pdf",width=4, height = 4, units = "in")
