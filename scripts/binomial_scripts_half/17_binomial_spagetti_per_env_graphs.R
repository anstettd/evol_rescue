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


binomial_table<-function(freq_A, freq_B, env_in){
  
  freq_A.melted <- reshape2::melt(freq_A, id.vars = c("Site", "Year"))
  colnames(freq_A.melted)<-c("Site", "Year", "SNP_ID", "Count_A")
  freq_B.melted <- reshape2::melt(freq_B, id.vars = c("Site", "Year"))
  colnames(freq_B.melted)<-c("Site", "Year", "SNP_ID", "Count_B")
  freq_table<-left_join(freq_A.melted, freq_B.melted, by=c("Site", "SNP_ID", "Year"))
  
  freq_table$Site <- as.numeric(freq_table$Site)
  freq_table <- left_join(freq_table,env_in,by=c("Site","SNP_ID")) %>% 
    filter(SE<5.5) %>% 
    select(-Slope, -ENV, -SE, -Type)
  
  
  freq_table_Final<-data.frame()
  
  for (i in 1:dim(freq_table)[1]){
    
    Binomial_A<-c(rep(1, freq_table[i, "Count_A"]), rep(0,freq_table[i, "Count_B"]))
    
    tmp_df<- as.data.frame(Binomial_A ) %>% mutate (Site=freq_table$Site[i], 
                                                    SNP_ID=freq_table$SNP_ID[i], 
                                                    Year=freq_table$Year[i])
    
    freq_table_Final<-rbind(freq_table_Final, tmp_df)
    
   # print(i)
  }
  
  
  freq_table_Final$Site <- as.factor(freq_table_Final$Site)
  freq_table_Final$Site <- factor(freq_table_Final$Site, levels = c(1,12,2,3,4,5,6,7,8,9,10,11))
  
  
  return(freq_table_Final)
  
}


############

#Import timeseries frequencies
freqA_env1 <- read_csv("data/binomial_data/freqA_env1.csv")
freqA_env2 <- read_csv("data/binomial_data/freqA_env2.csv")
freqA_env3 <- read_csv("data/binomial_data/freqA_env3.csv")
freqA_env4 <- read_csv("data/binomial_data/freqA_env4.csv")
freqA_env5 <- read_csv("data/binomial_data/freqA_env5.csv")
freqA_env6 <- read_csv("data/binomial_data/freqA_env6.csv")
freqA_env7 <- read_csv("data/binomial_data/freqA_env7.csv")
freqA_env8 <- read_csv("data/binomial_data/freqA_env8.csv")
freqA_env9 <- read_csv("data/binomial_data/freqA_env9.csv")

freqB_env1 <- read_csv("data/binomial_data/freqB_env1.csv")
freqB_env2 <- read_csv("data/binomial_data/freqB_env2.csv")
freqB_env3 <- read_csv("data/binomial_data/freqB_env3.csv")
freqB_env4 <- read_csv("data/binomial_data/freqB_env4.csv")
freqB_env5 <- read_csv("data/binomial_data/freqB_env5.csv")
freqB_env6 <- read_csv("data/binomial_data/freqB_env6.csv")
freqB_env7 <- read_csv("data/binomial_data/freqB_env7.csv")
freqB_env8 <- read_csv("data/binomial_data/freqB_env8.csv")
freqB_env9 <- read_csv("data/binomial_data/freqB_env9.csv")

#SE data
env_all <- read_csv("data/binomial_data_half/slope_obs_all_unique.csv")
colnames(env_all) <- c("Site","SNP_ID","Slope","SE","ENV","Type")

freq_table_1<-binomial_table(freqA_env1, freqB_env1,env_all)
freq_table_2<-binomial_table(freqA_env2, freqB_env2,env_all)
freq_table_3<-binomial_table(freqA_env3, freqB_env3,env_all)
freq_table_4<-binomial_table(freqA_env4, freqB_env4,env_all)
freq_table_5<-binomial_table(freqA_env5, freqB_env5,env_all)
freq_table_6<-binomial_table(freqA_env6, freqB_env6,env_all)
freq_table_7<-binomial_table(freqA_env7, freqB_env7,env_all)
freq_table_8<-binomial_table(freqA_env8, freqB_env8,env_all)
freq_table_9<-binomial_table(freqA_env9, freqB_env9,env_all)


#All Sites
ggplot(data=freq_table_1 ,aes(Year,Binomial_A,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~Site) + theme_spaghetti()
ggsave("graphs/spaghetii/spaghetii_obs_env1.pdf",width=8, height = 7, units = "in")

#All Sites
ggplot(data=freq_table_2 ,aes(Year,Binomial_A,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~Site) + theme_spaghetti()
ggsave("graphs/spaghetii/spaghetii_obs_env2.pdf",width=8, height = 7, units = "in")

#All Sites
ggplot(data=freq_table_3 ,aes(Year,Binomial_A,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~Site) + theme_spaghetti()
ggsave("graphs/spaghetii/spaghetii_obs_env3.pdf",width=8, height = 7, units = "in")

#All Sites
ggplot(data=freq_table_4 ,aes(Year,Binomial_A,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~Site) + theme_spaghetti()
ggsave("graphs/spaghetii/spaghetii_obs_env4.pdf",width=8, height = 7, units = "in")

#All Sites
ggplot(data=freq_table_5 ,aes(Year,Binomial_A,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~Site) + theme_spaghetti()
ggsave("graphs/spaghetii/spaghetii_obs_env5.pdf",width=8, height = 7, units = "in")

#All Sites
ggplot(data=freq_table_6 ,aes(Year,Binomial_A,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~Site) + theme_spaghetti()
ggsave("graphs/spaghetii/spaghetii_obs_env6.pdf",width=8, height = 7, units = "in")

#All Sites
ggplot(data=freq_table_7 ,aes(Year,Binomial_A,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~Site) + theme_spaghetti()
ggsave("graphs/spaghetii/spaghetii_obs_env7.pdf",width=8, height = 7, units = "in")

#All Sites
ggplot(data=freq_table_8 ,aes(Year,Binomial_A,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~Site) + theme_spaghetti()
ggsave("graphs/spaghetii/spaghetii_obs_env8.pdf",width=8, height = 7, units = "in")

#All Sites
ggplot(data=freq_table_9 ,aes(Year,Binomial_A,group=SNP_ID)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.09,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="Year") +  facet_wrap(.~Site) + theme_spaghetti()
ggsave("graphs/spaghetii/spaghetii_obs_env9.pdf",width=8, height = 7, units = "in")


