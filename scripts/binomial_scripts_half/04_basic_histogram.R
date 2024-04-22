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
env_all <- read_csv("data/binomial_data/slope_obs_all_unique.csv")
env_low <- env_all %>% filter(SE<5.5) 
#%>% filter(Slope<=5) %>% filter(Slope>=-5)

env_1 <- env_low %>% filter(Site==1)
env_2 <- env_low %>% filter(Site==2)
env_3 <- env_low %>% filter(Site==3)
env_4 <- env_low %>% filter(Site==4)
env_5 <- env_low %>% filter(Site==5)
env_6 <- env_low %>% filter(Site==6)
env_7 <- env_low %>% filter(Site==7)
env_8 <- env_low %>% filter(Site==8)
env_9 <- env_low %>% filter(Site==9)
env_10 <- env_low %>% filter(Site==10)
env_11 <- env_low %>% filter(Site==11)
env_12 <- env_low %>% filter(Site==12)






###################################################################################
## Low SE only, all pop, all env
###################################################################################

#all_slopes
env_slope<- ggplot(env_low,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection")+
  theme_classic()
#env_slope <- env_slope  + theme(
#  axis.text.x = element_text(size=12, face="bold"),
#  axis.text.y = element_text(size=12,face="bold"),
#  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
#  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
env_slope <- env_slope + facet_wrap(.~Site)
env_slope
ggsave("graphs/slope_obs/obs_slope_ab.pdf",width=10, height = 7.5, units = "in")



#all_SE
env_SE<- ggplot(env_low,aes(x=SE))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="SE")+
  theme_classic()
#env_SE <- env_SE  + theme(
#  axis.text.x = element_text(size=12, face="bold"),
#  axis.text.y = element_text(size=12,face="bold"),
#  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
# axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
env_SE <- env_SE + facet_wrap(.~ENV)
env_SE
ggsave("graphs/slope_obs/low_SE.pdf",width=10, height = 7.5, units = "in")

#SE versus Slope
env_SE<- ggplot(env_low,aes(x=Slope,y=SE))+
  geom_point(shape=4,size=0.8)+
  scale_y_continuous(name="SE")+
  scale_x_continuous(name="Strength of Selection")+
  theme_classic()
#env_SE <- env_SE  + theme(
#  axis.text.x = element_text(size=12, face="bold"),
#  axis.text.y = element_text(size=12,face="bold"),
#  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
# axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
env_SE <- env_SE + facet_wrap(.~ENV)
env_SE
ggsave("graphs/slope_obs/low_SE_vs_slope.pdf",width=10, height = 7.5, units = "in")






###################################################################################
## All high SE included, all pop, all env
###################################################################################

#all_slopes
env_slope<- ggplot(env_all,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection")+
  theme_classic()
#env_slope <- env_slope  + theme(
#  axis.text.x = element_text(size=12, face="bold"),
#  axis.text.y = element_text(size=12,face="bold"),
#  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
#  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
env_slope <- env_slope + facet_wrap(.~Site)
env_slope
ggsave("graphs/slope_obs/all_slope_high.pdf",width=10, height = 7.5, units = "in")

#all_SE
env_SE<- ggplot(env_all,aes(x=SE))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection")+
  theme_classic()
#env_SE <- env_SE  + theme(
#  axis.text.x = element_text(size=12, face="bold"),
#  axis.text.y = element_text(size=12,face="bold"),
#  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
# axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
env_SE <- env_SE + facet_wrap(.~Site)
env_SE
ggsave("graphs/slope_obs/all_SE_high.pdf",width=10, height = 7.5, units = "in")

#SE versus Slope
env_SE<- ggplot(env_all,aes(x=Slope,y=SE))+
  geom_point(shape=4,size=0.8)+
  scale_y_continuous(name="SE")+
  scale_x_continuous(name="Strength of Selection")+
  theme_classic()
#env_SE <- env_SE  + theme(
#  axis.text.x = element_text(size=12, face="bold"),
#  axis.text.y = element_text(size=12,face="bold"),
#  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
# axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
env_SE <- env_SE + facet_wrap(.~Site)
env_SE
ggsave("graphs/slope_obs/all_SE_vs_slope_high.pdf",width=10, height = 7.5, units = "in")










###################################################################################
## Full graphs per population
###################################################################################

#p1 
env_hist <- ggplot(env_1,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 1)",limits=c(-5,5))+
  theme_classic()
env_hist <- env_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
env_hist <- env_hist + facet_wrap(.~ENV)
env_hist
ggsave("graphs/slope_obs/01_slope_obs_p1.pdf",width=10, height = 7.5, units = "in")

#p12 
env_hist <- ggplot(env_12,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 12)",limits=c(-5,5))+
  theme_classic()
env_hist <- env_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
env_hist <- env_hist + facet_wrap(.~ENV)
env_hist
ggsave("graphs/slope_obs/02_slope_obs_p12.pdf",width=10, height = 7.5, units = "in")

#p2 
env_hist <- ggplot(env_2,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 2)",limits=c(-5,5))+
  theme_classic()
env_hist <- env_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
env_hist <- env_hist + facet_wrap(.~ENV)
env_hist
ggsave("graphs/slope_obs/03_slope_obs_p2.pdf",width=10, height = 7.5, units = "in")

#p3 
env_hist <- ggplot(env_3,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 3)",limits=c(-5,5))+
  theme_classic()
env_hist <- env_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
env_hist <- env_hist + facet_wrap(.~ENV)
env_hist
ggsave("graphs/slope_obs/04_slope_obs_p3.pdf",width=10, height = 7.5, units = "in")

#p4 
env_hist <- ggplot(env_4,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 4)",limits=c(-5,5))+
  theme_classic()
env_hist <- env_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
env_hist <- env_hist + facet_wrap(.~ENV)
env_hist
ggsave("graphs/slope_obs/05_slope_obs_p4.pdf",width=10, height = 7.5, units = "in")

#p5 
env_hist <- ggplot(env_5,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 5)",limits=c(-5,5))+
  theme_classic()
env_hist <- env_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
env_hist <- env_hist + facet_wrap(.~ENV)
env_hist
ggsave("graphs/slope_obs/06_slope_obs_p5.pdf",width=10, height = 7.5, units = "in")

#p6 
env_hist <- ggplot(env_6,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 6)",limits=c(-5,5))+
  theme_classic()
env_hist <- env_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
env_hist <- env_hist + facet_wrap(.~ENV)
env_hist
ggsave("graphs/slope_obs/07_slope_obs_p6.pdf",width=10, height = 7.5, units = "in")

#p7 
env_hist <- ggplot(env_7,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 7)",limits=c(-5,5))+
  theme_classic()
env_hist <- env_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
env_hist <- env_hist + facet_wrap(.~ENV)
env_hist
ggsave("graphs/slope_obs/08_slope_obs_p7.pdf",width=10, height = 7.5, units = "in")

#p8 
env_hist <- ggplot(env_8,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 8)",limits=c(-5,5))+
  theme_classic()
env_hist <- env_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
env_hist <- env_hist + facet_wrap(.~ENV)
env_hist
ggsave("graphs/slope_obs/09_slope_obs_p8.pdf",width=10, height = 7.5, units = "in")

#p9 
env_hist <- ggplot(env_9,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 9)",limits=c(-5,5))+
  theme_classic()
env_hist <- env_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
env_hist <- env_hist + facet_wrap(.~ENV)
env_hist
ggsave("graphs/slope_obs/10_slope_obs_p9.pdf",width=10, height = 7.5, units = "in")

#p10 
env_hist <- ggplot(env_10,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 10)",limits=c(-5,5))+
  theme_classic()
env_hist <- env_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
env_hist <- env_hist + facet_wrap(.~ENV)
env_hist
ggsave("graphs/slope_obs/11_slope_obs_p10.pdf",width=10, height = 7.5, units = "in")

#p11 
env_hist <- ggplot(env_11,aes(x=Slope))+
  geom_histogram(position="identity",bins=40,color="black",alpha=0.5)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Strength of Selection (Pop 11)",limits=c(-5,5))+
  theme_classic()
env_hist <- env_hist  + theme(
  axis.text.x = element_text(size=12, face="bold"),
  axis.text.y = element_text(size=12,face="bold"),
  axis.title.x = element_text(color="black", size=14, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=14,vjust = 2, face="bold",hjust=0.5))
env_hist <- env_hist + facet_wrap(.~ENV)
env_hist
ggsave("graphs/slope_obs/12_slope_obs_p11.pdf",width=10, height = 7.5, units = "in")




