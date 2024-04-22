##################################################################################
## Plot SNP frequency for all baseline populations
## Done for all SNP set that are in the timeseries
## 
## 
## Last Modified Sept 16, 2023
###################################################################################
#Import libraries
library(tidyverse)
library(egg)


###################################################################################

#Baseline Slopes
baseline_snp_slope <- read_csv("data/baseline_SNP_slope.csv")
colnames(baseline_snp_slope) <- c("ENV","chr_snp","Slope","Inv_Logit_Coeff","SE_base")



#Set env variables as factors and set order
baseline_snp_slope$ENV <- as.factor(baseline_snp_slope$ENV)
baseline_snp_slope$ENV <- factor(baseline_snp_slope$ENV, levels = c("MAT","MAP","PAS","EXT","CMD",
                                                    "Tave_wt","Tave_sm","PPT_wt","PPT_sm"))
#Make separate data frames for each climate variable
MAT<-baseline_snp_slope %>% filter(ENV=="MAT")
MAP<-baseline_snp_slope %>% filter(ENV=="MAP")
PAS<-baseline_snp_slope %>% filter(ENV=="PAS")
EXT<-baseline_snp_slope %>% filter(ENV=="EXT")
CMD<-baseline_snp_slope %>% filter(ENV=="CMD")
Tave_wt<-baseline_snp_slope %>% filter(ENV=="Tave_wt")
Tave_sm<-baseline_snp_slope %>% filter(ENV=="Tave_sm")
PPT_wt<-baseline_snp_slope %>% filter(ENV=="PPT_wt")
PPT_sm<-baseline_snp_slope %>% filter(ENV=="PPT_sm")


#Plot frequency vs climate
plot_env_1<- ggplot(MAT,aes(x=Slope))+
    geom_histogram(position="identity", color="black",alpha=0.5, binwidth = 0.1)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="MAT Slope")+
  theme_classic()

#Plot frequency vs climate
plot_env_2<- ggplot(MAP,aes(x=Slope))+
  geom_histogram(position="identity", color="black",alpha=0.5, binwidth = 0.1)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="MAP Slope")+
  theme_classic()

#Plot frequency vs climate
plot_env_3<- ggplot(PAS,aes(x=Slope))+
  geom_histogram(position="identity", color="black",alpha=0.5, binwidth = 0.1)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="PAS Slope")+
  theme_classic()

#Plot frequency vs climate
plot_env_4<- ggplot(EXT,aes(x=Slope))+
  geom_histogram(position="identity", color="black",alpha=0.5, binwidth = 0.1)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="EXT Slope")+
  theme_classic()

#Plot frequency vs climate
plot_env_5<- ggplot(CMD,aes(x=Slope))+
  geom_histogram(position="identity", color="black",alpha=0.5, binwidth = 0.1)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="CMD Slope")+
  theme_classic()

#Plot frequency vs climate
plot_env_6<- ggplot(Tave_wt,aes(x=Slope))+
  geom_histogram(position="identity", color="black",alpha=0.5, binwidth = 0.1)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Tave_wt Slope")+
  theme_classic()

#Plot frequency vs climate
plot_env_7<- ggplot(Tave_sm,aes(x=Slope))+
  geom_histogram(position="identity", color="black",alpha=0.5, binwidth = 0.1)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="Tave_sm Slope")+
  theme_classic()

plot_env_8<- ggplot(PPT_wt,aes(x=Slope))+
  geom_histogram(position="identity", color="black",alpha=0.5, binwidth = 0.1)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="PPT_wt Slope")+
  theme_classic()

plot_env_9<- ggplot(PPT_sm,aes(x=Slope))+
  geom_histogram(position="identity", color="black",alpha=0.5, binwidth = 0.1)+
  scale_y_continuous(name="Count")+
  scale_x_continuous(name="PPT_sm Slope")+
  theme_classic()


#Make single plot of all env variables vs SNP frequency
arrg_1<-ggarrange(plot_env_1, 
                  plot_env_2,
                  plot_env_3,                  
                  plot_env_4, 
                  plot_env_5, 
                  plot_env_6,
                  plot_env_7, 
                  plot_env_8, 
                  plot_env_9 , 
                  nrow = 3, ncol = 3)
arrg_1

#Export 8 X 7
ggsave("graphs/histograms/baseline_histogram/1_baseline_histogram.pdf",arrg_1,width=10, height = 9, units = "in")




