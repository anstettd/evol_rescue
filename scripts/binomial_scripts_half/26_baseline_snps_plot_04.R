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
theme_spaghetti <- function(){ 
  
  theme_classic() %+replace%    #replace elements we want to change
    theme(legend.position = "none",
          axis.text.x = element_text(size = 14, face = "bold", angle = 0,hjust = 1, vjust = 1), 
          axis.title = element_text(size =16, face = "bold"), 
          axis.text.y = element_text(size = 14, face = "bold"),legend.title = element_blank(),
          legend.text = element_text(size=12,face="bold"),
          strip.background = element_blank(), 
          strip.text.x=element_text(size=14,face="bold",hjust=0,vjust=-1.2))
}


###################################################################################
#Import SNP binomial data
abund_clim <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/abund_table_baseline_slope_SE_std_04.csv",
                       col_names = T)

###################################################################################

#Set env variables as factors and set order
abund_clim$ENV <- as.factor(abund_clim$ENV)
abund_clim$ENV <- factor(abund_clim$ENV, levels = c("MAT","MAP","PAS","EXT","CMD",
                                                     "Tave_wt","Tave_sm","PPT_wt","PPT_sm"))
#Make separate data frames for each climate variable
abund_MAT<-abund_clim %>% filter(ENV=="MAT")
abund_MAP<-abund_clim %>% filter(ENV=="MAP")
abund_PAS<-abund_clim %>% filter(ENV=="PAS")
abund_EXT<-abund_clim %>% filter(ENV=="EXT")
abund_CMD<-abund_clim %>% filter(ENV=="CMD")
abund_Tave_wt<-abund_clim %>% filter(ENV=="Tave_wt")
abund_Tave_sm<-abund_clim %>% filter(ENV=="Tave_sm")
abund_PPT_wt<-abund_clim %>% filter(ENV=="PPT_wt")
abund_PPT_sm<-abund_clim %>% filter(ENV=="PPT_sm")

############################################################################################################################
## Latitude
############################################################################################################################
#Plot frequency vs latitude vs
ggplot(data=abund_clim,aes(Latitude,Binomial_A,group=chr_snp)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.15,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="Latitude") +  facet_wrap(.~ENV) + theme_spaghetti()
ggsave("graphs/spaghetii/base_slope_cutoff/7_03_spaghetii_baseline_lat_04.pdf",width=8, height = 7, units = "in")

############################################################################################################################
## Climate
############################################################################################################################
#Plot frequency vs climate
plot_env_1<-ggplot(data=abund_MAT,aes(MAT,Binomial_A,group=chr_snp)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.1,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="MAT")  + theme_spaghetti()

plot_env_2<-ggplot(data=abund_MAP,aes(MAP,Binomial_A,group=chr_snp)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.1,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="MAP")  + theme_spaghetti()

plot_env_3<-ggplot(data=abund_PAS,aes(PAS,Binomial_A,group=chr_snp)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.1,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="PAS")  + theme_spaghetti()

plot_env_4<-ggplot(data=abund_EXT,aes(EXT,Binomial_A,group=chr_snp)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.1,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="EXT")  + theme_spaghetti()

plot_env_5<-ggplot(data=abund_CMD,aes(CMD,Binomial_A,group=chr_snp)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.1,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="CMD")  + theme_spaghetti()

plot_env_6<-ggplot(data=abund_Tave_wt,aes(Tave_wt,Binomial_A,group=chr_snp)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.1,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="Tave_wt")  + theme_spaghetti()

plot_env_7<-ggplot(data=abund_Tave_sm,aes(Tave_sm,Binomial_A,group=chr_snp)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.1,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="Tave_sm")  + theme_spaghetti()

plot_env_8<-ggplot(data=abund_PPT_wt,aes(PPT_wt,Binomial_A,group=chr_snp)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.1,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="PPT_wt")  + theme_spaghetti()

plot_env_9<-ggplot(data=abund_PPT_sm,aes(PPT_sm,Binomial_A,group=chr_snp)) + 
  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.1,cex=0.4,color="blue") + 
  labs(y="SNP Frequency",x="PPT_sm")  + theme_spaghetti()


#Make single plot of all env variables vs SNP frequency
arrg_1<-ggarrange(plot_env_1, 
                  plot_env_2 +theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_blank()),
                  plot_env_3 +theme(axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.title.y = element_blank()),                  
                  plot_env_4, 
                  plot_env_5 +theme(axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.title.y = element_blank()), 
                  plot_env_6 +theme(axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.title.y = element_blank()),
                  plot_env_7, 
                  plot_env_8 +theme(axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.title.y = element_blank()), 
                  plot_env_9 +theme(axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.title.y = element_blank()), 
                  nrow = 3, ncol = 3)

#Export 8 X 7

ggsave("graphs/spaghetii/base_slope_cutoff/8_04_spaghetii_baseline_env_std.pdf",arrg_1,width=8, height = 7, units = "in")







