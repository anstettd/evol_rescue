##################################################################################
## Plot SNP frequency for all baseline populations
## Done for all SNP set that are in the timeseries
## 
## This script exports the "strong" snp set.
## Has only SNPs with a linear correlation with climate variables
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
abund_clim <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/abund_table_baseline_slope_SE_std_03.csv",
                       col_names = T)

###################################################################################

#Set env variables as factors and set order
abund_clim$ENV <- as.factor(abund_clim$ENV)
abund_clim$ENV <- factor(abund_clim$ENV, levels = c("MAT","MAP","PAS","EXT","CMD",
                                                     "Tave_wt","Tave_sm","PPT_wt","PPT_sm"))
#Make separate data frames for each climate variable
abund_env1<-abund_clim %>% filter(ENV=="MAT")
abund_env2<-abund_clim %>% filter(ENV=="MAP")
abund_env3<-abund_clim %>% filter(ENV=="PAS")
abund_env4<-abund_clim %>% filter(ENV=="EXT")
abund_env5<-abund_clim %>% filter(ENV=="CMD")
abund_env6 <-abund_clim %>% filter(ENV=="Tave_wt")
abund_env7 <-abund_clim %>% filter(ENV=="Tave_sm")
abund_env8 <-abund_clim %>% filter(ENV=="PPT_wt")
abund_env9 <-abund_clim %>% filter(ENV=="PPT_sm")

#Filter for regressions close to 0 or 1 
abund_env1_rm <- abund_env1 %>% group_by(chr_snp) %>% summarize(avg_slope = mean(Binomial_A)) %>% 
  filter(avg_slope<0.97) %>% filter(avg_slope>0.03) 
abund_env2_rm <- abund_env2 %>% group_by(chr_snp) %>% summarize(avg_slope = mean(Binomial_A)) %>% 
  filter(avg_slope<0.97) %>% filter(avg_slope>0.03) 
abund_env3_rm <- abund_env3 %>% group_by(chr_snp) %>% summarize(avg_slope = mean(Binomial_A)) %>% 
  filter(avg_slope<0.97) %>% filter(avg_slope>0.03) 
abund_env4_rm <- abund_env4 %>% group_by(chr_snp) %>% summarize(avg_slope = mean(Binomial_A)) %>% 
  filter(avg_slope<0.97) %>% filter(avg_slope>0.03) 
abund_env5_rm <- abund_env5 %>% group_by(chr_snp) %>% summarize(avg_slope = mean(Binomial_A)) %>% 
  filter(avg_slope<0.97) %>% filter(avg_slope>0.03) 
abund_env6_rm <- abund_env6 %>% group_by(chr_snp) %>% summarize(avg_slope = mean(Binomial_A)) %>% 
  filter(avg_slope<0.97) %>% filter(avg_slope>0.03) 
abund_env7_rm <- abund_env7 %>% group_by(chr_snp) %>% summarize(avg_slope = mean(Binomial_A)) %>% 
  filter(avg_slope<0.97) %>% filter(avg_slope>0.03) 
abund_env8_rm <- abund_env8 %>% group_by(chr_snp) %>% summarize(avg_slope = mean(Binomial_A)) %>% 
  filter(avg_slope<0.97) %>% filter(avg_slope>0.03) 
abund_env9_rm <- abund_env9 %>% group_by(chr_snp) %>% summarize(avg_slope = mean(Binomial_A)) %>% 
  filter(avg_slope<0.97) %>% filter(avg_slope>0.03) 

#Make SNP list that has all weak basline slopes removed
strong_snp_set <- rbind(abund_env1_rm[,1],
                        abund_env2_rm[,1],
                        abund_env3_rm[,1],
                        abund_env4_rm[,1],
                        abund_env5_rm[,1],
                        abund_env6_rm[,1],
                        abund_env7_rm[,1],
                        abund_env8_rm[,1],
                        abund_env9_rm[,1]) %>% unique()
write_csv(abund_env1_rm[,1],"data/binomial_strong/strong_snp_set_env1.csv")
write_csv(abund_env2_rm[,1],"data/binomial_strong/strong_snp_set_env2.csv")
write_csv(abund_env3_rm[,1],"data/binomial_strong/strong_snp_set_env3.csv")
write_csv(abund_env4_rm[,1],"data/binomial_strong/strong_snp_set_env4.csv")
write_csv(abund_env5_rm[,1],"data/binomial_strong/strong_snp_set_env5.csv")
write_csv(abund_env6_rm[,1],"data/binomial_strong/strong_snp_set_env6.csv")
write_csv(abund_env7_rm[,1],"data/binomial_strong/strong_snp_set_env7.csv")
write_csv(abund_env8_rm[,1],"data/binomial_strong/strong_snp_set_env8.csv")
write_csv(abund_env9_rm[,1],"data/binomial_strong/strong_snp_set_env9.csv")
write_csv(strong_snp_set,"data/binomial_strong/strong_snp_set.csv")


#Filter SNPs close to 0 or 1
abund_MAT<-abund_env1 %>% filter(chr_snp %in% abund_env1_rm$chr_snp)
abund_MAP<-abund_env2 %>% filter(chr_snp %in% abund_env2_rm$chr_snp)
abund_PAS<-abund_env3 %>% filter(chr_snp %in% abund_env3_rm$chr_snp)
abund_EXT<-abund_env4 %>% filter(chr_snp %in% abund_env4_rm$chr_snp)
abund_CMD<-abund_env5 %>% filter(chr_snp %in% abund_env5_rm$chr_snp)
abund_Tave_wt<-abund_env6 %>% filter(chr_snp %in% abund_env6_rm$chr_snp)
abund_Tave_sm<-abund_env7 %>% filter(chr_snp %in% abund_env7_rm$chr_snp)
abund_PPT_wt<-abund_env8 %>% filter(chr_snp %in% abund_env8_rm$chr_snp)
abund_PPT_sm<-abund_env9 %>% filter(chr_snp %in% abund_env9_rm$chr_snp)


############################################################################################################################
## Latitude
############################################################################################################################
##Plot frequency vs latitude vs
#ggplot(data=abund_clim,aes(Latitude,Binomial_A,group=chr_snp)) + 
#  geom_line(stat="smooth",method = "glm", method.args = list(family = "binomial"), se = F, alpha=.15,cex=0.4,color="blue") + 
#  labs(y="SNP Frequency",x="Latitude") +  facet_wrap(.~ENV) + theme_spaghetti()
#ggsave("graphs/spaghetii/base_slope_cutoff/7_03_spaghetii_baseline_lat_03.pdf",width=8, height = 7, units = "in")

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

ggsave("graphs/spaghetii/base_slope_cutoff/8_03_spaghetii_baseline_env_03_extreme_cuttoff.pdf",arrg_1,width=8, height = 7, units = "in")







