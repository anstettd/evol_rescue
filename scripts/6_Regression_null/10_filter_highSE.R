##########################################################################################################
## Remove high SE values for all snp glm in each 100 k region
## 
## Author Daniel Anstett
## 
## 
## Last Modified May 3, 2023
##########################################################################################################
#Import libraries
library(tidyverse)

##########################################################################################################
#Import and filter for high SE

swiss_glm_1 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/swiss_glm_50_50_1.csv")
swiss_glm_2 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/swiss_glm_50_50_2.csv")
swiss_glm_3 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/swiss_glm_50_50_3.csv")
swiss_glm_4 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/swiss_glm_50_50_4.csv")
swiss_glm_5 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/swiss_glm_50_50_5.csv")
swiss_glm_6 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/swiss_glm_50_50_6.csv")
swiss_glm_7 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/swiss_glm_50_50_7.csv")
swiss_glm_8 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/swiss_glm_50_50_8.csv")
swiss_glm_9 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/swiss_glm_50_50_9.csv")
swiss_glm_10 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/swiss_glm_50_50_10.csv")
swiss_glm_11 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/swiss_glm_50_50_11.csv")
swiss_glm_12 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/swiss_glm_50_50_12.csv")
swiss_glm_13 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/swiss_glm_50_50_13.csv")
swiss_glm_14 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/swiss_glm_50_50_14.csv")
swiss_glm_15 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/swiss_glm_50_50_15.csv")
swiss_glm_16 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/swiss_glm_50_50_16.csv")
swiss_glm_17 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/swiss_glm_50_50_17.csv")


#Bind glms
swiss_glm <- rbind(swiss_glm_1,
                   swiss_glm_2,
                   swiss_glm_3,
                   swiss_glm_4,
                   swiss_glm_5,
                   swiss_glm_6,
                   swiss_glm_7,
                   swiss_glm_8,
                   swiss_glm_9,
                   swiss_glm_10,
                   swiss_glm_11,
                   swiss_glm_12,
                   swiss_glm_13,
                   swiss_glm_14,
                   swiss_glm_15,
                   swiss_glm_16,
                   swiss_glm_17)

#Filter for SE < 5
swiss_glm_filter <- swiss_glm %>% filter(SE<5)

#Export
write_csv(swiss_glm_filter,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/swiss_glm_50_50.csv")
#write_csv(swiss_glm,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/swiss_glm_50_50_all.csv")

#write_csv(swiss_glm_filter,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_filter_minor.csv")
#write_csv(swiss_glm,"/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/swiss_glm_ci_allSE_minor.csv")


#Discriptive Plots

#all_slopes
env_slope<- ggplot(swiss_glm_filter,aes(x=Slope))+
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


