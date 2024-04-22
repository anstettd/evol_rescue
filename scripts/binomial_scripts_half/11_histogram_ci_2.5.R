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


###################################################################################
## Histogram graphs with CI
###################################################################################

# p1
mat_p1_hist <- ggplot(env_obs_ci,aes(x=S,y=p1,ymin=p1_low,ymax=p1_up))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection (P1)", y = "Number of SNPs") +
  scale_y_continuous(limits=c(0,40))+ theme_ci() + facet_wrap(.~env)
mat_p1_hist 
ggsave("graphs/slopes_CI_2.5_half_95/01_slope_ci_p1.pdf",width=12, height = 8, units = "in")

# p12
mat_p12_hist <- ggplot(env_obs_ci,aes(x=S,y=p12,ymin=p12_low,ymax=p12_up))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection (p12)", y = "Number of SNPs") +
  scale_y_continuous(limits=c(0,40))+ theme_ci() + facet_wrap(.~env)
mat_p12_hist 
ggsave("graphs/slopes_CI_2.5_half_95/02_slope_ci_p12.pdf",width=12, height = 8, units = "in")


# p2
mat_p2_hist <- ggplot(env_obs_ci,aes(x=S,y=p2,ymin=p2_low,ymax=p2_up))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection (p2)", y = "Number of SNPs") +
  scale_y_continuous(limits=c(0,40))+ theme_ci() + facet_wrap(.~env)
mat_p2_hist 
ggsave("graphs/slopes_CI_2.5_half_95/03_slope_ci_p2.pdf",width=12, height = 8, units = "in")

# p3
mat_p3_hist <- ggplot(env_obs_ci,aes(x=S,y=p3,ymin=p3_low,ymax=p3_up))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection (p3)", y = "Number of SNPs") +
  scale_y_continuous(limits=c(0,40))+ theme_ci() + facet_wrap(.~env)
mat_p3_hist 
ggsave("graphs/slopes_CI_2.5_half_95/04_slope_ci_p3.pdf",width=12, height = 8, units = "in")

# p4
mat_p4_hist <- ggplot(env_obs_ci,aes(x=S,y=p4,ymin=p4_low,ymax=p4_up))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection (p4)", y = "Number of SNPs") +
  scale_y_continuous(limits=c(0,40))+ theme_ci() + facet_wrap(.~env)
mat_p4_hist 
ggsave("graphs/slopes_CI_2.5_half_95/05_slope_ci_p4.pdf",width=12, height = 8, units = "in")

# p5
mat_p5_hist <- ggplot(env_obs_ci,aes(x=S,y=p5,ymin=p5_low,ymax=p5_up))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection (p5)", y = "Number of SNPs") +
  scale_y_continuous(limits=c(0,40))+ theme_ci() + facet_wrap(.~env)
mat_p5_hist 
ggsave("graphs/slopes_CI_2.5_half_95/06_slope_ci_p5.pdf",width=12, height = 8, units = "in")

# p6
mat_p6_hist <- ggplot(env_obs_ci,aes(x=S,y=p6,ymin=p6_low,ymax=p6_up))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection (p6)", y = "Number of SNPs") +
  scale_y_continuous(limits=c(0,40))+ theme_ci() + facet_wrap(.~env)
mat_p6_hist 
ggsave("graphs/slopes_CI_2.5_half_95/07_slope_ci_p6.pdf",width=12, height = 8, units = "in")

# p7
mat_p7_hist <- ggplot(env_obs_ci,aes(x=S,y=p7,ymin=p7_low,ymax=p7_up))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection (p7)", y = "Number of SNPs") +
  scale_y_continuous(limits=c(0,40))+ theme_ci() + facet_wrap(.~env)
mat_p7_hist 
ggsave("graphs/slopes_CI_2.5_half_95/08_slope_ci_p7.pdf",width=12, height = 8, units = "in")

# p8
mat_p8_hist <- ggplot(env_obs_ci,aes(x=S,y=p8,ymin=p8_low,ymax=p8_up))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection (p8)", y = "Number of SNPs") +
  scale_y_continuous(limits=c(0,40))+ theme_ci() + facet_wrap(.~env)
mat_p8_hist 
ggsave("graphs/slopes_CI_2.5_half_95/09_slope_ci_p8.pdf",width=12, height = 8, units = "in")

# p9
mat_p9_hist <- ggplot(env_obs_ci,aes(x=S,y=p9,ymin=p9_low,ymax=p9_up))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection (p9)", y = "Number of SNPs") +
  scale_y_continuous(limits=c(0,40))+ theme_ci() + facet_wrap(.~env)
mat_p9_hist 
ggsave("graphs/slopes_CI_2.5_half_95/10_slope_ci_p9.pdf",width=12, height = 8, units = "in")

# p10
mat_p10_hist <- ggplot(env_obs_ci,aes(x=S,y=p10,ymin=p10_low,ymax=p10_up))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection (p10)", y = "Number of SNPs") +
  scale_y_continuous(limits=c(0,40))+ theme_ci() + facet_wrap(.~env)
mat_p10_hist 
ggsave("graphs/slopes_CI_2.5_half_95/11_slope_ci_p10.pdf",width=12, height = 8, units = "in")

# p11
mat_p11_hist <- ggplot(env_obs_ci,aes(x=S,y=p11,ymin=p11_low,ymax=p11_up))+
  geom_bar(colour = "black", stat = "identity", width = 0.1, fill = "lightblue1")+
  geom_errorbar(colour = "firebrick2", stat = "identity", width = 0.06) +
  geom_vline(xintercept=0) +
  labs(x = "Strength of Selection (p11)", y = "Number of SNPs") +
  scale_y_continuous(limits=c(0,40))+ theme_ci() + facet_wrap(.~env)
mat_p11_hist 
ggsave("graphs/slopes_CI_2.5_half_95/12_slope_ci_p11.pdf",width=12, height = 8, units = "in")






