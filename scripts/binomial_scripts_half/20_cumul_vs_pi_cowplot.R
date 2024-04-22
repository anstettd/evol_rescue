##################################################################################
## Plot cumulative slope against demography
## Author Daniel Anstett
## 
## 
## Last Modified April 19, 2023
###################################################################################

#Library install and import
library(tidyverse)
library(car)
library(RColorBrewer)
library(cowplot)

#Import data
offset_pop <- read_csv("data/binomial_data_half/time_cumul_beagle.csv")
pi.df <- read_csv("data/binomial_data/raw_pi.csv")

#Join and filter out for site 12 which has geographic sampling differneces between 2010 and 2011+
pi_cumul <- left_join(offset_pop,pi.df,by=c("Paper_ID"="Site")) %>% filter(Paper_ID!=12)

#Models
lm1 <- lm(cumul_pos~pi_snp_set,data=pi_cumul)
lm2 <- lm(cumul_pos~pi_all_snps,data=pi_cumul)
lm3 <- lm(cumul_all~pi_snp_set,data=pi_cumul)
lm4 <- lm(cumul_all~pi_all_snps,data=pi_cumul)

#Anova

Anova(lm1,type="III")
summary(lm1)
Anova(lm2,type="III")
summary(lm2)
Anova(lm3,type="III")
summary(lm3)
Anova(lm4,type="III")
summary(lm4)




###########################################################################################################

# N-S color gradient
lat_cols=colorRampPalette(brewer.pal(11,"Spectral"))
n.sites <- length(unique(offset_pop$Paper_ID))
color.list <- lat_cols(n.sites)


#directional selection vs PI Cliamte Associated
plot_1 <- ggplot(pi_cumul, aes(x=pi_snp_set, y=cumul_pos)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Directional Selection")+
  #scale_x_continuous(name="PI Climate Associated",limits=c(0.15,0.35),breaks=c(0.15,0.2,0.25,0.30,0.35))+
  scale_x_continuous(name="PI Climate Associated")+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=18, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=18,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.position = "none"
  )
#ggsave("graphs/bi_plot_pi/1_pi_cumul1.pdf",width=8, height = 6, units = "in")

#direcitonal selection vs PI Genome-Wide
plot_2 <- ggplot(pi_cumul, aes(x=pi_all_snps, y=cumul_pos)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Directional Selection")+
  scale_x_continuous(name="PI Genome-Wide",limits=c(0.13,0.25),breaks=c(0.15,0.2,0.25))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=18, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=18,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.position = "none"
  )
#ggsave("graphs/bi_plot_pi/2_pi_cumul2.pdf",width=8, height = 6, units = "in")



###########################################################################################################

#total selection vs PI Climate Associated
plot_3 <- ggplot(pi_cumul, aes(x=pi_snp_set, y=cumul_all)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Total Selection")+
  #scale_x_continuous(name="PI Climate Associated",limits=c(0.13,0.25),breaks=c(0.15,0.2,0.25))+
  scale_x_continuous(name="PI Climate Associated")+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=18, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=18,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.position = "none"
  )
#ggsave("graphs/bi_plot_pi/3_pi_cumul3.pdf",width=8, height = 6, units = "in")


#total selection vs PI Genome-Wide
plot_4 <- ggplot(pi_cumul, aes(x=pi_all_snps, y=cumul_all)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Total Selection")+
  scale_x_continuous(name="PI Genome-Wide",limits=c(0.13,0.25),breaks=c(0.15,0.2,0.25))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=18, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=18,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.position = "none"
  )
#ggsave("graphs/bi_plot_pi/4_pi_cumul4.pdf",width=8, height = 6, units = "in")



############################################################################################################
#Cowplot
plot_grid(plot_1,plot_2,plot_3,plot_4, labels = "AUTO",ncol = 2,label_x = 0.23) #export at 8 X 9


