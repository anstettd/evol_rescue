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
#Import data
offset_pop <- read_csv("data/binomial_data/time_cumul_beagle.csv")
#offset_pop_3 <- offset_pop %>% filter(Paper_ID!=11)


lm.cumul.1215_pos <- lm(cumul_pos~offset_1215,data=offset_pop)
lm.cumul.1215_all <- lm(cumul_all~offset_1215,data=offset_pop)

lm.cumul.ssp245_pos <- lm(cumul_pos~offset_SSP245,data=offset_pop)
lm.cumul.ssp245_all <- lm(cumul_all~offset_SSP245,data=offset_pop)

lm.cumul.ssp585_pos <- lm(cumul_pos~offset_SSP585,data=offset_pop)
lm.cumul.ssp585_all <- lm(cumul_all~offset_SSP585,data=offset_pop)

lm.cumul.clim_pos <- lm(cumul_pos~offset_climate,data=offset_pop)
lm.cumul.clim_all <- lm(cumul_all~offset_climate,data=offset_pop)

Anova(lm.cumul.1215_pos,type="III")
Anova(lm.cumul.1215_all,type="III")
summary(lm.cumul.1215_pos)
summary(lm.cumul.1215_all)

Anova(lm.cumul.ssp245_pos,type="III")
Anova(lm.cumul.ssp245_all,type="III")

Anova(lm.cumul.ssp585_pos,type="III")
Anova(lm.cumul.ssp585_all,type="III")

Anova(lm.cumul.clim_pos,type="III")
Anova(lm.cumul.clim_all,type="III")




#cumul slope plotted against 1215 offset
ggplot(offset_pop, aes(x=offset_1215, y=cumul_pos, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Cumulative Positive Selection")+
  scale_x_continuous(name="2012-2015 Genetic Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.65, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("graphs/bi_plots/1_rapid_offset_1215_pos.pdf",width=7, height = 6, units = "in")


#cumul slope plotted against 1215 offset
ggplot(offset_pop, aes(x=offset_1215, y=cumul_all, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Cumulative Selection")+
  scale_x_continuous(name="2012-2015 Genetic Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.65, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("graphs/bi_plots/2_rapid_offset_1215_all.pdf",width=7, height = 6, units = "in")


###########################################################################################################


#cumul slope plotted against 2040-2070 SSP245 Genetic Offset
ggplot(offset_pop, aes(x=offset_SSP245, y=cumul_pos, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Cumulative Positive Selection")+
  scale_x_continuous(name="2040-2070 SSP245 Genetic Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.65, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("graphs/bi_plots/3_rapid_offset_ssp245_pos.pdf",width=7, height = 6, units = "in")

#cumul slope plotted against 2040-2070 SSP245 Genetic Offset
ggplot(offset_pop, aes(x=offset_SSP245, y=cumul_all, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Cumulative Selection")+
  scale_x_continuous(name="2040-2070 SSP245 Genetic Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.65, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("graphs/bi_plots/4_rapid_offset_ssp245_all.pdf",width=7, height = 6, units = "in")


#cumul slope plotted against 2040-2070 SSP585 Genetic Offset
ggplot(offset_pop, aes(x=offset_SSP585, y=cumul_pos, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Cumulative Positive Selection")+
  scale_x_continuous(name="2040-2070 SSP585 Genetic Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.65, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("graphs/bi_plots/5_rapid_offset_ssp585_pos.pdf",width=7, height = 6, units = "in")


#cumul slope plotted against 2040-2070 SSP585 Genetic Offset
ggplot(offset_pop, aes(x=offset_SSP585, y=cumul_all, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Cumulative Selection")+
  scale_x_continuous(name="2040-2070 SSP585 Genetic Offset")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.65, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("graphs/bi_plots/6_rapid_offset_ssp585_all.pdf",width=7, height = 6, units = "in")


###########################################################################################################

#cumul slope plotted against 2012-2015 Climate Distance
ggplot(offset_pop, aes(x=offset_climate, y=cumul_pos, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Cumulative Positive Selection")+
  scale_x_continuous(name="2012-2015 Climate Distance")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.65, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("graphs/bi_plots/7_rapid_offset_climate_distance_all.pdf",width=7, height = 6, units = "in")


#cumul slope plotted against 2012-2015 Climate Distance
ggplot(offset_pop, aes(x=offset_climate, y=cumul_all, label=Paper_ID)) + 
  geom_point(aes(color=Region), size =4.5)+
  geom_smooth(method=lm,color="black")+
  #  geom_label_repel(aes(label = ID))+
  geom_text(hjust=-.15, vjust=-.2)+
  scale_y_continuous(name="Cumulative Selection")+
  scale_x_continuous(name="2012-2015 Climate Distance")+
  scale_color_manual(values= c("North"="#3399FF", "Center"="#FFCC00", "South"="#FF3333"))+
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.position = c(0.65, 0.85),legend.text=element_text(size=14),
    legend.title=element_text(size=16,face="bold"))
ggsave("graphs/bi_plots/8_rapid_offset_climate_distance_pos.pdf",width=7, height = 6, units = "in")





