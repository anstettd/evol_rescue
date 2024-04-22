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

#Import data
offset_pop <- read_csv("data/binomial_data/time_cumul_beagle.csv")
offset_pop_no4 <- offset_pop %>% filter(Paper_ID!=4)


lm.cumul.1215_pos <- lm(cumul_pos~offset_1215,data=offset_pop)
lm.cumul.1215_all <- lm(cumul_all~offset_1215,data=offset_pop)

lm.cumul.ssp245_pos <- lm(cumul_pos~offset_SSP245,data=offset_pop)
lm.cumul.ssp245_all <- lm(cumul_all~offset_SSP245,data=offset_pop)

lm.cumul.ssp585_pos <- lm(cumul_pos~offset_SSP585,data=offset_pop)
lm.cumul.ssp585_all <- lm(cumul_all~offset_SSP585,data=offset_pop)

lm.cumul.clim_pos <- lm(cumul_pos~offset_climate,data=offset_pop)
lm.cumul.clim_all <- lm(cumul_all~offset_climate,data=offset_pop)

#lm.cumul.1215_pos_no4 <- lm(cumul_pos~offset_1215,data=offset_pop_no4)
#lm.cumul.1215_all_no4 <- lm(cumul_all~offset_1215,data=offset_pop_no4)


Anova(lm.cumul.1215_pos,type="III")
Anova(lm.cumul.1215_all,type="III")
summary(lm.cumul.1215_pos)
summary(lm.cumul.1215_all)

Anova(lm.cumul.ssp245_pos,type="III")
Anova(lm.cumul.ssp245_all,type="III")
summary(lm.cumul.ssp245_pos)
summary(lm.cumul.ssp245_all)


Anova(lm.cumul.ssp585_pos,type="III")
Anova(lm.cumul.ssp585_all,type="III")
summary(lm.cumul.ssp585_pos)
summary(lm.cumul.ssp585_all)

Anova(lm.cumul.clim_pos,type="III")
Anova(lm.cumul.clim_all,type="III")
summary(lm.cumul.clim_pos)
summary(lm.cumul.clim_all)


#Anova(lm.cumul.1215_pos_no4,type="III")
#Anova(lm.cumul.1215_all_no4,type="III")
#summary(lm.cumul.1215_pos_no4)
#summary(lm.cumul.1215_all_no4)




###########################################################################################################

# N-S color gradient
lat_cols=colorRampPalette(brewer.pal(11,"Spectral"))
n.sites <- length(unique(offset_pop$Paper_ID))
color.list <- lat_cols(n.sites)


#cumul slope plotted against 1215 offset
ggplot(offset_pop, aes(x=offset_1215, y=cumul_pos)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Cumulative Positive Selection")+
  scale_x_continuous(name="2012-2015 Genetic Offset")+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
ggsave("graphs/bi_plots_gradient/1_rapid_offset_1215_pos.pdf",width=8, height = 6, units = "in")


#cumul slope plotted against 1215 offset
ggplot(offset_pop, aes(x=offset_1215, y=cumul_all)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Cumulative Selection")+
  scale_x_continuous(name="2012-2015 Genetic Offset")+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
ggsave("graphs/bi_plots_gradient/2_rapid_offset_1215_all.pdf",width=8, height = 6, units = "in")


###########################################################################################################


#cumul slope plotted against 2040-2070 SSP245 Genetic Offset
ggplot(offset_pop, aes(x=offset_SSP245, y=cumul_pos)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Cumulative Positive Selection")+
  scale_x_continuous(name="2040-2070 SSP245 Genetic Offset")+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
ggsave("graphs/bi_plots_gradient/3_rapid_offset_ssp245_pos.pdf",width=8, height = 6, units = "in")

#cumul slope plotted against 2040-2070 SSP245 Genetic Offset
ggplot(offset_pop, aes(x=offset_SSP245, y=cumul_all)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Cumulative Selection")+
  scale_x_continuous(name="2040-2070 SSP245 Genetic Offset")+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
ggsave("graphs/bi_plots_gradient/4_rapid_offset_ssp245_all.pdf",width=8, height = 6, units = "in")


#cumul slope plotted against 2040-2070 SSP585 Genetic Offset
ggplot(offset_pop, aes(x=offset_SSP585, y=cumul_pos)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Cumulative Positive Selection")+
  scale_x_continuous(name="2040-2070 SSP585 Genetic Offset")+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
ggsave("graphs/bi_plots_gradient/5_rapid_offset_ssp585_pos.pdf",width=8, height = 6, units = "in")


#cumul slope plotted against 2040-2070 SSP585 Genetic Offset
ggplot(offset_pop, aes(x=offset_SSP585, y=cumul_all)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Cumulative Selection")+
  scale_x_continuous(name="2040-2070 SSP585 Genetic Offset")+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
ggsave("graphs/bi_plots_gradient/6_rapid_offset_ssp585_all.pdf",width=8, height = 6, units = "in")


###########################################################################################################

#cumul slope plotted against 2012-2015 Climate Distance
ggplot(offset_pop, aes(x=offset_climate, y=cumul_pos)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Cumulative Positive Selection")+
  scale_x_continuous(name="2012-2015 Climate Distance")+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
ggsave("graphs/bi_plots_gradient/7_rapid_offset_climate_distance_all.pdf",width=8, height = 6, units = "in")


#cumul slope plotted against 2012-2015 Climate Distance
ggplot(offset_pop, aes(x=offset_climate, y=cumul_all)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Cumulative Selection")+
  scale_x_continuous(name="2012-2015 Climate Distance")+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines")  # Increase the size of the legend dots
  )
ggsave("graphs/bi_plots_gradient/8_rapid_offset_climate_distance_pos.pdf",width=8, height = 6, units = "in")



###########################################################################################################


#cumul slope plotted against 1215 offset
ggplot(offset_pop_no4, aes(x=offset_1215, y=cumul_pos)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Cumulative Positive Selection")+
  scale_x_continuous(name="2012-2015 Genetic Offset")+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank())
#ggsave("graphs/bi_plots_gradient/9_rapid_offset_1215_pos.pdf",width=8, height = 6, units = "in")


#cumul slope plotted against 1215 offset
ggplot(offset_pop_no4, aes(x=offset_1215, y=cumul_all)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Cumulative Selection")+
  scale_x_continuous(name="2012-2015 Genetic Offset")+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank())
#ggsave("graphs/bi_plots_gradient/10_rapid_offset_1215_all.pdf",width=8, height = 6, units = "in")



