##################################################################################
## Plot cumulative slope against demography
## Author Daniel Anstett
## 
## 
## Last Modified April 19, 2023
###################################################################################

#Library install and import
library(tidyverse)
library(RColorBrewer)
library(cowplot)

#Import data
offset_pop <- read_csv("data/binomial_data/time_cumul_beagle.csv")


###########################################################################################################

# N-S color gradient
lat_cols=colorRampPalette(brewer.pal(11,"Spectral"))
n.sites <- length(unique(offset_pop$Paper_ID))
color.list <- lat_cols(n.sites)


#cumul slope plotted against 1215 offset
plot_1 <- ggplot(offset_pop, aes(x=offset_1215, y=cumul_pos)) + 
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


#cumul slope plotted against 1215 offset
plot_2 <-ggplot(offset_pop, aes(x=offset_1215, y=cumul_all)) + 
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
    legend.title = element_blank(), legend.position="none")


###########################################################################################################


#cumul slope plotted against 2040-2070 SSP245 Genetic Offset
plot_3 <-ggplot(offset_pop, aes(x=offset_SSP245, y=cumul_pos)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Cumulative Positive Selection")+
  scale_x_continuous(name="SSP245 Genetic Offset")+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=18,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(), legend.position="none")


#cumul slope plotted against 2040-2070 SSP245 Genetic Offset
plot_4 <-ggplot(offset_pop, aes(x=offset_SSP245, y=cumul_all)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Cumulative Selection")+
  scale_x_continuous(name="SSP245 Genetic Offset")+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank())


#cumul slope plotted against 2040-2070 SSP585 Genetic Offset
plot_5 <-ggplot(offset_pop, aes(x=offset_SSP585, y=cumul_pos)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Cumulative Positive Selection")+
  scale_x_continuous(name="SSP585 Genetic Offset")+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=18,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(), legend.position="none")


#cumul slope plotted against 2040-2070 SSP585 Genetic Offset
plot_6 <-ggplot(offset_pop, aes(x=offset_SSP585, y=cumul_all)) + 
  geom_point(aes(fill=as.factor(round(Lat, 1))),shape=21,size =4.5)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Cumulative Selection")+
  scale_x_continuous(name="SSP585 Genetic Offset")+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5),
    legend.title = element_blank(), legend.position="none")


#cumul slope plotted against 2012-2015 Climate Distance
plot_7 <- ggplot(offset_pop, aes(x=offset_climate, y=cumul_pos)) + 
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
    legend.title = element_blank())



#cumul slope plotted against 2012-2015 Climate Distance
plot_8 <-ggplot(offset_pop, aes(x=offset_climate, y=cumul_all)) + 
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
    legend.title = element_blank())


###########################################################################################################
###########################################################################################################
#Make Cowplots

plot_grid(plot_2,plot_3,plot_4,plot_5,plot_6, labels = "AUTO",ncol = 3,label_x = 0.23) #export at 6 X 8 
























