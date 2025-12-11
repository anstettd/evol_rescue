#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Produce demography graphs for publication
#### AUTHOR: Amy Angert & Daniel Anstett
#### DATE LAST MODIFIED: 20251211


#*******************************************************************************
#### 0. Clean workspace and load required packages
#*******************************************************************************

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

# Make vector of packages needed
packages_needed <- c("tidyverse", "RColorBrewer","cowplot")

# Install packages needed (if not already installed)
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
}

# Load packages needed
for (i in 1:length(packages_needed)){
  library(packages_needed[i], character.only = TRUE)
}


#*******************************************************************************
### 1. Read in lambda estimates for each time period 
#*******************************************************************************
r_means <- read.csv("data/demography data/siteYear.lambda_responses_2010-2019.csv") %>% 
  dplyr::select(Site,Paper_ID,Latitude,Longitude,Lat.Color,mean.r.pre,mean.r.drought,mean.r.recovery) %>%  
  arrange(Latitude)

r_se <- read.csv("data/demography data/siteYear.lambda_responses_2010-2019.csv") %>% 
  dplyr::select(Site,Paper_ID,Latitude,Longitude,Lat.Color,se.r.pre,se.r.drought,se.r.recovery) %>%  
  arrange(Latitude)

# pivot longer for plotting
r_means_long <- r_means %>% 
  pivot_longer(cols=mean.r.pre:mean.r.recovery, names_to="time", values_to="mean_r")
r_means_long$time <- gsub("mean.r.", "", r_means_long$time)

r_se_long <- r_se %>% 
  pivot_longer(cols=se.r.pre:se.r.recovery, names_to="time", values_to="se_r")
r_se_long$time <- gsub("se.r.", "", r_means_long$time)

r_long <- left_join(r_means_long, r_se_long) %>% 
  mutate(ymin=mean_r-se_r,
         ymax=mean_r+se_r)

# set level orders to reflect time
level_order_all = c("pre", "drought", "recovery")

dodge <- position_dodge(width=0.2)

colours_in<-r_long %>% 
  dplyr::select(Latitude, Lat.Color) 

colours_in_to_plot<-unique(colours_in) #%>% filter(Lat.Color=="#FDB567") #remove this population because it doesn't have drought r value


#*******************************************************************************
### 2. Visualize mean estimates over each period for all sites
#*******************************************************************************

# all three time periods (alternative Fig S1, option 1)
ggplot(r_long, aes(x=factor(time, level=level_order_all), y=mean_r, group=Latitude, fill=as.factor(Latitude))) +
  facet_wrap(~Latitude, scale="free") +
  geom_line(aes(colour=as.factor(Latitude)), position=dodge, size=1.5) +
  geom_point(shape=21, size=3, position=dodge) +
  #geom_errorbar(aes(ymin=ymin, ymax=ymax, colour=as.factor(Latitude)), width=0.1, position=dodge) +
  scale_fill_manual(values=colours_in_to_plot$Lat.Color, 
                     labels = unique(r_long$Paper_ID) ) +
  scale_color_manual(values=colours_in_to_plot$Lat.Color, 
                    labels = unique(r_long$Paper_ID )) +
  #scale_color_manual(values=unique(r_long$Lat.Color), aesthetics=c("color", "fill")) +
  scale_y_continuous(name="Mean Population Growth Rate")+ 
  scale_x_discrete(name="Time Period", labels=c("Pre-drought", "Drought", "Recovery")) + 
  geom_hline(yintercept=0) +
  theme_classic() + theme(
    axis.text.x=element_text(size=16,face="bold"),
    axis.text.y=element_text(size=16,face="bold"),
    axis.title.x=element_text(color="black", size=20, vjust=0.5, face="bold"),
    axis.title.y=element_text(color="black", size=20,vjust=2, face="bold", hjust=0.5),
    strip.background=element_blank(), 
    strip.text.x=element_blank(),
    legend.text = element_text(size = 14),
    legend.key.size = unit(2, "lines"),
    legend.key.height = unit(1.6, "lines"), #Reduce height
    legend.title=element_blank())+
  guides(color = guide_legend(reverse = TRUE, override.aes = list(linetype = 0)),
    fill  = guide_legend(reverse = TRUE))
ggsave("Graphs/Demography/01b_decline_recovery_rmeans_faceted.pdf",width=8, height = 6.5, units = "in")

r_long <- r_long %>% mutate(time_num = ifelse(time=="pre", 1, ifelse(time=="drought", 2, 3)))

U_mod_raw_lin <- lm(mean_r ~ time_num, data=r_long)
U_mod_raw_quad <- lm(mean_r ~ poly(time_num,2), data=r_long)
AIC(U_mod_raw_lin, U_mod_raw_quad)

#Pre-drought to drought period
dodge <- position_dodge(width=0.2)
r_long_early <- r_long %>% filter(time!="recovery") %>% droplevels()
level_order_pre = c("pre", "drought")

colours_in<-r_long_early %>% 
  dplyr::select(Latitude, Lat.Color) 

colours_in_to_plot<-unique(colours_in) #%>% filter(Lat.Color!="#FDB567") #remove this population because it doesn't have drought r value


a <- ggplot(r_long_early, aes(x=factor(time, level=level_order_pre), y=mean_r, group=Latitude, fill=as.factor(Latitude))) +
  geom_point(shape=21, size=4, position=dodge) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax, colour=as.factor(Latitude)), width=0.1, position=dodge) +
  geom_line(aes(colour=as.factor(Latitude)), position=dodge, size=1.5) +
  #scale_color_manual(values=colours_in_to_plot$Lat.Color, 
  #                   labels=round(colours_in_to_plot$Latitude,1),aesthetics=c("color", "fill")) +
  scale_fill_manual(values=colours_in_to_plot$Lat.Color, 
                    labels = unique(r_long_early$Paper_ID) ) +
  scale_color_manual(values=colours_in_to_plot$Lat.Color, 
                     labels = unique(r_long_early$Paper_ID )) +
  scale_y_continuous(name="Mean Population Growth Rate")+ 
  scale_x_discrete(name="Time Period", labels=c("Pre-drought", "Drought")) + 
  geom_hline(yintercept=0) +
  theme_classic() + theme(
    axis.text.x = element_text(size=16, face="bold"),
    axis.text.y = element_text(size=16, face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust=0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust=2, face="bold", hjust=0.5),
    strip.background=element_blank(), 
    strip.text.x=element_blank(),
    legend.title=element_blank())+
    guides(color = guide_legend(reverse = TRUE, override.aes = list(linetype = 0)),
       fill  = guide_legend(reverse = TRUE))
a
ggsave("Graphs/Demography/01c_decline.pdf",a, width=6, height = 8, units = "in")

#Drought to recovery period
r_long_late <- r_long %>% filter(time!="pre") %>% droplevels()
level_order_post = c("drought", "recovery")

colours_in<-r_long_late %>% 
  dplyr::select(Latitude, Lat.Color) 

colours_in_to_plot<-unique(colours_in) #%>% filter(Lat.Color!="#FDB567") #remove this population because it doesn't have drought r value


b <- ggplot(r_long_late, aes(x=factor(time, level=level_order_post), y=mean_r, group=Latitude, fill=as.factor(Latitude))) +
  geom_point(shape=21, size=4, position=dodge) +
  geom_line(aes(colour=as.factor(Latitude)), position=dodge, size=1.5) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax, colour=as.factor(Latitude)), width=0.1, position=dodge) +
  scale_y_continuous(name="Mean Population Growth Rate")+ 
  scale_x_discrete(name="Time Period", labels=c("Drought", "Recovery")) + 
  scale_fill_manual(values=colours_in_to_plot$Lat.Color, 
                    labels = unique(r_long_early$Paper_ID) ) +
  scale_color_manual(values=colours_in_to_plot$Lat.Color, 
                     labels = unique(r_long_early$Paper_ID )) +
  geom_hline(yintercept=0) +
  theme_classic() + theme(
    axis.text.x = element_text(size=16, face="bold"),
    axis.text.y = element_text(size=16,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust=0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24, vjust=2, face="bold", hjust=0.5),
    strip.background = element_blank(), 
    strip.text.x = element_blank(),
    legend.title = element_blank())+
  guides(color = guide_legend(reverse = TRUE, override.aes = list(linetype = 0)),
         fill  = guide_legend(reverse = TRUE))
b
ggsave("Graphs/Demography/01d_recovery.pdf",b, width=6, height = 8, units = "in")


# sensitivity test: remove 3 extirpated populations
r_long_late_cull <- filter(r_long_late, mean_r>-2)

colours_in<-r_long_late_cull %>% 
  dplyr::select(Latitude, Lat.Color) 

colours_in_to_plot<-unique(colours_in) #%>% filter(Lat.Color!="#FDB567") #remove this population because it doesn't have drought r value



ggplot(r_long_late_cull, aes(x=factor(time, level=level_order_post), y=mean_r, group=Latitude, fill=as.factor(Latitude))) +
  geom_point(shape=21, size=3, position=dodge) +
  geom_line(aes(colour=as.factor(Latitude)), position=dodge, size=1.5) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax, colour=as.factor(Latitude)), width=0.1, position=dodge) +
  scale_y_continuous(name="Mean r")+ 
  scale_x_discrete(name="Time Period", labels=c("Drought", "Recovery")) + 
  scale_color_manual(values=colours_in_to_plot$Lat.Color, 
                     labels=round(colours_in_to_plot$Latitude,1),aesthetics=c("color", "fill")) +
  geom_hline(yintercept=0) +
  theme_classic() + theme(
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(size=11,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust=0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20, vjust=2, face="bold", hjust=0.5),
    strip.background = element_blank(), 
    strip.text.x = element_blank(),
    legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(linetype = 0)))
ggsave("Graphs/Demography/01e_recovery_cull.pdf",width=14, height = 8, units = "in")



#*******************************************************************************
### 3. Visualize decline over time for exemplar site (Old Fig 1D)
#*******************************************************************************

#Little Jamison 
dat_little <- dat %>% 
  filter(Site=="Little Jameson Creek") %>% 
  filter(Year<2015)
plot_3 <- ggplot(dat_little, aes(x=Year, y=r)) + 
  geom_point(size=4,shape=21,fill="#9BD7A4") +
  scale_y_continuous(name="r") + 
  scale_x_continuous(name="") +
  geom_hline(yintercept=0, linetype="dotted",size=0.8) + 
  theme_classic() + 
  theme(
    axis.text.x=element_text(size=18,face="bold"),
    axis.text.y=element_text(size=18,face="bold"),
    axis.title.x=element_text(color="black", size=20, vjust=0.5, face="bold"),
    axis.title.y=element_text(color="black", size=22, vjust=2, face="bold", hjust=0.5))
plot_3
ggsave("Graphs/Demography/sites/05_little_jamison.pdf",width=5, height=4, units="in")



#*******************************************************************************
### 4. Try expressing changes normalized to starting lambda per reviewer 2
#*******************************************************************************

r_means_norm <- r_means %>% 
  mutate(delta.r.drought.norm = (mean.r.drought-mean.r.pre)/abs(mean.r.pre),
#         delta.r.recovery.norm = (mean.r.recovery-mean.r.drought)/abs(mean.r.drought),
         delta.r.recovery.norm = (mean.r.recovery-mean.r.drought)/abs(mean.r.pre),
         rel.start = mean.r.pre/mean.r.pre)
write_csv(r_means_norm, "data/demography data/site_r_means.csv")

r_means_norm_long <- r_means_norm %>% 
  pivot_longer(cols=delta.r.drought.norm:rel.start, names_to="time", values_to="mean_r_norm")

colours_in<-r_means_norm_long %>% 
  dplyr::select(Latitude, Lat.Color) 

colours_in_to_plot<-unique(colours_in) #%>% filter(Lat.Color!="#FDB567") #remove this population because it doesn't have drought r value


r_means_norm_long_drought <- r_means_norm_long %>% filter(time!="delta.r.recovery.norm") %>% droplevels()
level_order_drought = c("rel.start", "delta.r.drought.norm")
dodge <- position_dodge(width=0.2)

colours_in<-r_means_norm_long_drought %>% 
  dplyr::select(Latitude, Lat.Color) 

colours_in_to_plot<-unique(colours_in) #%>% filter(Lat.Color!="#FDB567") #remove this population because it doesn't have drought r value

level_order_all = c("rel.start", "delta.r.drought.norm", "delta.r.recovery.norm")
dodge <- position_dodge(width=0.2)

# normalized trajectories (all relative to starting lambda pre-drought) = attempt at normalized U-shape
e <- ggplot(filter(r_means_norm_long, Site!="South Fork Middle Fork Tule"), aes(x=factor(time, level=level_order_all), y=mean_r_norm, group=Latitude, fill=as.factor(Latitude))) +
  facet_wrap(~Latitude, scale="free")+
  geom_point(shape=21, size=3, position=dodge) +
  #geom_errorbar(aes(ymin=ymin, ymax=ymax, colour=as.factor(Latitude)), width=0.1, position=dodge) +
  geom_line(aes(colour=as.factor(Latitude)), position=dodge, size=1.5) +
  scale_fill_manual(values=colours_in_to_plot$Lat.Color, 
                    labels = unique(r_means_norm_long$Paper_ID) ) +
  scale_color_manual(values=colours_in_to_plot$Lat.Color, 
                     labels = unique(r_means_norm_long$Paper_ID )) +
  scale_y_continuous(name="Normalized r")+ 
  scale_x_discrete(name="Time Period", labels=c("Pre", "Drght", "Recov")) + 
  geom_hline(yintercept=1) +
  theme_classic() + theme(
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(size=11, face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust=0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust=2, face="bold", hjust=0.5),
    strip.background=element_blank(), 
    strip.text.x=element_blank(),
    legend.title=element_blank())+
  guides(color = guide_legend(reverse = TRUE, override.aes = list(linetype = 0)),
         fill  = guide_legend(reverse = TRUE))
e
ggsave("Graphs/Demography/r_U_norm_faceted.pdf",width=5, height=4, units="in")

r_means_norm_long <- r_means_norm_long %>% 
  mutate(time_num = ifelse(time=="delta.r.drought.norm", 2, ifelse(time=="delta.r.recovery.norm", 3, 1)))

U_mod_norm_lin <- lm(mean_r_norm ~ time_num, data=r_means_norm_long)
summary(U_mod_lin)
U_mod_norm_quad <- lm(mean_r_norm ~ poly(time_num,2), data=r_means_norm_long)
summary(U_mod_quad)
AIC(U_mod_norm_lin, U_mod_norm_quad)

#normalized decline during drought (new Fig 1D)
c <- ggplot(r_means_norm_long_drought, aes(x=factor(time, level=level_order_drought), y=mean_r_norm, group=Latitude, fill=as.factor(Latitude))) +
  geom_point(shape=21, size=3, position=dodge) +
  #geom_errorbar(aes(ymin=ymin, ymax=ymax, colour=as.factor(Latitude)), width=0.1, position=dodge) +
  geom_line(aes(colour=as.factor(Latitude)), position=dodge, size=1.5) +
  scale_fill_manual(values=colours_in_to_plot$Lat.Color, 
                    labels = unique(r_means_norm_long_drought$Paper_ID) ) +
  scale_color_manual(values=colours_in_to_plot$Lat.Color, 
                     labels = unique(r_means_norm_long_drought$Paper_ID )) +
  scale_y_continuous(name="Normalized r")+ 
  scale_x_discrete(name="Time Period", labels=c("Pre-drought", "Drought")) + 
  geom_hline(yintercept=1) +
  theme_classic() + theme(
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(size=11, face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust=0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust=2, face="bold", hjust=0.5),
    strip.background=element_blank(), 
    strip.text.x=element_blank(),
    legend.title=element_blank())+
  guides(color = guide_legend(reverse = TRUE, override.aes = list(linetype = 0)),
         fill  = guide_legend(reverse = TRUE))
c
ggsave("Graphs/Demography/r_drought_norm.pdf",width=5, height=4, units="in")


#############################################################################################################

r_means_norm_long_recovery <- r_means_norm_long %>% filter(time!="delta.r.drought.norm") %>% droplevels()
level_order_drought = c("rel.start", "delta.r.recovery.norm")
dodge <- position_dodge(width=0.2)

colours_in<-r_means_norm_long_recovery %>% 
  dplyr::select(Latitude, Lat.Color) 

colours_in_to_plot<-unique(colours_in) #%>% filter(Lat.Color!="#FDB567") #remove this population because it doesn't have drought r value

d <- ggplot(r_means_norm_long_recovery, aes(x=factor(time, level=level_order_drought), y=mean_r_norm, group=Latitude, fill=as.factor(Latitude))) +
  geom_point(shape=21, size=3, position=dodge) +
  #geom_errorbar(aes(ymin=ymin, ymax=ymax, colour=as.factor(Latitude)), width=0.1, position=dodge) +
  geom_line(aes(colour=as.factor(Latitude)), position=dodge, size=1.5) +
  scale_fill_manual(values=colours_in_to_plot$Lat.Color, 
                    labels = unique(r_means_norm_long_recovery$Paper_ID) ) +
  scale_color_manual(values=colours_in_to_plot$Lat.Color, 
                     labels = unique(r_means_norm_long_recovery$Paper_ID )) +
  #scale_fill_manual(name = "Latitude (°N)",labels=round(r_means_norm_long_recovery$Latitude,1),
  #                  values=as.character(r_means_norm_long_recovery$Lat.Color)) +
  scale_y_continuous(name="Normalized r")+ 
  scale_x_discrete(name="Time Period", labels=c("Drought", "Recovery")) + 
  geom_hline(yintercept=1) +
  theme_classic() + theme(
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(size=11, face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust=0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust=2, face="bold", hjust=0.5),
    strip.background=element_blank(), 
    strip.text.x=element_blank(),
    legend.title=element_blank())+
  guides(color = guide_legend(reverse = TRUE, override.aes = list(linetype = 0)),
         fill  = guide_legend(reverse = TRUE))
d

# New Fig S1, option 2
plot_grid(a, b, c, d, labels = c("A", "B", "C", "D"), label_size = 16,nrow=2)
ggsave("Graphs/Demography/r_means_4panel.pdf",width=10, height=10, units="in")
