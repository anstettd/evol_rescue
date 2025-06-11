#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Produce demography graphs for publication
#### AUTHOR: Amy Angert & Daniel Anstett
#### DATE LAST MODIFIED: 20250610


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
### 1. Read in lambda estimates for each site and year 
#*******************************************************************************
dat <- read.csv("data/demography data/siteYear.lambda_2010-2019.csv") %>% 
  mutate(Site.Name=paste(Latitude, Site, sep="_")) %>% 
  mutate(r = log(lambda+0.01)) %>%  #add small value so lambdas=0 doesn't turn to NA)
  filter(Year!=2018) %>% 
  filter(Site!="Mill Creek")

# add color ramp column
color <- read_csv("data/genomic_data/pop_meta_data.csv") %>% 
  dplyr::select(Site, Lat.Color) %>% 
  unique()

dat <- left_join(dat, color)

dat_pre <- dat %>% filter(Year<2012)
dat_drought <- dat %>% filter(between(Year,2012,2014))
dat_recovery <- dat %>% filter(between(Year,2015,2017))
dat_post <- dat %>% filter(Year>2017)
dat_decline <- dat %>% filter(Year<2015)


#*******************************************************************************
### 1. Visualize estimates over time for all sites (Fig S1)
#*******************************************************************************

#Decline & Recovery
ggplot(dat, aes(x=Year, y=r)) +
  geom_point(data=filter(dat, Year<2012), col="grey") +
  geom_point(data=filter(dat, between(Year,2012,2014)), col="red") +
  geom_point(data=filter(dat, between(Year,2015,2017)), col="blue") +
  geom_point(data=filter(dat, Year>2017), col="grey") +
  scale_y_continuous(name="r")+ 
  scale_x_continuous(name="Year", limits=c(2010,2017))+
  geom_hline(yintercept=0, linetype="dotted") +
  facet_wrap(~Site.Name, scale="free") + theme_classic() + theme(
    axis.text.x=element_text(face="bold"),
    axis.text.y=element_text(size=11, face="bold"),
    axis.title.x=element_text(color="black", size=20, vjust=0.5, face="bold"),
    axis.title.y=element_text(color="black", size=20,vjust=2, face="bold", hjust=0.5))

ggsave("Graphs/Demography/01_decline_recovery.pdf",width=14, height = 8, units = "in")

#*******************************************************************************
### 2. Visualize mean estimates over each period for all sites
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

# all three time periods (alternative Fig S1)
dodge <- position_dodge(width=0.3)
ggplot(r_long, aes(x=factor(time, level=level_order_all), y=mean_r, group=Latitude, fill=as.factor(Latitude))) +
  geom_point(shape=21, size=3, position=dodge) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax, colour=as.factor(Latitude)), width=0.1, position=dodge) +
  geom_line(aes(colour=as.factor(Latitude)), position=dodge) +
  scale_color_manual(values=unique(r_long$Lat.Color), aesthetics=c("color", "fill")) +
  scale_y_continuous(name="Mean population growth rate")+ 
  scale_x_discrete(name="Time Period", labels=c("Pre-drought", "Drought", "Recovery")) + 
  geom_hline(yintercept=0) +
  theme_classic() + theme(
    axis.text.x=element_text(face="bold"),
    axis.text.y=element_text(size=11,face="bold"),
    axis.title.x=element_text(color="black", size=20, vjust=0.5, face="bold"),
    axis.title.y=element_text(color="black", size=20,vjust=2, face="bold", hjust=0.5),
    strip.background=element_blank(), 
    strip.text.x=element_blank(),
    legend.title=element_blank())

#Pre-drought to drought period
r_means_long_early <- r_means_long %>% filter(time!="mean.r.recovery") %>% droplevels()
level_order_pre = c("mean.r.pre", "mean.r.drought")
ggplot(r_means_long_early, aes(x=time, y=mean_r)) +
  geom_point(aes(x=factor(time, level=level_order_pre), fill=as.factor(Latitude)), shape=21, size=3) +
  geom_line(aes(group=Latitude), color=r_means_long_early$Lat.Color) +
  scale_color_manual(values=unique(r_means_long_early$Lat.Color), aesthetics=c("color", "fill")) +
  scale_y_continuous(name="Mean population growth rate")+ 
  scale_x_discrete(name="Time Period", labels=c("Pre-drought", "Drought")) + 
  geom_hline(yintercept=0) +
  theme_classic() + theme(
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(size=11, face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust=0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust=2, face="bold", hjust=0.5),
    strip.background=element_blank(), 
    strip.text.x=element_blank(),
    legend.title=element_blank())

#Drought to recovery period
r_means_long_late <- r_means_long %>% filter(time!="mean.r.pre") %>% droplevels()
level_order_post = c("mean.r.drought", "mean.r.recovery")
ggplot(r_means_long_late, aes(x=time, y=mean_r)) +
  geom_point(aes(x=factor(time, level=level_order_post), fill=as.factor(Latitude)), shape=21, size=3) +
  geom_line(aes(group=Latitude), color=r_means_long_late$Lat.Color) +
  scale_y_continuous(name="Mean population growth rate")+ 
  scale_x_discrete(name="Time Period", labels=c("Drought", "Recovery")) + 
  scale_color_manual(values=unique(r_means_long_late$Lat.Color), aesthetics=c("color", "fill")) +
  geom_hline(yintercept=0) +
  theme_classic() + theme(
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(size=11,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust=0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20, vjust=2, face="bold", hjust=0.5),
    strip.background = element_blank(), 
    strip.text.x = element_blank(),
    legend.title = element_blank())


# sensitivity test: remove 3 extirpated populations
r_means_cull <- r_means %>% filter(mean.r.recovery>-2.5)
r_means_cull_long <- r_means_cull %>% pivot_longer(cols=mean.r.pre:mean.r.recovery, names_to="time", values_to="mean_r")
r_means_cull_long_early <- r_means_cull_long %>% filter(time!="mean.r.recovery") %>% droplevels()
r_means_cull_long_late <- r_means_cull_long %>% filter(time!="mean.r.pre") %>% droplevels()

ggplot(r_means_cull_long_early, aes(x=time, y=mean_r)) +
  geom_point(aes(x=factor(time, level=level_order_pre), fill=as.factor(Latitude)), shape=21, size=2) +
  geom_line(aes(group=Latitude), color=r_means_cull_long_early$Lat.Color) +
  scale_y_continuous(name="Mean population growth rate")+ 
  scale_color_manual(values=unique(r_means_cull_long_early$Lat.Color), aesthetics=c("color", "fill")) +
  scale_x_discrete(name="Time Period", labels=c("Pre-Drought", "Drought")) + 
  scale_color_manual(values=unique(r_means_long_late$Lat.Color), aesthetics=c("color", "fill")) +
  geom_hline(yintercept=0) +
  theme_classic() + theme(
    axis.text.x=element_text(face="bold"),
    axis.text.y=element_text(size=11, face="bold"),
    axis.title.x=element_text(color="black", size=20, vjust=0.5, face="bold"),
    axis.title.y=element_text(color="black", size=20, vjust=2, face="bold", hjust=0.5),
    strip.background=element_blank(), 
    strip.text.x=element_blank(),
    legend.title = element_blank())

ggplot(r_means_cull_long_late, aes(x=time, y=mean_r)) +
  geom_point(aes(x=factor(time, level=level_order_post), fill=as.factor(Latitude)), shape=21, size=2) +
  geom_line(aes(group=Latitude), color=r_means_cull_long_late$Lat.Color) +
  scale_color_manual(values=unique(r_means_cull_long_late$Lat.Color), aesthetics=c("color", "fill")) +
  scale_y_continuous(name="Mean population growth rate")+ 
  scale_x_discrete(name="Time Period", labels=c("Drought", "Recovery")) + 
  scale_color_manual(values=unique(r_means_long_late$Lat.Color), aesthetics=c("color", "fill")) +
  geom_hline(yintercept=0) +
  theme_classic() + theme(
    axis.text.x=element_text(face="bold"),
    axis.text.y=element_text(size=11,face="bold"),
    axis.title.x=element_text(color="black", size=20, vjust=0.5, face="bold"),
    axis.title.y=element_text(color="black", size=20,vjust=2, face="bold", hjust=0.5),
    strip.background=element_blank(), 
    strip.text.x=element_blank(),
    legend.title=element_blank())


#*******************************************************************************
### 3. Visualize decline over time for exemplar site (Fig 1D)
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



