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

#Decline & Recovery (Old Figure S1)
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

# Alternative depiction in same style as climate-year graphs
r_spag_plot<-ggplot() + 
  geom_line(stat="smooth",data=dat, aes(y=r, x=Year, group=Latitude, col=as.factor(Latitude)), alpha=0.75, linewidth = 1.5, se=FALSE) + 
  xlab("Year") + 
  ylab("Mean Population Growth rate") +
  geom_hline(yintercept=0) +
  theme_classic() +
  scale_x_continuous(breaks=c(2010,2012,2014,2016,2018))+
  scale_color_manual(values=rev(unique(dat$Lat.Color))) 
r_spag_plot
ggsave("Graphs/Demography/01_decline_recovery_spaghetti.pdf",width=6, height = 8, units = "in")

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

# all three time periods (alternative Fig S1, option 1)
dodge <- position_dodge(width=0.2)
ggplot(r_long, aes(x=factor(time, level=level_order_all), y=mean_r, group=Latitude, fill=as.factor(Latitude))) +
  geom_line(aes(colour=as.factor(Latitude)), position=dodge, size=1.5) +
  geom_point(shape=21, size=3, position=dodge) +
  #geom_errorbar(aes(ymin=ymin, ymax=ymax, colour=as.factor(Latitude)), width=0.1, position=dodge) +
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
ggsave("Graphs/Demography/01b_decline_recovery_rmeans.pdf",width=14, height = 8, units = "in")

#Pre-drought to drought period
dodge <- position_dodge(width=0.2)
r_long_early <- r_long %>% filter(time!="recovery") %>% droplevels()
level_order_pre = c("pre", "drought")
a <- ggplot(r_long_early, aes(x=factor(time, level=level_order_pre), y=mean_r, group=Latitude, fill=as.factor(Latitude))) +
  geom_point(shape=21, size=3, position=dodge) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax, colour=as.factor(Latitude)), width=0.1, position=dodge) +
  geom_line(aes(colour=as.factor(Latitude)), position=dodge, size=1.5) +
  scale_color_manual(values=unique(r_long_early$Lat.Color), aesthetics=c("color", "fill")) +
  scale_y_continuous(name="Mean r")+ 
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
ggsave(a, "Graphs/Demography/01c_decline.pdf",width=14, height = 8, units = "in")

#Drought to recovery period
r_long_late <- r_long %>% filter(time!="pre") %>% droplevels()
level_order_post = c("drought", "recovery")
b <- ggplot(r_long_late, aes(x=factor(time, level=level_order_post), y=mean_r, group=Latitude, fill=as.factor(Latitude))) +
  geom_point(shape=21, size=3, position=dodge) +
  geom_line(aes(colour=as.factor(Latitude)), position=dodge, size=1.5) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax, colour=as.factor(Latitude)), width=0.1, position=dodge) +
  scale_y_continuous(name="Mean r")+ 
  scale_x_discrete(name="Time Period", labels=c("Drought", "Recovery")) + 
  scale_color_manual(values=unique(r_long_late$Lat.Color), aesthetics=c("color", "fill")) +
  geom_hline(yintercept=0) +
  theme_classic() + theme(
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(size=11,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust=0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20, vjust=2, face="bold", hjust=0.5),
    strip.background = element_blank(), 
    strip.text.x = element_blank(),
    legend.title = element_blank())
ggsave("Graphs/Demography/01d_recovery.pdf",width=14, height = 8, units = "in")


# sensitivity test: remove 3 extirpated populations
r_long_late_cull <- filter(r_long_late, mean_r>-2)
ggplot(r_long_late_cull, aes(x=factor(time, level=level_order_post), y=mean_r, group=Latitude, fill=as.factor(Latitude))) +
  geom_point(shape=21, size=3, position=dodge) +
  geom_line(aes(colour=as.factor(Latitude)), position=dodge, size=1.5) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax, colour=as.factor(Latitude)), width=0.1, position=dodge) +
  scale_y_continuous(name="Mean r")+ 
  scale_x_discrete(name="Time Period", labels=c("Drought", "Recovery")) + 
  scale_color_manual(values=unique(r_long_late_cull$Lat.Color), aesthetics=c("color", "fill")) +
  geom_hline(yintercept=0) +
  theme_classic() + theme(
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(size=11,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust=0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20, vjust=2, face="bold", hjust=0.5),
    strip.background = element_blank(), 
    strip.text.x = element_blank(),
    legend.title = element_blank())
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
         delta.r.recovery.norm = (mean.r.recovery-mean.r.drought)/abs(mean.r.drought),
         rel.start = mean.r.pre/mean.r.pre)
write_csv(r_means_norm, "data/demography data/site_r_means.csv")

r_means_norm_long <- r_means_norm %>% 
  pivot_longer(cols=delta.r.drought.norm:rel.start, names_to="time", values_to="mean_r_norm")

r_means_norm_long_drought <- r_means_norm_long %>% filter(time!="delta.r.recovery.norm") %>% droplevels()
level_order_drought = c("rel.start", "delta.r.drought.norm")
dodge <- position_dodge(width=0.2)

#normalized decline during drought (new Fig 1D)
c <- ggplot(r_means_norm_long_drought, aes(x=factor(time, level=level_order_drought), y=mean_r_norm, group=Latitude, fill=as.factor(Latitude))) +
  geom_point(shape=21, size=3, position=dodge) +
  #geom_errorbar(aes(ymin=ymin, ymax=ymax, colour=as.factor(Latitude)), width=0.1, position=dodge) +
  geom_line(aes(colour=as.factor(Latitude)), position=dodge, size=1.5) +
  scale_color_manual(values=unique(r_means_norm_long_drought$Lat.Color), aesthetics=c("color", "fill")) +
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
    legend.title=element_blank())
ggsave("Graphs/Demography/r_drought_norm.pdf",width=5, height=4, units="in")

r_means_norm_long_recovery <- r_means_norm_long %>% filter(time!="delta.r.drought.norm") %>% droplevels()
level_order_drought = c("rel.start", "delta.r.recovery.norm")
dodge <- position_dodge(width=0.2)
d <- ggplot(r_means_norm_long_recovery, aes(x=factor(time, level=level_order_drought), y=mean_r_norm, group=Latitude, fill=as.factor(Latitude))) +
  geom_point(shape=21, size=3, position=dodge) +
  #geom_errorbar(aes(ymin=ymin, ymax=ymax, colour=as.factor(Latitude)), width=0.1, position=dodge) +
  geom_line(aes(colour=as.factor(Latitude)), position=dodge, size=1.5) +
  scale_color_manual(values=unique(r_means_norm_long_recovery$Lat.Color), aesthetics=c("color", "fill")) +
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
    legend.title=element_blank())

# New Fig S1, option 2
plot_grid(a, b, c, d, nrow=2)
ggsave("Graphs/Demography/r_means_4panel.pdf",width=10, height=8, units="in")
