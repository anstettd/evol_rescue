#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Produce demography graphs for publication
#### AUTHOR: Amy Angert & Daniel Anstett
#### DATE LAST MODIFIED: 20240627


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
  mutate(Site.Name=paste(Latitude,Site, sep = "_")) 

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
### 1. Visualize declines over time for all sites (Fig S1)
#*******************************************************************************

#Decline only
ggplot(dat_decline, aes(x=Year, y=lambda)) +
  geom_point(aes(fill=as.factor(round(Latitude, 1))), shape=21, size=4) +
  geom_smooth(method=lm, se=FALSE, aes(color=as.factor(round(Latitude, 1)))) + 
  scale_color_manual(values=rev(unique(dat_decline$Lat.Color)), aesthetics = c("color", "fill")) +
  scale_y_continuous(name="Lambda")+ scale_x_continuous(name="Year")+
  geom_hline(yintercept=1, linetype="dotted") +
  facet_wrap(~Site.Name, scale="free") + theme_classic() + theme(
    #facet_wrap(~Site.Lat, scale="free") + theme_classic() + theme(
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(size=11,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5)) + theme(strip.background = element_blank(), strip.text.x = element_blank(), legend.title = element_blank())

ggsave("Graphs/Demography/01_decline.pdf",width=14, height = 8, units = "in")


#*******************************************************************************
### 3. Visualize estimates over time for all sites (Fig S9)
#*******************************************************************************

#Decline & Recovery
ggplot(dat, aes(x=Year, y=lambda)) +
  geom_point(data=filter(dat, Year<2012), col="grey") +
  geom_point(data=filter(dat, between(Year,2012,2014)), col="red") +
  geom_point(data=filter(dat, between(Year,2015,2017)), col="blue") +
  geom_point(data=filter(dat, Year>2017), col="grey") +
  scale_y_continuous(name="Lambda")+ scale_x_continuous(name="Year")+
  geom_hline(yintercept=1, linetype="dotted") +
  facet_wrap(~Site.Name, scale="free") + theme_classic() + theme(
  #facet_wrap(~Site.Lat, scale="free") + theme_classic() + theme(
    axis.text.x = element_text(face="bold"),
    axis.text.y = element_text(size=11,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
#+
#theme(strip.background = element_blank(), strip.text.x = element_blank(),
#     legend.title = element_blank())

ggsave("Graphs/Demography/01_decline_recovery.pdf",width=14, height = 8, units = "in")


#*******************************************************************************
### 4. Visualize estimates over time for select sites
#*******************************************************************************

#Kitchen creek
dat_kitchen <- dat %>% filter(Paper_ID==15)
plot_1 <- ggplot(dat_kitchen, aes(x=Year, y=lambda)) + geom_point(size=2) +
  geom_smooth(data=filter(dat_kitchen, Year<2015), method="lm", se=FALSE, col="red") +
  geom_smooth(data=filter(dat_kitchen, Year>2014), method="lm", se=FALSE, col="blue") +
  scale_y_continuous(name="Lambda")+ scale_x_continuous(name="")+
  geom_hline(yintercept=1, linetype="dotted") + theme_classic() + theme(
    #facet_wrap(~Site.Lat, scale="free") + theme_classic() + theme(
    axis.text.x = element_text(size=18,face="bold"),
    axis.text.y = element_text(size=18,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=22,vjust = 2, face="bold",hjust=0.5))
ggsave("Graphs/Demography/sites/01_kitchen.pdf",width=6, height = 4, units = "in")

#Wawona
dat_Wawona <- dat %>% filter(Paper_ID==5)
plot_2 <- ggplot(dat_Wawona, aes(x=Year, y=lambda)) + geom_point(size=2) +
  geom_smooth(data=filter(dat_Wawona, Year<2015), method="lm", se=FALSE, col="red") +
  geom_smooth(data=filter(dat_Wawona, Year>2014), method="lm", se=FALSE, col="blue") +
  scale_y_continuous(name="Lambda")+ scale_x_continuous(name="")+
  geom_hline(yintercept=1, linetype="dotted") + theme_classic() + theme(
    #facet_wrap(~Site.Lat, scale="free") + theme_classic() + theme(
    axis.text.x = element_text(size=18,face="bold"),
    axis.text.y = element_text(size=18,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=22,vjust = 2, face="bold",hjust=0.5))
ggsave("Graphs/Demography/sites/02_wawona.pdf",width=6, height = 4, units = "in")

#Little Jamison
dat_little <- dat %>% filter(Paper_ID==7)
plot_3 <- ggplot(dat_little, aes(x=Year, y=lambda)) + geom_point(size=2) +
  geom_smooth(data=filter(dat_little, Year<2015), method="lm", se=FALSE, col="red") +
  geom_smooth(data=filter(dat_little, Year>2014), method="lm", se=FALSE, col="blue") +
  scale_y_continuous(name="Lambda")+ scale_x_continuous(name="")+
  geom_hline(yintercept=1, linetype="dotted") + theme_classic() + theme(
    #facet_wrap(~Site.Lat, scale="free") + theme_classic() + theme(
    axis.text.x = element_text(size=18,face="bold"),
    axis.text.y = element_text(size=18,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=22,vjust = 2, face="bold",hjust=0.5))
ggsave("Graphs/Demography/sites/03_little_jamison.pdf",width=6, height = 4, units = "in")

#Coast Fork of Willamette
dat_coast <- dat %>% filter(Paper_ID==55)
plot_4 <- ggplot(dat_coast, aes(x=Year, y=lambda)) + geom_point(size=2) +
  geom_smooth(data=filter(dat_coast, Year<2015), method="lm", se=FALSE, col="red") +
  geom_smooth(data=filter(dat_coast, Year>2014), method="lm", se=FALSE, col="blue") +
  scale_y_continuous(name="Lambda")+ scale_x_continuous(name="")+
  geom_hline(yintercept=1, linetype="dotted") + theme_classic() + theme(
    #facet_wrap(~Site.Lat, scale="free") + theme_classic() + theme(
    axis.text.x = element_text(size=18,face="bold"),
    axis.text.y = element_text(size=18,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=22,vjust = 2, face="bold",hjust=0.5))
ggsave("Graphs/Demography/sites/04_coast_fork.pdf",width=6, height = 4, units = "in")


#*******************************************************************************
### 5. Visualize decline over time for select sites (Fig 1D)
#*******************************************************************************

#Little Jamison
dat_little <- dat %>% filter(Paper_ID==7) %>% filter(Year<2015)
plot_3 <- ggplot(dat_little, aes(x=Year, y=lambda)) + geom_point(size=4,shape=21,fill="#9BD7A4") +
  geom_smooth(data=filter(dat_little, Year<2015), method="lm", se=FALSE, col="#9BD7A4",size=1.5) +
  #geom_smooth(data=filter(dat_little, Year>2014), method="lm", se=FALSE, col="blue") +
  scale_y_continuous(name="Lambda")+ scale_x_continuous(name="")+
  geom_hline(yintercept=1, linetype="dotted",size=0.8) + theme_classic() + theme(
    #facet_wrap(~Site.Lat, scale="free") + theme_classic() + theme(
    axis.text.x = element_text(size=18,face="bold"),
    axis.text.y = element_text(size=18,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=22,vjust = 2, face="bold",hjust=0.5))
plot_3
ggsave("Graphs/Demography/sites/05_little_jamison.pdf",width=5, height = 4, units = "in")
