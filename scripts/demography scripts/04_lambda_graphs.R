#### PROJECT: Genomic offsets and demographic trajectories of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Produce demography graphs for publication
#### AUTHOR: Amy Ang & Daniel Anstett
#### DATE LAST MODIFIED: 20230815


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
### 1. Read in lambda estimates for each site and year & manupulate
#*******************************************************************************
dat <- read.csv("data/demography data/siteYear.lambda_2010-2019.csv")
pop_key <- read_csv("data/genomic_data/offset_pop_beagle.csv") %>% dplyr::select(Site, Paper_ID) #just to get translation of pop names <--> numbers
pop_key[20,1] <- "Mill Creek"
pop_key[20,2] <- 12
dat<-left_join(dat,pop_key,by="Site") %>% mutate(lat.2 = round(Latitude, 2), 
                                                 Site.Lat=paste(lat.2,Paper_ID, sep = "_"),
                                                 Site.Name=paste(lat.2,Site, sep = "_")
                                                 ) %>% select(-lat.2) %>% filter(Site!="Mill Creek")
dat_decline <- dat %>% filter(Year<2015)
dat_recovery <- dat %>% filter(Year>2014)


#*******************************************************************************
### 2. Visualize estimates over time for all sites
#*******************************************************************************

#Decline & Recovery
ggplot(dat, aes(x=Year, y=lambda)) + geom_point() +
  geom_smooth(data=filter(dat, Year<2015), method="lm", se=FALSE, col="red") +
  geom_smooth(data=filter(dat, Year>2014), method="lm", se=FALSE, col="blue") +
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

#Decline
ggplot(dat_decline, aes(x=Year, y=lambda)) + geom_point() +
  geom_smooth(method="lm", se=FALSE, col="red") +
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

ggsave("Graphs/Demography/02_decline.pdf",width=14, height = 8, units = "in")

#Recovery
ggplot(dat_recovery, aes(x=Year, y=lambda)) + geom_point() +
  geom_smooth(method="lm", se=FALSE, col="blue") +
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

ggsave("Graphs/Demography/03_recovery.pdf",width=14, height = 8, units = "in")


#*******************************************************************************
### 2. Visualize estimates over time for select sites
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





