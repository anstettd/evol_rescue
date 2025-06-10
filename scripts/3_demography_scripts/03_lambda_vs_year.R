#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Calculate metrics of the rate of demographic decline during drought and rate of recovery after drought
#### AUTHOR: Amy Angert
#### DATE LAST MODIFIED: 20250608


#*******************************************************************************
#### 0. Clean workspace and load required packages
#*******************************************************************************

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

# Make vector of packages needed
packages_needed <- c("tidyverse", "RColorBrewer", "FSA")

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
dat <- read.csv("data/demography data/siteYear.lambda_2010-2019.csv") %>% dplyr::select(-Region, -RegionRank) %>% filter(Site!="Mill Creek")
 
# check distributions
hist(dat$lambda) #extreme right skew
hist(log(dat$lambda)) #much better

# express as population growth as intrinsic rate of increase (log lambda)
dat <- dat %>% mutate(r = log(lambda+0.01)) #adding small value lambdas of 0 don't turn to NA

#*******************************************************************************
### 2. Visualize estimates over time for each site
#*******************************************************************************

ggplot(dat, aes(x=Year, y=r)) + 
  geom_point() +
  geom_smooth(data=filter(dat, Year<2015), method="lm", se=FALSE, col="red") +
  geom_smooth(data=filter(dat, Year>2014), method="lm", se=FALSE, col="blue") +
  ylab("Intrinsic rate of increase") +
  #ylim(0,2) +
  geom_hline(yintercept=0, linetype="dotted") +
  facet_wrap(~Site, scale="free", nrow=3) +
  theme_classic() #+
  #theme(strip.background = element_blank(), strip.text.x = element_blank(),
   #     legend.title = element_blank())

# Note: some sites' slopes (e.g., Buck Meadows, Mill) are affected by 2015-16, which had very high recruitment and we assume indicated relatively early recovery. This is part of the rationale for calculating decline only until 2014-15.

# Note: more generally, slopes are often pulled up or down by single years with extremely high lambdas. High lambda values are dramatically higher with addition of 2016-2019 data and its effects on vital rate model selection and model fits

# Note: Canton, NFMF Tule, SFMF Tule, & Redwood have only one annual transition estimate during drought recovery (because of fire and flood closures that prevented site access in 2016 or 2017)


#*******************************************************************************
### 3A. Calculate geometric mean population growth rate BEFORE DROUGHT for each site
# Before drought = 2010-11, 2011-12 transitions
#*******************************************************************************

dat.mean.pre <- dat %>%
  group_by(Latitude, Site) %>%
  filter(Year==2010|Year==2011) %>%
  na.omit() %>%
  dplyr::summarize(mean.r.pre = mean(r), #arithmetic mean of logs
                   geomean.l.pre = exp(mean(log(lambda))),#geometric mean of unlogged values
                   sd.r.pre = sd(r),
                   se.r.pre = se(r)) 
hist(dat.mean.pre$mean.r.pre) # better 
hist(dat.mean.pre$geomean.l.pre) #long right tail

#*******************************************************************************
### 3B. Calculate mean lambda DURING CORE DROUGHT for each site 
# Core drought = 2012-13, 2013-14, 2014-15 transitions
#*******************************************************************************

dat.mean.drought <- dat %>% 
  group_by(Latitude, Site) %>% 
  filter(Year==2012|Year==2013|Year==2014) %>% 
  na.omit() %>% 
  dplyr::summarize(mean.r.drought = mean(r), 
                   geomean.l.drought = exp(mean(log(lambda))),
                   sd.r.drought = sd(r),
                   se.r.drought = se(r)) 
hist(dat.mean.drought$mean.r.drought) #left-skewed
hist(dat.mean.drought$geomean.l.drought) #better

#*******************************************************************************
### 3C. Calculate mean lambda DURING POST-DROUGHT for each site
# Post drought = 2015-16, 2016-17, 2017-18 transitions
#*******************************************************************************

dat.mean.recovery <- dat %>% 
  group_by(Latitude, Site) %>% 
  filter(Year==2015|Year==2016|Year==2017) %>% 
  na.omit() %>% 
  dplyr::summarize(mean.r.recovery = mean(r), 
            geomean.l.recovery = exp(mean(log(lambda))),
            sd.r.recovery = sd(r),
            se.r.recovery = se(r))
hist(dat.mean.recovery$mean.r.recovery) #better
hist(dat.mean.recovery$geomean.l.recovery) #right-skewed

# Add Latitude and other covariates back in
covar <- read_csv("data/genomic_data/pop_meta_data.csv") %>% 
  dplyr::select(Site, Paper_ID, Latitude, Longitude, Lat.Color) %>% 
  unique()

mean.lambda <- left_join(dat.mean.drought,dat.mean.pre) %>% left_join(dat.mean.recovery) %>% left_join(covar) # Join to slopes

####### TABLE S1 #########
# Save to .csv file 
write.csv(mean.lambda,"data/demography data/siteYear.lambda_responses_2010-2019.csv",row.names=FALSE)


#*******************************************************************************
### 4. Compare population growth rates during different periods

#*******************************************************************************

r_means <- mean.lambda %>% 
  dplyr::select(Site,Paper_ID,Latitude,Longitude,Lat.Color,mean.r.pre,mean.r.drought,mean.r.recovery) %>%  
  arrange(Latitude)

# pre to drought
t.test(x=r_means$mean.r.pre, y=r_means$mean.r.drought, paired=TRUE)

# drought to recovery
t.test(x=r_means$mean.r.recovery, y=r_means$mean.r.drought, paired=TRUE)

# sensitivity test: remove 3 extirpated populations
r_means_cull <- r_means %>% filter(mean.r.recovery>-2.5)
t.test(x=r_means_cull$mean.r.pre, y=r_means_cull$mean.r.drought, paired=TRUE)
t.test(x=r_means_cull$mean.r.recovery, y=r_means_cull$mean.r.drought, paired=TRUE)
