#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Calculate slopes of lambda versus year as a metric of the rate of demographic decline during drought
#### AUTHOR: Amy Angert
#### DATE LAST MODIFIED: 20240627


#*******************************************************************************
#### 0. Clean workspace and load required packages
#*******************************************************************************

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

# Make vector of packages needed
packages_needed <- c("tidyverse", "RColorBrewer")

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
dat <- read.csv("data/demography data/siteYear.lambda_2010-2019.csv")


#*******************************************************************************
### 2. Visualize estimates over time for each site
#*******************************************************************************

# N-S color gradient
lat_cols=colorRampPalette(brewer.pal(11,"Spectral"))
n.sites <- length(unique(dat$Site))
color.list <- lat_cols(n.sites)

ggplot(dat, aes(x=Year, y=lambda)) + #, color=as.factor(round(Latitude, 1))
  geom_point() +
  geom_smooth(data=filter(dat, Year<2015), method="lm", se=FALSE, col="red") +
  geom_smooth(data=filter(dat, Year>2014), method="lm", se=FALSE, col="blue") +
  #scale_color_manual(values=color.list) +
  ylab("Lambda") +
  #ylim(0,2) +
  geom_hline(yintercept=1, linetype="dotted") +
  facet_wrap(~Site, scale="free", nrow=3) +
  theme_classic() #+
  #theme(strip.background = element_blank(), strip.text.x = element_blank(),
   #     legend.title = element_blank())

# Note: some sites' slopes (e.g., Buck Meadows, Mill) are affected by 2015-16, which had very high recruitment and we assume indicated relatively early recovery. This is part of the rationale for calculating rates of decline as slope until 2014-15.

# Note: more generally, slopes are often pulled up or down by single years with extremely high lambdas. High lambda values are dramatically higher with addition of 2016-2019 data and its effects on vital rate model selection and model fits

# Note: Mill Creek only has two annual transition estimates during drought (because of 100% plot wash-out in 2010 and flooding that prevented site access in 2013), so Mill Creek should be removed from calculations of demographic declines  

# Note: Canton, NFMF Tule, SFMF Tule, & Redwood have only one annual transition estimate during drought recovery (because of fire and flood closures in 2016 and 2017), so they should be removed from calculations of demographic recovery


#*******************************************************************************
### 3A. Calculate mean lambda DURING CORE DROUGHT for each site 
# Core drought = 2012-13, 2013-14, 2014-15 transitions
#*******************************************************************************

dat.mean.drought <- dat %>% 
  group_by(Latitude, Site) %>% 
  filter(Year==2012|Year==2013|Year==2014) %>% 
  na.omit() %>% 
  summarize(mean.lambda.drought = exp(mean(log(lambda)))) #GEOMETRIC mean

#*******************************************************************************
### 3B. Calculate mean lambda DURING POST-DROUGHT for each site
# Post drought = 2015-16, 2016-17, 2017-18 transitions
#*******************************************************************************

dat.mean.recovery <- dat %>% 
  group_by(Latitude, Site) %>% 
  filter(Year==2015|Year==2016|Year==2017) %>% 
  na.omit() %>% 
  summarize(mean.lambda.recovery = exp(mean(log(lambda)))) #GEOMETRIC mean

# Add Latitude and other covariates back in
covar <- dat %>% 
  dplyr::select(Site, Latitude, Longitude, Elevation, Region, RegionRank) %>% 
  unique()

mean.lambda <- left_join(dat.mean.drought, dat.mean.recovery) %>% left_join(covar) # Join to slopes

# Save to .csv file 
write.csv(mean.lambda,"data/demography data/siteYear.lambda_responses_2010-2019.csv",row.names=FALSE)

