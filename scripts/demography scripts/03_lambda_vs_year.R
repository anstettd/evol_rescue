#### PROJECT: Genomic offsets and demographic trajectories of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Calculate slopes of lambda versus year as a metric of the rate of demographic decline during drought
#### AUTHOR: Amy Angert
#### DATE LAST MODIFIED: 20230621


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

# Note: Canton, NFMF Tule, SFMF Tule, & Redwood have only two annual transition estimates during drought recovery (because of fire and flood closures in 2016 and 2017), so they should be removed from calculations of demographic recovery

#*******************************************************************************
### 3A. Calculate slopes of lambda over time DURING DROUGHT-INDUCED DECLINE for each site
#*******************************************************************************

# with cleaned lambdas, 2010-11 through 2014-15
site.vec <- unique(dat$Site)
slopes.lam.decline <- c()
site.lam <- c()

for (i in 1:length(site.vec)) {
  dat.site <- dat %>% filter(Site==site.vec[i],
                             Year<2015) 
  mod <- lm(lambda ~ Year, dat.site)
  slopes.lam.decline[i] <- coefficients(mod)[2]
  site.lam[i] <- site.vec[i]
}

slopes.lambda.decline <- bind_cols(site.lam, slopes.lam.decline) %>% 
  dplyr::select(Site=...1, Lambda.Slope.1011.1415=...2) 


#*******************************************************************************
### 3B. Calculate slopes of lambda over time DURING POST-DROUGHT RECOVERY for each site
#*******************************************************************************

# with cleaned lambdas, 2015-16 through 2018-19
site.vec <- unique(dat$Site)
slopes.lam.recovery <- c()
site.lam <- c()

for (i in 1:length(site.vec)) {
  dat.site <- dat %>% filter(Site==site.vec[i],
                             Year>2014) 
  mod <- lm(lambda ~ Year, dat.site)
  slopes.lam.recovery[i] <- coefficients(mod)[2]
  site.lam[i] <- site.vec[i]
}

slopes.lambda.recovery <- bind_cols(site.lam, slopes.lam.recovery) %>% 
  dplyr::select(Site=...1, Lambda.Slope.1516.1819=...2) 
# some of these are negative because of an extremely high early value

# Add Latitude and other covariates back in
covar <- dat %>% 
  dplyr::select(Site, Latitude, Longitude, Elevation, Region, RegionRank) %>% 
  unique()

slopes.lambda <- left_join(slopes.lambda.decline, slopes.lambda.recovery) %>% left_join(covar, by="Site") 


#*******************************************************************************
### 4A. Calculate mean lambda DURING CORE DROUGHT for each site
#*******************************************************************************

dat.mean.drought <- dat %>% 
  group_by(Latitude, Site) %>% 
  filter(Year==2012|Year==2013|Year==2014) %>% 
  na.omit() %>% 
  summarize(mean.lambda.drought = exp(mean(log(lambda)))) #GEOMETRIC mean

#*******************************************************************************
### 4B. Calculate mean lambda DURING POST-DROUGHT for each site
#*******************************************************************************

dat.mean.recovery <- dat %>% 
  group_by(Latitude, Site) %>% 
  filter(Year==2015|Year==2016|Year==2017|Year==2018|Year==2019) %>% 
  # NOTE: should we restrict this to a core 3 years, as above?
  na.omit() %>% 
  summarize(mean.lambda.recovery = exp(mean(log(lambda)))) #GEOMETRIC mean

# Join to slopes
demog.response <- left_join(slopes.lambda, dat.mean.drought) %>%
  left_join(dat.mean.recovery) %>% 
  dplyr::select(Site, Latitude, Longitude, Elevation, Region, lambda.slope.decline=Lambda.Slope.1011.1415, lambda.slope.recovery=Lambda.Slope.1516.1819, lambda.mean.drought=mean.lambda.drought, lambda.mean.recovery=mean.lambda.recovery)

# Save to .csv file 
write.csv(demog.response,"data/demography data/siteYear.lambda_responses_2010-2019.csv",row.names=FALSE)

