#### PROJECT: Genomic offsets and demographic trajectories of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Create data frame of vital rate parameters and build integral projection models to obtain estimates of annual lambdas for each population
#### AUTHOR: Seema Sheth and Amy Angert
#### DATE LAST MODIFIED: 20230619


#*******************************************************************************
#### 0. Clean workspace and load required packages
#*******************************************************************************

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

# Make vector of packages needed
packages_needed <- c("lme4", "glmmTMB", "tidyverse")

# Install packages needed (if not already installed)
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
}

# if there are errors related to glmmTMB, try: install.packages("glmmTMB", type="source")

# Load packages needed
for (i in 1:length(packages_needed)){
  library(packages_needed[i], character.only = TRUE)
}


#*******************************************************************************
#### 2. Read in vital rate data frames ###
#*******************************************************************************

data <- read.csv("data/demography data/Mcard_demog_data_2010-2019_cleanindivs.csv")
data$Site = factor(data$Site)
data$Year = factor(data$Year)
data$Year = factor(data$Year)
data$SiteYear = factor(data$SiteYear)

site_fruit_count_data <- read.csv("data/demography data/Mcard_demog_data_2010-2019_seedinput.csv")
site_fruit_count_data$Site = factor(site_fruit_count_data$Site)
site_fruit_count_data$Year = factor(site_fruit_count_data$Year)
site_fruit_count_data$Year = factor(site_fruit_count_data$Year)
site_fruit_count_data$SiteYear = factor(site_fruit_count_data$SiteYear)

#*******************************************************************************
#### 1. Create global survival, growth and fecundity models using data from all sites ###
#*******************************************************************************

# Set up data frames of model parameters
params=c()
growth=c()
flower=c()
fruit=c()


  #*******************************************************************************
  ### 1A. Survival ###
  #*******************************************************************************

  # Read in top survival model output (Formula: Surv ~ logSize + (1|Year/Site))
  surv.reg <- load("data/demography data/surv.reg.rda")

  # Store model coefficients for each Site:Year
  params$SiteYear = rownames(coefficients(s4)$'Site:Year')
  params$surv.int=coefficients(s4)$'Site:Year'[,1] 
  params$surv.slope=coefficients(s4)$'Site:Year'[,2] 

  # make a data frame
  params <- data.frame(params)
  
  #*******************************************************************************
  ### 1B. Growth ###
  #*******************************************************************************
  
  # Read in top growth model output (Formula: logSizeNext ~ logSize + (logSize|Year/Site))
  growth.reg <- load("data/demography data/growth.reg.rda")
  
  # Store model coefficients
  growth$SiteYear = rownames(coefficients(g3)$'Site:Year')
  growth$growth.int=coefficients(g3)$'Site:Year'[,1] 
  growth$growth.slope=coefficients(g3)$'Site:Year'[,2] 
  growth$growth.sd=rep(sigma(g3), dim(coefficients(g3)$'Site:Year')[1])
 
  # Note: script "dnorm.R" determines that we do not need to switch to the cdf estimation for the growth function
 
  # Make a data frame
  growth=data.frame(growth)
  
  #*******************************************************************************
  ### 1C. Flowering ###
  #*******************************************************************************
  
  # Read in top flowering model output (Formula: Fec0 ~ logSize + (1|Year/Site))
  flowering.reg=load("data/demography data/flowering.reg.rda")

  # Store model coefficients
  flower$SiteYear = rownames(coefficients(fl4)$'Site:Year')
  flower$flowering.int=coefficients(fl4)$'Site:Year'[,1] 
  flower$flowering.slope=coefficients(fl4)$'Site:Year'[,2] 
  
  # make a data frame
  flower=data.frame(flower)

    #*******************************************************************************
  ### 1D. Fruit number (untransformed) using negative binomial regression ###
  #*******************************************************************************
  # Read in top model output for fruit.reg (Formula: Fec1 ~ logSize + (logSize|Year) + (logSize|Site)   
  fruit.reg=load("data/demography data/fruit.reg.rda")

  # Store model coefficients (fr5 from glmmTMB)
  # Store model coefficients for each Site:Year
  Site=rownames(ranef(fr5)$cond$'Site')
  Year=rownames(ranef(fr5)$cond$'Year')
  counter=0
  for (i in 1:length(Site)) {
    for (j in 1:length(Year)) {
      counter=counter+1
      fruit$fruit.int[counter]=fixef(fr5)$cond[1] + ranef(fr5)$cond$'Site'[i,1] + ranef(fr5)$cond$'Year'[j,1]
      fruit$fruit.slope[counter]=fixef(fr5)$cond[2] + ranef(fr5)$cond$'Site'[i,2] + ranef(fr5)$cond$'Year'[j,2]
      fruit$SiteYear[counter]=paste(rownames(ranef(fr5)$cond$'Site')[i],":",rownames(ranef(fr5)$cond$'Year')[j], sep="")
    }
  }
  
  # make data frame
  fruit=data.frame(fruit) 
  
  #*******************************************************************************
  ### 1E. Size distribution of recruits ###
  #*******************************************************************************
  
  # Calculate mean  and standard deviation of size in year t+1 of plants that were new that year (unrecorded at time t)
  recruit_size <- data %>% 
    filter(is.na(logSize)) %>% 
    group_by(SiteYear) %>% 
    summarise(recruit.size.mean = mean(logSizeNext),
              recruit.size.sd = sd(logSizeNext))
  
   #*******************************************************************************
  ### 1F. Create data frame of site-specific parameter estimates and join all estimates ###
  #*******************************************************************************
  
  params <- full_join(params,growth) %>% full_join(flower) %>% full_join(fruit) %>% full_join(recruit_size) %>% 
    replace_na(list(recruit.size.mean=0, recruit.size.sd=0)) # assign 0s to sites with no recruitment (mean & SD = NA) or only 1 recruit (mean = value but SD=NA)

  
  #*******************************************************************************
  ### 1G. Number of seeds per fruit ###
  #*******************************************************************************
  
  # Note: these are site-specific, but not site-by-year-specific, values, but we want to copy the site-specific values into each site-year level
  
  # obtain mean seed counts per fruit per site
  seeds.per.siteyear <- data %>% 
    group_by(SiteYear, Site) %>% 
    summarize(SeedCtSY = mean(SeedCt, na.rm=T)) # this yields NaN for some site-year combos

  seeds.per.site <- data %>% 
    group_by(Site) %>% 
    summarize(seed.ct = mean(SeedCt, na.rm=T)) 
  
  seeds.per.site <- left_join(seeds.per.siteyear, seeds.per.site) %>% dplyr::select(SiteYear, seed.ct)
  
  # join seeds per fruit to other parameters
  params <- left_join(params, seeds.per.site)
                      
  
  #*******************************************************************************
  ### 1H. Establishment probability ###
  #*******************************************************************************
  
  # Obtain number of new recruits per site at year = t+1
  recruit.number <- data %>% 
    dplyr::filter(is.na(logSize)) %>% 
    group_by(SiteYear) %>% 
    summarise(recruit.number = n()) 
    
  # Obtain total fruit count per site at year = t
  fruits.per.site <- site_fruit_count_data %>% 
    filter(!is.na(Fec1)) %>% 
    group_by(SiteYear) %>% 
    summarise(fruits.per.site = sum(Fec1))
  
  # Join recruit # and fruits per site into one frame
  # Obtain total seed count per site (= # fruits per site at year t * # seeds per fruit per site)
  # Obtain establishment probability (= # of new recruits at year t+1/total seed count per site at year t)
  establishment <- full_join(recruit.number, fruits.per.site) %>% full_join(seeds.per.site) %>% 
    replace_na(list(recruit.number=0, fruits.per.site=0)) %>% # assign 0s to sites with no recruitment or no fertility
    # NOTE: Amy added this in 2023. 
    mutate(total.seeds.per.site = fruits.per.site*seed.ct,
           establishment.prob = ifelse(recruit.number==0, 0, 
                                  ifelse(recruit.number>0 & fruits.per.site==0, recruit.number/seed.ct,recruit.number/total.seeds.per.site))) 
    # NOTE: there are a handful of site-years where recruit.number>0 but fruits.per.site=0, resulting in Inf for establishment.prob. Since establishment is non-zero in these site-years, we need a reasonable number for fruits.per.site. Solution used here is to replace fruits.per.site with value of 1 because fecundity was very low that year.
  
  # Join with params frame
    params <- left_join(params, establishment)

  # Add separate columns for year and site
  params <- params %>% separate(SiteYear, c("Site","Year"), ":", remove=FALSE)
  
  
#*******************************************************************************
### 2. Create site-specific IPMs parameterized by site-specific parameters derived from global vital rates models 
#*******************************************************************************

  #*******************************************************************************
  ### 2A. Subset data for site f
  #*******************************************************************************
  
  # Which site-years are missing parameter estimates?
  params.missing <- params[!complete.cases(params),]
  params.missing$SiteYear

  # Canton Creek:2011 --> growth slopes, intercepts, and sd are inestimable because of sparse data, but site was visited and censused and all other parameters are estimable; use mean across other years at this site for this year's growth parameters
  params$growth.int[params$SiteYear=="Canton Creek:2011"] = mean(params$growth.int[params$Site=="Canton Creek"], na.rm=T)
  params$growth.slope[params$SiteYear=="Canton Creek:2011"] = mean(params$growth.slope[params$Site=="Canton Creek"], na.rm=T)
  params$growth.sd[params$SiteYear=="Canton Creek:2011"] = mean(params$growth.sd[params$Site=="Canton Creek"], na.rm=T)
  # Canton Creek:2016 --> site inaccessible in 2017 due to high water; lambda is NA
  # Canton Creek 2017 --> site inaccessible in 2017 due to high water; lambda is NA
  
  # Rock Creek:2012 --> 2013 data folder lost; lambda is NA
  # Rock Creek:2013 --> 2013 data folder lost; lambda is NA
  
  # Buck Meadows:2012 --> site inaccessible in 2013 due to fire; lambda is NA
  # Buck Meadows:2013 --> site inaccessible in 2013 due to fire; lambda is NA
  
  # Rainbow Pool:2012 --> site inaccessible in 2013 due to fire; lambda is NA
  # Rainbow Pool:2013 --> site inaccessible in 2013 due to fire; lambda is NA
  
  # Carlon:2012 --> site inaccessible in 2013 due to fire; lambda is NA
  # Carlon:2013 --> site inaccessible in 2013 due to fire; lambda is NA
  # Carlon:2016 --> growth slopes, intercepts, and sd are inestimable because of sparse data, but site was visited and censused and all other parameters are estimable; use mean across other years at this site for this year's growth parameters
  params$growth.int[params$SiteYear=="Carlon:2016"] = mean(params$growth.int[params$Site=="Carlon"], na.rm=T)
  params$growth.slope[params$SiteYear=="Carlon:2016"] = mean(params$growth.slope[params$Site=="Carlon"], na.rm=T)
  params$growth.sd[params$SiteYear=="Carlon:2016"] = mean(params$growth.sd[params$Site=="Carlon"], na.rm=T)
  
  # Hauser Creek: 2011 --> all plants dead on plots 1-2, but many plants on plot 3 that were indistinguishable from one another, making growth parameters too sparse to estimate; use mean across other years
  params$growth.int[params$SiteYear=="Hauser Creek:2011"] = mean(params$growth.int[params$Site=="Hauser Creek"], na.rm=T)
  params$growth.slope[params$SiteYear=="Hauser Creek:2011"] = mean(params$growth.slope[params$Site=="Hauser Creek"], na.rm=T)
  params$growth.sd[params$SiteYear=="Hauser Creek:2011"] = mean(params$growth.sd[params$Site=="Hauser Creek"], na.rm=T)
  # Hauser Creek:2012 --> all plants dead; manually set lambda to 0 below
  # Hauser Creek:2013 --> all plants dead; manually set lambda to 0 below
  # Hauser Creek:2014 --> all plants dead; manually set lambda to 0 below
  # Hauser Creek:2015 --> all plants dead; manually set lambda to 0 below
  # Hauser Creek:2016 --> all plants dead; manually set lambda to 0 below
  # Hauser Creek:2017 --> all plants dead; manually set lambda to 0 below
  
  # Kitchen Creek: 2013 --> only remaining plant died; manually set lambda to 0 below
  # Kitchen Creek:2014 --> all plants dead; manually set lambda to 0 below
  # Kitchen Creek:2015 --> all plants dead; manually set lambda to 0 below
  # Kitchen Creek:2016 --> all plants dead; manually set lambda to 0 below
  # Kitchen Creek:2017 --> all plants dead; manually set lambda to 0 below
  # Kitchen Creek:2018 --> new recruits appeared but all other transitions inestimable because of no survivors; lambda is NA
  
  # Mill Creek:2010 --> all 2010 plots washed out and new plots established in 2011; lambda is NA
  # Mill Creek:2013 --> site inaccessible in 2013 due to flood; lambda is NA
  # Mill Creek:2014 --> site inaccessible in 2013 due to flood; lambda is NA
  
  # West Fork Mojave River:2013 --> existing plots 1-5 all dead, so some parameters inestimable. But new plot 6 established in 2014, so the entire site was not dead, only the main area where we were observing 2010-2013. Keep as NA because the entire site was not dead (in contrast to Hauser, Kitchen, Whitewater, where we set lambda to 0 when all plants died).

  # Whitewater Canyon:2014 --> all plants died; manually set lambda to 0 below
  # Whitewater Canyon:2015 --> site visited but no plants; manually set lambda to 0 below
  # Whitewater Canyon:2016 --> site visited but no plants; manually set lambda to 0 below
  # Whitewater Canyon:2017 --> site visited but no plants; manually set lambda to 0 below
  # Whitewater Canyon:2018 --> site visited but no plants; manually set lambda to 0 below

  # North Fork Middle Fork Tule:2016 --> fire closure in 2017; lambda is NA
  # North Fork Middle Fork Tule:2017 --> fire closure in 2017; lambda is NA
  
  # Redwood Creek:2015 --> fire closure in 2016; lambda is NA
  # Redwood Creek:2016 --> fire closure in 2016; lambda is NA
  
  # South Fork Middle Fork Tule:2010 --> fire closure in 2011; lambda is NA
  # South Fork Middle Fork Tule:2011 --> fire closure in 2011; lambda is NA
  # South Fork Middle Fork Tule:2016 --> fire closure in 2017; lambda is NA
  # South Fork Middle Fork Tule:2017 --> fire closure in 2017; lambda is NA
  

  
# Store parameters in .csv file for later use
write.csv(params,"data/demography data/vital_rate_coefficients.csv", row.names=FALSE)

# Remove rows of parameters with NA
  params <- params[complete.cases(params), ]
  
  
  # Remove site x year combinations without parameter estimates
  siteYear <- params$SiteYear
  
  # Create empty vectors for lambda and site to be filled
  lambda=c()
  SiteYear=character()
  
  for (f in 1:length(siteYear)) {
    data1 = subset(data, SiteYear==siteYear[f])
    params1 = subset(params, SiteYear==siteYear[f])
    params1 = subset(params1, select=-SiteYear)
    
    #*******************************************************************************
    ### 2B. Create survival, growth, and fecundity functions and build IPM by running integral_projection_model.R script
    #*******************************************************************************
    
    source("scripts/demography scripts/integral_projection_model.R")
    
    #*******************************************************************************
    ### 2C. Obtain lambda estimate for site f
    #*******************************************************************************
    
    lambda[f] <- Re(eigen(K)$values[1])
    SiteYear[f]=as.character(siteYear[f])
    } # end loop to run IPMs and estimate lambdas for each site
    
    # make data frame of site and lambda
    siteYear.lambda=data.frame(SiteYear,lambda)
    
#*******************************************************************************
### 3. Merge site information with lambda estimates and save to .csv file
#*******************************************************************************

# Create data frame of Site, Latitude, Longitude, Region, and Elevation 
site.info <- subset(data, select=c(Site,Year,SiteYear,Latitude,Longitude,Elevation,Region,RegionRank)) %>% unique() %>% arrange(-Latitude)
    
# Merge site info with lambda estimates
site.info <- full_join(site.info, siteYear.lambda)

# How many site-year combos are possible? 
20*9 #180
focal.sites <- c(rep("Coast Fork of Williamette",9),
                 rep("Canton Creek",9),
                 rep("Rock Creek",9),
                 rep("O'Neil Creek",9),
                 rep("Deep Creek",9),
                 rep("Little Jameson Creek",9),
                 rep("Oregon Creek",9),
                 rep("Rainbow Pool",9),
                 rep("Carlon",9),
                 rep("Buck Meadows",9),
                 rep("Wawona",9),
                 rep("Redwood Creek",9),
                 rep("North Fork Middle Fork Tule",9),
                 rep("South Fork Middle Fork Tule",9),
                 rep("West Fork Mojave River",9),
                 rep("Mill Creek",9),
                 rep("Whitewater Canyon",9),
                 rep("Sweetwater River",9),
                 rep("Kitchen Creek",9),
                 rep("Hauser Creek",9)) 
years <- rep(c("2010","2011","2012","2013","2014","2015","2016","2017","2018"), 20) 
site.years.max <- as.data.frame(cbind(focal.sites, years)) %>% mutate(SiteYear = paste(focal.sites,":", years, sep="")) %>% dplyr::select(SiteYear)

# How many site-years have some observed data but are missing lambda estimates?
length(unique(data$SiteYear)) #153
site.years.obs <- as.data.frame(unique(data$SiteYear))
colnames(site.years.obs) = "SiteYear"
   
# Which site-years are missing, and why?
missing.site.years <- anti_join(site.years.max, site.years.obs)

# Canton Creek:2016 --> site inaccessible in 2017 due to high water; this is a real NA
# Rock Creek:2012 --> 2013 data folder lost; this is a real NA 
# Rock Creek:2013 --> 2013 data folder lost; this is a real NA 
# Rainbow Pool:2012 --> site inaccessible in 2013 due to fire; this is a real NA
# Rainbow Pool:2013 --> site inaccessible in 2013 due to fire; this is a real NA
# Carlon:2012 --> site inaccessible in 2013 due to fire; this is a real NA
# Carlon:2013 --> site inaccessible in 2013 due to fire; this is a real NA
# Buck Meadows:2012 --> site inaccessible in 2013 due to fire; this is a real NA
# Buck Meadows:2013 --> site inaccessible in 2013 due to fire; this is a real NA
# Redwood Creek:2015 --> site inaccessible in 2016 due to fire; this is a real NA
# North Fork Middle Fork Tule: 2016 --> site inaccessible in 2017 due to fire; this is a real NA
# North Fork Middle Fork Tule: 2017 --> site inaccessible in 2017 due to fire; this is a real NA
# South Fork Middle Fork Tule:2010 --> site inaccessible in 2011 due to fire; this is a real NA
# South Fork Middle Fork Tule:2011 --> site inaccessible in 2011 due to fire; this is a real NA
# South Fork Middle Fork Tule:2016 --> site inaccessible in 2017 due to fire; this is a real NA
# South Fork Middle Fork Tule:2017 --> site inaccessible in 2017 due to fire; this is a real NA
# Mill Creek: 2013 --> site inaccessible in 2013 due to flood; this is a real NA
# Whitewater Canyon:2015 --> site visited but no plants so manually set lambda=0
Whitewater=filter(site.info, Site=="Whitewater Canyon" & Year==2014)
Whitewater$Year=2015 %>% factor()
Whitewater$SiteYear="Whitewater Canyon:2015" %>% factor()
Whitewater$lambda=0
site.info=bind_rows(site.info, Whitewater) 

# Kitchen Creek:2014 --> site visited but no plants so manually set lambda=0
Kitchen=filter(site.info, Site=="Kitchen Creek" & Year==2013)
Kitchen$Year=2014 %>% factor()
Kitchen$SiteYear="Kitchen Creek:2014" %>% factor()
Kitchen$lambda=0
site.info=bind_rows(site.info, Kitchen) 

# Kitchen Creek:2015 --> site visited but no plants so manually set lambda=0
Kitchen$Year=2015 %>% factor()
Kitchen$SiteYear="Kitchen Creek:2015" %>% factor()
Kitchen$lambda=0
site.info=bind_rows(site.info, Kitchen) 

# Kitchen Creek:2017 --> site visited but no plants so manually set lambda=0
Kitchen$Year=2017 %>% factor()
Kitchen$SiteYear="Kitchen Creek:2017" %>% factor()
Kitchen$lambda=0
site.info=bind_rows(site.info, Kitchen) 

# Hauser Creek:2012 --> site visited but no plants so manually set lambda=0
Hauser=filter(site.info, Site=="Hauser Creek" & Year==2011)
Hauser$Year=2012 %>% factor()
Hauser$SiteYear="Hauser Creek:2012" %>% factor()
Hauser$lambda=0
site.info=bind_rows(site.info, Hauser) 

# Hauser Creek:2013 --> site visited but no plants so manually set lambda=0
Hauser$Year=2013 %>% factor()
Hauser$SiteYear="Hauser Creek:2013" %>% factor()
Hauser$lambda=0
site.info=bind_rows(site.info, Hauser) 

# Hauser Creek:2014 --> site visited but no plants so manually set lambda=0
Hauser$Year=2014 %>% factor()
Hauser$SiteYear="Hauser Creek:2014" %>% factor()
Hauser$lambda=0
site.info=bind_rows(site.info, Hauser) 

# Hauser Creek:2015 --> site visited but no plants so manually set lambda=0
Hauser$Year=2015 %>% factor()
Hauser$SiteYear="Hauser Creek:2015" %>% factor()
Hauser$lambda=0
site.info=bind_rows(site.info, Hauser) 

# Hauser Creek:2016 --> site visited but no plants so manually set lambda=0
Hauser$Year=2016 %>% factor()
Hauser$SiteYear="Hauser Creek:2016" %>% factor()
Hauser$lambda=0
site.info=bind_rows(site.info, Hauser) 

# Hauser Creek:2017 --> site visited but no plants so manually set lambda=0
Hauser$Year=2017 %>% factor()
Hauser$SiteYear="Hauser Creek:2017" %>% factor()
Hauser$lambda=0
site.info=bind_rows(site.info, Hauser) 


# Which site-years have lambda=NA, and why?
lambda.calc.failed <- site.info %>% dplyr::filter(is.na(lambda)) %>% dplyr::select(SiteYear)

# Canton Creek:2017 --> site skipped due to high water; real NA
# Redwood Creek:2016 --> site skipped due to fire closure; real NA
# West Fork Mojave River:2013 --> existing plots 1-5 all dead, so some parameters inestimable. But new plot 6 established in 2014, so the entire site was not dead, only the main area where we were observing 2010-2013. Keep as NA because the entire site was not dead (in contrast to Hauser, Kitchen, Whitewater, where we set lambda to 0 when all plants died).
# Mill Creek:2010 --> all 2010 plots washed out and new plots established in 2011; lambda is NA
# Mill Creek:2014 --> site inaccessible in 2013 due to flood; lambda is NA
# Whitewater Canyon:2014 --> all plants died, so set lambda to 0 
site.info$lambda[site.info$SiteYear=="Whitewater Canyon:2014"] = 0
# Whitewater Canyon:2016 --> all plants died, so set lambda to 0 
site.info$lambda[site.info$SiteYear=="Whitewater Canyon:2016"] = 0
# Whitewater Canyon:2017 --> all plants died, so set lambda to 0 
site.info$lambda[site.info$SiteYear=="Whitewater Canyon:2017"] = 0
# Whitewater Canyon:2018 --> all plants died, so set lambda to 0 
site.info$lambda[site.info$SiteYear=="Whitewater Canyon:2018"] = 0
# Kitchen Creek: 2013 --> all plants died, so set lambda to 0
site.info$lambda[site.info$SiteYear=="Kitchen Creek:2013"] = 0
# Kitchen Creek: 2016 --> all plants died, so set lambda to 0
site.info$lambda[site.info$SiteYear=="Kitchen Creek:2016"] = 0
# Kitchen Creek: 2018 --> new recruits, but no survivors to estimate surv and growth functions; keep as NA

# View final data frame
str(site.info)

# Save to lambdas to .csv file 
write.csv(site.info,"data/demography data/siteYear.lambda_2010-2019.csv",row.names=FALSE)

# Merge params and lambdas for supplementary table
full.table <- left_join(site.info, params) %>% mutate_at(9:20, round, 2) %>% arrange(Latitude)
write.csv(full.table,"data/demography data/siteYear.paramslambda_2010-2019.csv",row.names=FALSE)

#*******************************************************************************
### 4. Compare to preliminary estimates 
#*******************************************************************************
# 
old.lambda <- read.csv("data/demography data/siteYear.lambda_2010-2016_old.csv") %>% 
  dplyr::select(SiteYear, lambda.old=lambda)
# Note: this file is copied in from demography_analyses repo

site.info <- site.info %>% 
  mutate(SiteYear = gsub(" ", "", SiteYear))

all <- left_join(site.info, old.lambda) %>% 
  mutate(lambda.diff = lambda - lambda.old,
         check = ifelse(abs(lambda.diff)>1, "YES", "no"),
         check.prop = ifelse(abs(lambda.diff/lambda)>0.50, "YES", "no"))

summary(all$lambda.diff)
table(all$check)
table(all$check.prop)

check.list = all$SiteYear[all$check=="YES"]

ggplot(data=all, aes(x=lambda, y=lambda.old, label=SiteYear)) +
  geom_point() + geom_text() + geom_abline(x=y) + xlim(0,70)
# Note: Buck Meadows:2015 and Sweetwater:2016 are absurd - what happened??
ggplot(data=all, aes(x=lambda, y=lambda.old, label=SiteYear)) +
  geom_point() + geom_text() + xlim(0,0.5) + ylim(0,0.5) + geom_abline(x=y)
# looks great - symmetrical and very close to 1:1
ggplot(data=all, aes(x=lambda, y=lambda.old, label=SiteYear)) +
  geom_point() + geom_text() + xlim(0.5,1) + ylim(0.5,1) + geom_abline(x=y)
# looks good - symmetrical and fairly close to 1:1
ggplot(data=all, aes(x=lambda, y=lambda.old, label=SiteYear)) +
  geom_point() + geom_text() + xlim(1,2) + ylim(1,2) + geom_abline(x=y)
# looks good - symmetrical and mostly close to 1:1
ggplot(data=all, aes(x=lambda, y=lambda.old, label=SiteYear)) +
  geom_point() + geom_text() + xlim(2,3) + ylim(2,3) + geom_abline(x=y)
ggplot(data=all, aes(x=lambda, y=lambda.old, label=SiteYear)) +
  geom_point() + geom_text() + xlim(3,10) + ylim(3,10) + geom_abline(x=y)
ggplot(data=all, aes(x=lambda, y=lambda.old, label=SiteYear)) +
  geom_point() + geom_text() + xlim(10,50) + ylim(10,50) + geom_abline(x=y)


