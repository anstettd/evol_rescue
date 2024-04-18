#### PROJECT: Genomic offsets and demographic trajectories of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Tidy raw demographic vital rate data in preparation for IPM analyses 
#### AUTHOR: Seema Sheth and Amy Angert
#### DATE LAST MODIFIED: 20230619

#*******************************************************************************
#### 1. Clean workspace and load required packages
#*******************************************************************************

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

# Easy code for installing packages in R (if not installed) and calling their libraries
# From: https://gist.github.com/DrK-Lo/a945a29d6606b899022d0f03109b9483

# Make vector of packages needed
packages_needed <- c("tidyverse")

# Install packages needed (if not already installed)
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
}

# Load packages needed
for (i in 1:length(packages_needed)){
  library( packages_needed[i], character.only = TRUE)
}

#*******************************************************************************
#### 2. Bring in and combine M. cardinalis vital rate data from 2010-2019
#*******************************************************************************

# Read in seed count per fruit data, select and rename relevant columns, and round to nearest integer
seed.ct <- read.csv("data/demography data/fec2.seed.per.fruit.2010.2011.2012.csv")
seed.ct <- subset(seed.ct,select=c(site,newgrandmean))
colnames(seed.ct) = c("Site","SeedCt")
seed.ct$SeedCt = round(seed.ct$SeedCt,digits=0) # round seed count to nearest integer
seed.ct$Site = factor(seed.ct$Site) # make site column a factor to streamline joining
# Note: could try to update with more recent years' collections if available

# Read in vital rate data for 2010-11, 2011-12, 2012-13, 2013-14 transitions 
# Note: PY=previous year (time t), CY=current year (time t+1); ignore PPY
data_2010.2014=read.csv("data/demography data/Mcard_demog_data_2010-2014.csv") %>% dplyr::select(-Reasoning, -Reasoning.1) %>% 
  mutate(PlotID=NA) #remove unwanted columns and add plot column for merging
# Note: this file was created for the analyses published in Sheth and Angert 2018 PNAS 
# It results from Amy Angert's work in July 2016 (original file: "Mcard_demog_data_2010-2013_ALA.xlsx") to scan datasheet notes to identify individuals to exclude, based on these columns:
# Column 'NotAnIndividual': 
#   0 = ok, definitely include (includes "scattered" as long as nothing else noted)
#   1 = not ok, definitely exclude from survival, growth, and fecundity but ok for seed input denominator for recruitment (history of lumping/splitting/relumping; redundant IDs)
#   2 = maybe (notes about difficult to distinguish, or merged once)
# Column 'NotARecruit':
#   0 = ok, definitely include
#   1 = wrong, definitely exclude (reasons include new plot, site not visited in prior year, ID within prior years' ranges, coordinates well outside of prior year's search
#   2 = plant noted as possibly missed in prior year (looks old, missed?, J-20xx?, could be [old ID], etc)
#   3 = plant not noted as unusually old-looking, but is within size range of or larger than most plants that are in category 2
#   NA = size measures at year t

# Read in vital rate data for 2014-15 transition 
# Note: PY=previous year (time t), CY=current year (time t+1); ignore PPY
data_2014.2015=read.csv("data/demography data/SS_Horizontal_2015_Notes.csv") %>% 
  rename(Site = SiteID, #rename to match 2010-2014 data
         Class = Class_PY,
         ClassNext = Class_CY,
         Fec1 = TotFr_PY,
         Size = TotStLn_PY,
         SizeNext = TotStLn_CY,
         Surv = SurvPYCY) %>% 
  mutate(Year = 2014) #add year at time t (=2014) column
# Note: this file was queried from the database on 2023-02-01
# We will tidy it below based on the decision rules as above for 2010-14 data, but scripted here instead of done manually in excel. 

# Read in vital rate data for 2015-16 transition
# Note: PY=previous year (time t), CY=current year (time t+1); ignore PPY
data_2015.2016 = read.csv("data/demography data/SS_Horizontal_2016_Notes.csv") %>% 
  rename(Site = SiteID, #rename to match 2010-2014 data
         Class = Class_PY,
         ClassNext = Class_CY,
         Fec1 = TotFr_PY,
         Size = TotStLn_PY,
         SizeNext = TotStLn_CY,
         Surv = SurvPYCY) %>% 
  mutate(Year = 2015) #add year at time t (=2015) column)
# Note: this file was queried from the database on 2023-02-01
# We script the creation of NotARecruit and NotAnIndividual columns below (instead of  manually in excel)

# Read in vital rate data for 2016-17 transition
# Note: PY=previous year (time t), CY=current year (time t+1); ignore PPY
data_2016.2017 = read.csv("data/demography data/SS_Horizontal_2017_Notes.csv") %>% 
  rename(Site = SiteID, #rename to match 2010-2014 data
         Class = Class_PY,
         ClassNext = Class_CY,
         Fec1 = TotFr_PY,
         Size = TotStLn_PY,
         SizeNext = TotStLn_CY,
         Surv = SurvPYCY) %>% 
  mutate(Year = 2016) #add year at time t (=2016) column)
# Note: this file was queried from the database on 2023-06-19
# We script the creation of NotARecruit and NotAnIndividual columns below (instead of  manually in excel)

# Read in vital rate data for 2017-18 transition
# Note: PY=previous year (time t), CY=current year (time t+1); ignore PPY
data_2017.2018 = read.csv("data/demography data/SS_Horizontal_2018_Notes.csv") %>% 
  rename(Site = SiteID, #rename to match 2010-2014 data
         Class = Class_PY,
         ClassNext = Class_CY,
         Fec1 = TotFr_PY,
         Size = TotStLn_PY,
         SizeNext = TotStLn_CY,
         Surv = SurvPYCY) %>% 
  mutate(Year = 2017) #add year at time t (=2017) column)
# Note: this file was queried from the database on 2023-06-19
# We script the creation of NotARecruit and NotAnIndividual columns below (instead of  manually in excel)

# Read in vital rate data for 2018-19 transition
# Note: PY=previous year (time t), CY=current year (time t+1); ignore PPY
data_2018.2019 = read.csv("data/demography data/SS_Horizontal_2019_Notes.csv") %>% 
  rename(Site = SiteID, #rename to match 2010-2014 data
         Class = Class_PY,
         ClassNext = Class_CY,
         Fec1 = TotFr_PY,
         Size = TotStLn_PY,
         SizeNext = TotStLn_CY,
         Surv = SurvPYCY) %>% 
  mutate(Year = 2018) #add year at time t (=2018) column)
# Note: this file was queried from the database on 2023-06-19
# We script the creation of NotARecruit and NotAnIndividual columns below (instead of  manually in excel)

# Combine 2014-19 data into one data frame 
data_2014.2019 = rbind(data_2014.2015,data_2015.2016,data_2016.2017,data_2017.2018,data_2018.2019)

# Add 'NotAnIndividual' column
data_2014.2019 <- data_2014.2019 %>% 
  mutate(NotAnIndividual = ifelse(str_detect(OtherNotesCY, "lump"), 1, 
                           ifelse(str_detect(OtherNotesCY, "Lump"), 1,
                           ifelse(str_detect(OtherNotesCY, "split"), 1,
                           ifelse(str_detect(OtherNotesCY, "Split"), 1,
                           ifelse(str_detect(OtherNotesCY, "same as"), 1,
                           ifelse(str_detect(OtherNotesCY, "merge"), 1,
                           ifelse(str_detect(OtherNotesCY, "Merge"), 1,
                           ifelse(str_detect(OtherNotesPY, "part of"), 1, 
                           ifelse(str_detect(OtherNotesPY, "lump"), 1, 
                           ifelse(str_detect(OtherNotesPY, "Lump"), 1,
                           ifelse(str_detect(OtherNotesPY, "split"), 1,
                           ifelse(str_detect(OtherNotesPY, "Split"), 1,
                           ifelse(str_detect(OtherNotesPY, "same as"), 1,
                           ifelse(str_detect(OtherNotesPY, "merge"), 1,
                           ifelse(str_detect(OtherNotesPY, "Merge"), 1,
                           ifelse(str_detect(OtherNotesPY, "part of"), 1,       
                                0))))))))))))))))) 
# Note: this is lacking level 2 (=maybe), so is possibly more restrictive than 2010-14 filter
# TO DO: Repeat this automatic coding for 2010-2014 as a sensitivity analysis

# Read in list of skipped plots and sites in certain years. These cannot be used for recruitment the following year, because these plots are missing from the seed input denominator at time t, and at time t+1 we cannot distinguish 1-yr-old and 2-yr-old plants 
skipped <- read.csv("data/demography data/SkippedPlots.csv") %>% 
  separate(Squawk, sep=" ", c("Status", NA, "Year")) 
skipped$Year = as.numeric(skipped$Year)
### TO DO: ADD 2018

data_2014.2019 <- left_join(data_2014.2019, skipped, by=c("Site"="Site","PlotID"="PlotID","Year"="Year")) 
# Note: This file results from the "Skipped In" query, to which Amy added several plots after consulting field notes for 2014-2019
# For example, Deer Creek plot 4 line 1 was totally reset in 2015. 2014 plants were recorded as dead, but some of them could have survived and not been mapped onto new coordinate system. It looks as if this site had 100% mortality and high recruitment of large individuals, but in reality there are survivors that could not be identified and translated from old to new coordinate system. Exclude entire line for 2014-15 transitions and 2015 recruitment.

# Add 'NotARecruit' column
data_2014.2019 <- data_2014.2019 %>% 
  mutate(NotARecruit = ifelse(!is.na(Size), NA, 
                       ifelse(str_detect(OtherNotesCY, "old"), 2, 
                       ifelse(str_detect(OtherNotesCY, "Old"), 2,
                       ifelse(str_detect(OtherNotesCY, "missed"), 2,
                       ifelse(str_detect(OtherNotesCY, "14?"), 2, 
                       ifelse(str_detect(OtherNotesCY, "15?"), 2, 
                       ifelse(str_detect(OtherNotesCY, "16?"), 2, 
                       ifelse(str_detect(OtherNotesCY, "17?"), 2, 
                       ifelse(str_detect(OtherNotesCY, "18?"), 2, 
                       ifelse(str_detect(OtherNotesCY, "19?"), 2, 
                       ifelse(NewPlot_CY==TRUE, 1,
                       ifelse(Status=="Skipped", 1, 0)))))))))))))
# Note: this is lacking level 3 (=size range of other recruits), which is not reliable


# Create columns of log-transformed sizes
# first, replace with NA observations that erroneously have size=0 
data_2014.2019[68460,"Size"] <- NA
data_2014.2019[55399,"SizeNext"] <- NA
data_2014.2019[65542,"SizeNext"] <- NA
data_2014.2019[65547,"SizeNext"] <- NA
data_2014.2019$logSize = log(data_2014.2019$Size)
data_2014.2019$logSizeNext = log(data_2014.2019$SizeNext)

# Add a column ranking regions from south to north
data_2014.2019$RegionRank[data_2014.2019$Region=="S1"]=1
data_2014.2019$RegionRank[data_2014.2019$Region=="S2"]=2
data_2014.2019$RegionRank[data_2014.2019$Region=="C1"]=3
data_2014.2019$RegionRank[data_2014.2019$Region=="C2"]=4
data_2014.2019$RegionRank[data_2014.2019$Region=="C3"]=5
data_2014.2019$RegionRank[data_2014.2019$Region=="N1"]=6
data_2014.2019$RegionRank[data_2014.2019$Region=="N2"]=7

# Create column for probability of flowering
data_2014.2019$Fec0 = (ifelse(data_2014.2019$Class=="A", 1, 
                              ifelse(data_2014.2019$Class=="J",0 , NA))) 

# Merge transition data with seed count data
data_2014.2019 <- merge(data_2014.2019,seed.ct,by="Site",all.x=TRUE,all.y=FALSE)

data_2014.2019 <- data_2014.2019 %>% dplyr::select(colnames(data_2010.2014))

# Combine all years
data <- rbind(data_2010.2014, data_2014.2019)

# Make site x year variable
data$SiteYear = paste(data$Site, data$Year, sep=":") %>% factor()

#*******************************************************************************
#### 3. Remove unwanted data
#*******************************************************************************

# Remove sites that were not monitored after 2014 and one site with problematic/missing data (Deer Creek, where frequent plot washouts cause identifications to be unreliable)
# Retain focal 20 populations
focal.sites <- c("Coast Fork of Williamette",
                 "Canton Creek",
                 "Rock Creek",
                 "O'Neil Creek",
                 "Deep Creek",
                 "Little Jameson Creek",
                 "Oregon Creek",
                 "Rainbow Pool",
                 "Carlon",
                 "Buck Meadows",
                 "Wawona",
                 "Redwood Creek",
                 "North Fork Middle Fork Tule",
                 "South Fork Middle Fork Tule",
                 "West Fork Mojave River",
                 "Mill Creek",
                 "Whitewater Canyon",
                 "Sweetwater River",
                 "Kitchen Creek",
                 "Hauser Creek")

data <- data %>% filter(Site %in% focal.sites) %>% droplevels()

# Remove plants that were dead at time t (Class=="D") 
data <- subset(data, Class!="D" | is.na(Class))

# Remove plants with Class=? or Class=excluded
data <- subset(data, Class!="E" & Class!="?" | is.na(Class))

# Double check that Deer Creek 2013 Plot 4 Line 1, where existing plants were all marked D in 2014 but should have been marked E, are not in the data frame
# deer2013P04 <- subset(data, Site=="Deer Creek" & Year==2013 & PlotID==238) #no rows

# Remove rows for which size at time t AND size at t+1 is NA
data <- data[!(is.na(data$logSize) & is.na(data$logSizeNext)),]

# Remove plants that were recorded as having a class at time t but have no size measurements in that year
data <- subset(data, !(!is.na(Class) & is.na(logSize)))

# Remove plants that were recorded as having a class at time t but have no survival recorded from t to t+1 (were either excluded at time t+1, recorded as "?" in ClassNext field, or recorded as "NA" in ClassNext field) 
data <- subset(data, !(!is.na(Class) & is.na(Surv)))

# Make Fec1 numeric and round fruit # to nearest integer
data$Fec1 <- round(as.numeric(data$Fec1, digits=0)) 

# Only include seed counts for plants that produced at least one fruit
data$SeedCt[data$Fec1 < 1 | is.na(data$Fec1)]=NA

# For most vital rate transitions, remove monster plants where individuals were not distinguishable 
# Note: these should still be included when calculating seed input denominator for recruitment probability
data.indivs = subset(data, NotAnIndividual != 1 | is.na(NotAnIndividual))

# For estimating size distribution of new recruits, also remove individuals that should not be scored as new recruits
data.indivs = subset(data.indivs, NotARecruit != 1 | is.na(NotARecruit))
# Note: because these rows are NA at time t, their exclusion has no bearing on survival, growth, and fecundity transitions. Confirmed:
summary(data.indivs$logSize[data.indivs$NotARecruit==1])


#*******************************************************************************
#### 4. Prepare data for IPMs and write to new .csv files
#*******************************************************************************

# Obtain total fruit and seed counts for each individual at each site in each year, including monster plants
site_fruit_count_data = subset(data, select=c(Site,Year,SiteYear,Region,Fec1,SeedCt)) 

# Examine column names and classes of data
names(data.indivs)
str(data.indivs)

# Check site names to be sure no discrepancies
unique(site_fruit_count_data$Site)
unique(data.indivs$Site)

# Check classes to be sure no discrepancies
unique(data.indivs$Class)
unique(data.indivs$ClassNext)

# Make appropriate columns of data frame a factor
site_fruit_count_data$Site=factor(site_fruit_count_data$Site)
site_fruit_count_data$Year=factor(site_fruit_count_data$Year)
site_fruit_count_data$Region=factor(site_fruit_count_data$Region)
data.indivs$Site=factor(data.indivs$Site)
data.indivs$Year=factor(data.indivs$Year)
data.indivs$Region=factor(data.indivs$Region)
data.indivs$ID=factor(data.indivs$ID)
data.indivs$NotARecruit=factor(data.indivs$NotARecruit)
data.indivs$NotAnIndividual=factor(data.indivs$NotAnIndividual)

# Examine data
names(data.indivs)
head(data.indivs)
tail(data.indivs)

# Sort data by latitude
data.indivs=data.indivs[order(-data.indivs$Latitude,data.indivs$Year),]

# Summarize for methods
length(unique(data.indivs$Site))
length(unique(data.indivs$ID))
site.n <- data.indivs %>% 
  group_by(Site) %>% 
  summarise(n = n_distinct(ID))
site.n
min(site.n$n)
max(site.n$n)
mean(site.n$n)

# write to .csv
write.csv(data.indivs,"data/demography data/Mcard_demog_data_2010-2019_cleanindivs.csv",row.names=FALSE)
write.csv(site_fruit_count_data,"data/demography data/Mcard_demog_data_2010-2019_seedinput.csv",row.names=FALSE)


#*******************************************************************************
#### 5. Description of columns in _cleanindivs .csv file
#*******************************************************************************

# Site: population
# ID: unique identifier for each individual
# Region: latitudinal region that population is nested within
# Latitude: latitude of population
# Longitude: longitude of population
# Elevation: elevation of population in meters
# Class: stage class (juvenile, adult, or NA) of plant at time = t (PY)
# Fec1: TotFr (Total number of fruits per individual, rounded to nearest integer)   
# logSize: total stem length of the individual, in log-transformed cm
# ClassNext: stage class (juvenile, adult, dead, or NA) of plant at time = t+1 (CY)
# logSizeNext: same as "size" above, for t+1
# Surv: survival (1) or not (0) of individuals between time = t (PY) and time = t+1 (CY)
# Year: annual transition of the long-term data at time = t (2010-2015)
# Fec0: Probability of flowering (1 if Class=="A" for adult, 0 if Class=="J" for juvenile)
# RegionRank: ordinal rank of regions from south to north
# SeedCt: mean seed count for each site
# NotAnIndividual: see above
# NotARecruit: see above
