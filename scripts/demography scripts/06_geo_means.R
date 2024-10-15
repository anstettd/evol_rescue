##################################################################################
#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Make table of Geometric Means
#### AUTHOR: Daniel Anstett and Amy Angert
#### DATE LAST MODIFIED: 20241015
###################################################################################

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

#Library install and import
library(tidyverse) 

#Import demography data (includes metadata)
demog_recovery <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv")

geo_means <- demog_recovery %>% select(Site,Paper_ID,Latitude,Longitude,mean.lambda.drought,mean.lambda.recovery) %>%
  arrange(Latitude)

write_csv(geo_means,"data/demography data/geo_means.csv")
