##################################################################################
#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Make table of Geometric Means
#### AUTHOR: Daniel Anstett and Amy Angert
#### DATE LAST MODIFIED: 20250217
###################################################################################

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

#Library install and import
library(tidyverse) 

#Import demography data (includes metadata)
demog_recovery <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv")

geo_means <- demog_recovery %>% select(Site,Paper_ID,Latitude,Longitude,mean.lambda.pre,mean.lambda.drought,mean.lambda.recovery) %>%
  arrange(Latitude)

write_csv(geo_means,"data/demography data/geo_means.csv")

# add color ramp column
color <- read_csv("data/genomic_data/pop_meta_data.csv") %>%
  dplyr::select(Site, Lat.Color) %>%
  unique()

# pivot longer for plotting
geo_means_long <- left_join(geo_means, color) %>% pivot_longer(cols=mean.lambda.pre:mean.lambda.recovery, names_to="time", values_to="lambda")

ggplot(geo_means_long, aes(x=time,y=lambda)) +
  geom_point(aes(fill=as.factor(Latitude)), shape=21, size=2) +
  geom_line(aes(group=Latitude), color=geo_means_long$Lat.Color) +
  scale_color_manual(values=unique(geo_means_long$Lat.Color), aesthetics = c("color", "fill")) +
  geom_hline(yintercept=1) +
  theme_classic()

# pre to drought
t.test(x=geo_means$mean.lambda.pre, y=geo_means$mean.lambda.drought, paired=TRUE)

# drought to recovery
t.test(x=geo_means$mean.lambda.recovery, y=geo_means$mean.lambda.drought, paired=TRUE)


