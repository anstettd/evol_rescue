#############################################################################################################
## Get climate NA files for climate anomaly
## Author Daniel Anstett
## 
## 
## Last Modified June 17, 2023
#############################################################################################################
#Import libraries
#install.packages('C:/Users/anstett3/Documents/Genomics/Mimulus_offset/ClimateNAr.zip', repos=NULL, type='source')

library(ClimateNAr)
library(tidyverse)

#Import
demo_pop <- read_csv("data/climate_data/demography_pop.csv")

demo_clim_1979 <- ClimateNA_API2(ClimateBC_NA='NA', inputFile=demo_pop,period='Year_1979',MSY='M')



