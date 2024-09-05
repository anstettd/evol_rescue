##################################################################################
## Make population mapes of baseline and timeseires
## Author Daniel Anstett
## 
## 
## Last Modified Jan 20, 2023
###################################################################################

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

#Library install and import
library(tidyverse) # for data manipulation
library(raster) # for handeling rasters
library(tidyverse)
library(cowplot)
library(sf)
library(tmap)
library(rnaturalearth)
library(rnaturalearthdata)
#library(devtools)
#devtools::install_github("ropensci/rnaturalearthhires")
library(rnaturalearthhires)
library(rgeos)
library(RColorBrewer)
###################################################################################
#Define CRS
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS

#Baseline & Timeseries
all_pop <- read_csv("data/genomic_data/Baseline_Timeseries_pops_final2.csv")%>% filter(Paper_ID<56)
all_pop_sf <- st_as_sf(all_pop,coords=c("Long","Lat"), crs=EPSG4326)
timeseries_pop <- all_pop %>% filter(Paper_ID<13) %>% dplyr::select(Long,Lat)
timeseries_pop_sf <- st_as_sf(timeseries_pop,coords=c("Long","Lat"), crs=EPSG4326)
baseline_pop <- all_pop %>% filter(Paper_ID>12) %>% dplyr::select(Long,Lat)
baseline_pop_sf <- st_as_sf(baseline_pop,coords=c("Long","Lat"), crs=EPSG4326)

#Demography
demo_pop <- read_csv("data/genomic_data/pop_meta_data.csv")
demo_pop_sf <- st_as_sf(demo_pop ,coords=c("Longitude","Latitude"), crs=EPSG4326)

#Pi
pi_df <- read_csv("data/genomic_data/baseline_pi.csv")
pi_all_pop <-left_join(all_pop,pi_df,by=c("Paper_ID"="Site"))
pi_all_pop_sf <- st_as_sf(pi_all_pop,coords=c("Long","Lat"), crs=EPSG4326)

# California & Oregon Map Setup
states<-ne_states(country=c("canada","united states of america"),returnclass= "sf")
calo <- states %>%
  filter(name_en=="Oregon" | name_en=="California") #| name_en=="Nevada")

###################################################################################

#Baseline + Timeseries
tmap_mode("plot")
#tmap_mode("view")
mim <-
  tm_shape(calo)+
  tm_borders()+
  tm_shape(baseline_pop_sf)+
  tm_dots(size=0.4,shape=21,col="blue",border.col ="black")+
  tm_shape(timeseries_pop_sf)+
  tm_dots(size=0.6,shape=21,col="red",border.col ="black")+
  tm_layout(legend.position = c(1.03, 0.73),legend.title.size = 0.001,frame.lwd = NA)
mim
tmap_save(mim, filename = "Graphs/Maps/base_time.png",width=4, height=7)


#Baseline + Timeseries, Greyscale
tmap_mode("plot")
#tmap_mode("view")
mim <-
  tm_shape(calo)+
  tm_borders()+
  tm_shape(baseline_pop_sf)+
  tm_dots(size=0.6,shape=4)+
  tm_shape(timeseries_pop_sf)+
  tm_dots(size=0.8,shape=24)+
  tm_layout(legend.position = c(1.03, 0.73),legend.title.size = 0.001,frame.lwd = NA)
mim
tmap_save(mim, filename = "Graphs/Maps/base_time_grey.png",width=4, height=7)




#Plot Baseline
tmap_mode("plot")
#tmap_mode("view")
mim_base <-
  tm_shape(calo)+
  tm_borders()+
  tm_shape(all_pop_sf)+
  tm_dots(size=0.6,shape=21,col="red2",border.col ="black")+
  tm_layout(legend.position = c(1.03, 0.73),legend.title.size = 0.001,frame.lwd = NA)
mim_base
tmap_save(mim_base, filename = "Graphs/Maps/baseline.png",width=4, height=7)



#Plot Timeseries
tmap_mode("plot")
#tmap_mode("view")
mim_time <-
  tm_shape(calo)+
  tm_borders()+
  tm_shape(timeseries_pop_sf)+
  tm_dots(size=0.9,shape=21,col="blue",border.col ="black")+
  tm_layout(legend.position = c(1.03, 0.73),legend.title.size = 0.001,frame.lwd = NA)
mim_time
tmap_save(mim_time, filename = "Graphs/Maps/timeseries.png",width=4, height=7)


#Plot Demography
tmap_mode("plot")
#tmap_mode("view")
mim_time <-
  tm_shape(calo)+
  tm_borders()+
  tm_shape(demo_pop_sf)+
  tm_dots(size=0.9,shape=21,col="magenta2",border.col ="black")+
  tm_layout(legend.position = c(1.03, 0.73),legend.title.size = 0.001,frame.lwd = NA)
mim_time
#tmap_save(mim_time, filename = "Graphs/Maps/demography.png",width=4, height=7)


##############################################################################################
#Plot pi on maps

snp_pallet <- c("#2166AC","#67A9CF","#D1E5F0","#FDDBC7","#EF8A62","#B2182B")

#Plot SNP set pi
tmap_mode("plot")
#tmap_mode("view")
mim_base <-
  tm_shape(calo)+
  tm_borders()+
  tm_shape(pi_all_pop_sf)+
  tm_dots(size=0.4,shape=21, col="pi_snp_set",palette = snp_pallet, border.col="black" )+
  tm_layout(frame = FALSE,legend.position = c(0.59, 0.39),legend.title.size = 0.001,legend.text.size=1.28)
mim_base
tmap_save(mim_base, filename = "Graphs/Maps/pi_snp_set.png",width=4, height=7)


# Load the RdBu palette

global_pallet <- c("#08417b","#2166AC","#67A9CF","#D1E5F0")

#Plot Global pi
tmap_mode("plot")
#tmap_mode("view")
mim_base <-
  tm_shape(calo)+
  tm_borders()+
  tm_shape(pi_all_pop_sf)+
  tm_dots(size=0.4,shape=21, col="pi_all_snps",palette = global_pallet, border.col="black" )+
  tm_layout(frame = FALSE,legend.position = c(0.59, 0.45),legend.title.size = 0.001,legend.text.size=1.28)
mim_base
tmap_save(mim_base, filename = "Graphs/Maps/pi_global.png",width=4, height=7)
