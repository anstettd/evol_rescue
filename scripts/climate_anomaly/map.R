#### PROJECT: Mimulus cardinalis evolutionary rescue
#### PURPOSE: Plot maps for M. cardinalis study population locations and nucleotide diversity (pi)
#### AUTHOR: Seema Sheth 
#### DATE LAST MODIFIED: 20240906


###################################################################################
### Install and load packages 
# Easy code for installing packages in R (if not installed) and calling their libraries
# From: https://gist.github.com/DrK-Lo/a945a29d6606b899022d0f03109b9483

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

# make vector of packages needed
packages_needed <- c("tidyverse","cowplot","RColorBrewer","ggmap","maps","mapdata","mapproj","raster","rnaturalearth","rnaturalearthdata")

# install packages needed (if not already installed)
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
}

# load packages needed
for (i in 1:length(packages_needed)){
  library( packages_needed[i], character.only = TRUE)
}

###################################################################################
### Read in population locations from Baseline & Time Series datasets

# Pops used in genome sequencing
all_pop <- read_csv("data/genomic_data/Baseline_Timeseries_pops_final2.csv") %>% filter(Paper_ID<56)
#all_pop_sf <- st_as_sf(all_pop,coords=c("Long","Lat"), crs=EPSG4326)
timeseries_pop <- all_pop %>% filter(Paper_ID<12) %>% dplyr::select(Long,Lat)
#timeseries_pop_sf <- st_as_sf(timeseries_pop,coords=c("Long","Lat"), crs=EPSG4326)
baseline_pop <- all_pop %>% filter(Paper_ID>11) %>% dplyr::select(Long,Lat)
#baseline_pop_sf <- st_as_sf(baseline_pop,coords=c("Long","Lat"), crs=EPSG4326)


###################################################################################
### Map of focal populations

# get world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# create map of states
states <- map_data("state")

# make bounding box
ymin = min(all_pop$Lat) - 1.2
ymax = max(all_pop$Lat) + 1.2
xmin = min(all_pop$Long) - 1.2 
xmax = max(all_pop$Long) + 1.2 
#e = extent(xmin, xmax,ymin, ymax)

# create a map of focal populations and occurrences without a legend  
card_map_focal_pops = ggplot(data=world,fill="lightgrey",col="black",size=0.3) + 
  geom_sf() +
  geom_polygon(aes(x = long, y = lat, group = group), data=states,fill="transparent",col="black",size=0.2) +
  theme_minimal() +
  coord_sf(xlim = c(xmin,xmax), ylim = c(ymin+0.5,ymax), expand = FALSE) +
  geom_point(aes(x=Long,y=Lat),data=baseline_pop,shape=4,size=3) +
  geom_point(aes(x=Long,y=Lat),data=timeseries_pop,shape=24,size=4, fill="black", color="black", alpha=0.6) +
  # scale_fill_manual(values=c("#9E0142","red2","goldenrod1","gold","dodgerblue4","#5E4FA2"), guide = "none") +
  #  scale_fill_manual(values=c("#9E0142","red2","goldenrod1","gold","dodgerblue4","#5E4FA2"), guide = guide_legend(reverse = TRUE)) +
  labs(x="Longitude",y="Latitude") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=16),legend.position = c(.99, .80),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        plot.title=element_text(hjust=0,size=18),
        legend.title = element_blank()) 

card_map_focal_pops

# save map 
save_plot(filename = "Graphs/Maps/base_time_map.png",
          plot = card_map_focal_pops,
          bg = "transparent", 
          base_width = 4, base_height = 8)



###################################################################################
# Demography
demo_pop <- read_csv("data/genomic_data/pop_meta_data.csv")
#demo_pop_sf <- st_as_sf(demo_pop ,coords=c("Longitude","Latitude"), crs=EPSG4326)

# Pi
pi_df <- read_csv("data/genomic_data/baseline_pi.csv")
pi_all_pop <-left_join(all_pop,pi_df,by=c("Paper_ID"="Site"))
#pi_all_pop_sf <- st_as_sf(pi_all_pop,coords=c("Long","Lat"), crs=EPSG4326)

# Define CRS
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS

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
  tm_layout(frame = FALSE,legend.position = c(0.59, 0.42),legend.title.size = 0.001,legend.text.size=1.28)
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
  tm_layout(frame = FALSE,legend.position = c(0.59, 0.5),legend.title.size = 0.001,legend.text.size=1.28)
mim_base
tmap_save(mim_base, filename = "Graphs/Maps/pi_global.png",width=4, height=7)
