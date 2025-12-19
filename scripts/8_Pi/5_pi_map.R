#### PROJECT: Mimulus cardinalis evolutionary rescue
#### PURPOSE: Plot maps for M. cardinalis study population locations and nucleotide diversity (pi)
#### AUTHOR: Seema Sheth 
#### DATE LAST MODIFIED: 20241211


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
timeseries_pop <- all_pop %>% filter(Paper_ID<12) %>% dplyr::select(Long,Lat)
baseline_pop <- all_pop %>% filter(Paper_ID>11) %>% dplyr::select(Long,Lat)

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

# save genomics sampling map 
save_plot(filename = "Graphs/Maps/base_time_map_2B.png",
          plot = card_map_focal_pops,
          bg = "transparent", 
          base_width = 4, base_height = 8)
save_plot(filename = "Graphs/Maps/base_time_map_2B.pdf",
          plot = card_map_focal_pops,
          bg = "transparent", 
          base_width = 4, base_height = 8)



###################################################################################
### Read in nucleutide diversity estimates (Pi)
pi_df <- read_csv("data/genomic_data/baseline_pi_clump.csv")
pi_all_pop <-left_join(all_pop, pi_df,by=c("Paper_ID"="Site")) %>% filter(Paper_ID!=12)

############################################################################################ Plot pi on maps

# pi of climate SNP
pi_climate_map = ggplot(data=world,fill="lightgrey",col="black",size=0.3) + 
  geom_sf() +
  geom_polygon(aes(x = long, y = lat, group = group), data=states, fill="transparent", col="black",size=0.2) +
  theme_minimal() +
  coord_sf(xlim = c(xmin,xmax), ylim = c(ymin+0.5,ymax), expand = FALSE) +
  geom_point(data=pi_all_pop, aes(x=Long, y=Lat, fill=pi_snp_set), shape=21, size=4, alpha=0.8) +
  scale_fill_gradient(low="blanchedalmond", high="firebrick3") +
  #scale_fill_gradient(limits = range(pi_all_pop$pi_snp_set, pi_all_pop$pi_all_snps)) +
  labs(x="Longitude",y="Latitude") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=16),legend.position = c(.93, .79),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        plot.title=element_text(hjust=0,size=18),
        legend.title = element_blank()) 
pi_climate_map

# save pi map
save_plot(filename = "Graphs/Maps/pi_climate_snps_map_S5A.png",
          plot = pi_climate_map,
          bg = "transparent", 
          base_width = 4, base_height = 8)
save_plot(filename = "Graphs/Maps/pi_climate_snps_mapS5A.pdf",
          plot = pi_climate_map,
          bg = "transparent", 
          base_width = 4, base_height = 8)

# pi of all SNP
pi_global_map = ggplot(data=world,fill="lightgrey",col="black",size=0.3) + 
  geom_sf() +
  geom_polygon(aes(x = long, y = lat, group = group), data=states, fill="transparent", col="black",size=0.2) +
  theme_minimal() +
  coord_sf(xlim = c(xmin,xmax), ylim = c(ymin+0.5,ymax), expand = FALSE) +
  geom_point(data=pi_all_pop, aes(x=Long, y=Lat, fill=pi_all_snps), shape=21, size=4, alpha=0.8) +
  scale_fill_gradient(low="white", high="grey30") +
  labs(x="Longitude",y="Latitude") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=16),legend.position = c(.95, .79),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        plot.title=element_text(hjust=0,size=18),
        legend.title = element_blank()) 
pi_global_map

# save pi map
save_plot(filename = "Graphs/Maps/pi_all_snps_map_S5B.png",
          plot = pi_global_map,
          bg = "transparent", 
          base_width = 4, base_height = 8)
save_plot(filename = "Graphs/Maps/pi_all_snps_map_S5B.pdf",
          plot = pi_global_map,
          bg = "transparent", 
          base_width = 4, base_height = 8)

