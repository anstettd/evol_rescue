##################################################################################
## Daniel Anstett
## Maps showing 19 demography populations
## 
## Last Modified January 22, 2020
###################################################################################
# Remove objects and clear workspace
rm(list = ls(all=TRUE))

# make vector of packages needed
packages_needed <- c("tidyverse","cowplot","RColorBrewer","ggmap","maps","mapdata","mapproj","raster","rnaturalearth","rnaturalearthdata")
library(ggrepel)
library(sf)
# load packages needed
for (i in 1:length(packages_needed)){
  library( packages_needed[i], character.only = TRUE)
}

###################################################################################

#Import demography estimates + metadata
demog_means <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv") %>% 
  dplyr::select(Site,Paper_ID,Lat.Color)

# Pops used in genome sequencing
all_pop <- read_csv("data/genomic_data/all_pops_mimulus.csv") %>% 
  filter(Paper_ID %in% demog_means$Paper_ID) %>% dplyr::select(Site_Name,Lat,Long,Paper_ID)

#Merge
df <- demog_means <- left_join(demog_means,all_pop,by="Paper_ID")

###################################################################################
### Map of focal populations

# get world map
#world <- ne_countries(scale = "medium", returnclass = "sf")
states_sf <- st_as_sf(map("state", regions = c("california", "oregon"), plot = FALSE, fill = TRUE))


# create map of states
states <- map_data("state")

# make bounding box
ymin = min(all_pop$Lat) - 1.2
ymax = max(all_pop$Lat) + 1.2
xmin = min(all_pop$Long) - 1.2 
xmax = max(all_pop$Long) + 1.2 

bbox <- st_bbox(states_sf)

xmin <- bbox["xmin"] #- 1
xmax <- bbox["xmax"] + 1
ymin <- bbox["ymin"] - 1
ymax <- bbox["ymax"] #+ 1

# Map with numbers 
card_map_focal_pops = ggplot(data=states_sf,fill="white",col="black",size=0.3) + 
  geom_sf() +
  geom_polygon(aes(x = long, y = lat, group = group), data=states,fill="white",col="black",size=0.2) +
  theme_minimal() +
  coord_sf(xlim = c(xmin,xmax), ylim = c(ymin+0.5,ymax), expand = FALSE) +
  #geom_point(aes(x=Long,y=Lat),fill = df$Lat.Color, data=df,shape=21,size=4) +
  #geom_text(aes(x = Long, y = Lat, label = Paper_ID, color = Lat.Color), data = df, size = 4,fontface = "bold")+
  geom_text_repel(
    aes(x = Long, y = Lat, label = Paper_ID, color = Lat.Color),
    data = df,
    size = 5,
    force = 0.001,          # default is 1, try lowering to 0.5 or 0.3
    max.overlaps = Inf,  # prevent dropping any labels
    fontface = "bold"
    ) +
  scale_color_identity()+
  #scale_fill_identity()+
  labs(x="Longitude",y="Latitude") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=16),
        legend.position = "none",
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        plot.title=element_text(hjust=0,size=18),
        legend.title = element_blank()) 

card_map_focal_pops


#Use Labels instead
map_demo = ggplot(data=states_sf,fill="white",col="black",size=0.3) + 
  geom_sf() +
  geom_polygon(aes(x = long, y = lat, group = group), data=states,fill="white",col="black",size=0.2) +
  theme_minimal() +
  coord_sf(xlim = c(xmin,xmax), ylim = c(ymin+0.5,ymax), expand = FALSE) +
  
  geom_point(aes(x = Long, y = Lat, fill = Lat.Color),
             data = df,
             shape = 21,
             color = "black",
             size = 4,
             stroke = 1)+#,
             #position = position_jitter(width = 0.2, height = 0.2))+
  
  # Labels next to points, repel to avoid overlap
  geom_text_repel(aes(x = Long, y = Lat, label = Paper_ID),
                  data = df,
                  size = 4,
                  fontface = "bold",
                  color = "black",
                  nudge_y = 0.1,   # small vertical nudge so labels don’t cover points
                  segment.color = "grey50") +
  scale_fill_identity()+
  labs(x="Longitude",y="Latitude") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=16),
        legend.position = "none",
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        plot.title=element_text(hjust=0,size=18),
        legend.title = element_blank()) 

map_demo


# save genomics sampling map 
save_plot(filename = "Graphs/Maps/demo_map.png",
          plot = map_demo,
          bg = "transparent", 
          base_width = 5, base_height = 8)
save_plot(filename = "Graphs/Maps/demo_map.pdf",
          plot = map_demo,
          bg = "transparent", 
          base_width = 5, base_height = 8)

