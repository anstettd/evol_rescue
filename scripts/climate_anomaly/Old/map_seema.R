##################################################################################
#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Generate map of baseline and time series populations
#### AUTHOR: Modified from Seema Sheth's script from 2018 PNAS paper
#### DATE LAST MODIFIED: 4 Sep 2024 BUT THIS CODE IS NOT WORKING YET - DO NOT USE
###################################################################################

# remove objects and clear workspace
rm(list = ls(all=TRUE))

#load packages
library(sp) #no longer actively developed --> best to switch to sf
#library(rgdal) #not available for this version of R
#library(maptools) #not available for this version of R
#library(rgeos) #not available for this version of R
library(sf)
library(maps)
library(mapdata)
library(raster)
library(rworldmap)
library(rworldxtra)
library(RColorBrewer)
#library(plyr)
#library(dplyr)
library(tidyverse)


#*************************************
# Read in population locations
#*************************************

#read in M. cardinalis demography data and extract lat/lon for each population
#data=read.csv("Data/Mcard_demog_data_2010-2013.csv") %>% select(Latitude,Longitude) %>% unique() %>% arrange(-Latitude)

#Baseline & Timeseries
all_pop <- read_csv("data/genomic_data/Baseline_Timeseries_pops_final2.csv")%>% filter(Paper_ID<56)
#all_pop_sf <- st_as_sf(all_pop,coords=c("Long","Lat"), crs=EPSG4326)
timeseries_pop <- all_pop %>% filter(Paper_ID<12) %>% dplyr::select(Long,Lat)
baseline_pop <- all_pop %>% filter(Paper_ID>12) %>% dplyr::select(Long,Lat)

#Demography
#demo_pop <- read_csv("data/genomic_data/pop_meta_data.csv")
#demo_pop_sf <- st_as_sf(demo_pop ,coords=c("Longitude","Latitude"), crs=EPSG4326)


#read in all occurrence data from occupancy MS + Baja records and merge into one data frame
#localities=read.csv("data/raw_data/all.records.aug.31.csv")
#localities=subset(localities,DATASET=="herb"&PRESABS==1)
#localities=subset(localities,select=c("Latitude","Longitude"))
#baja=read.csv("data/raw_data/Baja.csv")
#baja=subset(baja,select=c("Lat","Lon"))
#baja=baja[1:14,] # Remove AZ and Cedros Island localities
#colnames(baja)=c("Latitude","Longitude")
#localities=rbind(localities,baja)

#create numeric codes for sites
#data$Site.code=rev(seq(1,32,1))
#lat_cols=colorRampPalette(brewer.pal(11,"Spectral"))
#data$Col=lat_cols(32)[as.numeric(cut(data$Latitude,breaks = 32))]

#define projections
prj.wgs = "+proj=longlat +ellps=WGS84"
prj.aea = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0"

#convert to spatial data
timeseries_pop_sf_wgs <- st_as_sf(timeseries_pop,coords=c("Long","Lat"), crs=prj.wgs)
baseline_pop_sf_wgs <- st_as_sf(baseline_pop,coords=c("Long","Lat"), crs=prj.wgs)

#coordinates(timeseries_pop) <- ~Long+Lat
#class(timeseries_pop)
#coordinates(baseline_pop)<- ~Long+Lat

#re-project pop locations
timeseries_pop_sf_aea = st_transform(timeseries_pop_sf_wgs, CRS=CRS(prj.aea))
baseline_pop_sf_aea = st_transform(baseline_pop_sf_wgs, CRS=CRS(prj.aea))

# crop to buffered bounding box
ymin = min(baseline_pop$Lat) - 1
ymax = max(baseline_pop$Lat) + 1
xmin = min(baseline_pop$Long) - 1 
xmax = max(baseline_pop$Long) + 1
e = c(xmin, ymin, xmax, ymax)
#e = extent(xmin, xmax, ymin, ymax)
bbox.wgs = st_bbox(e)
bbox.aea = st_transform(bbox.wgs, CRS(prj.aea)) #cannot st_transform a bbox

# ecoregions (why? for gridlines?) 
ecoreg = read_sf(dsn="data/map_data/ecoregions.shp", layer="us_eco_l3_no_st") 

# project to match sites
ecoreg.wgs = st_transform(ecoreg, CRS(prj.wgs))
ecoreg.aea = st_transform(ecoreg, CRS(prj.aea))

# crop polygon to buffered bounding box
bbox = st_cast(bbox.wgs, "POLYGON")
#as(e, "SpatialPolygons")
#proj4string(bbox) = CRS(prj.wgs)
#bbox.aea = st_transform(bbox, CRS=CRS(prj.aea))
ecoreg.crop.aea <- st_crop(ecoreg.aea, bbox.aea) 
ecoreg.crop.wgs <- st_crop(ecoreg.wgs, bbox.wgs) 

## world map polygons
data(countriesLow)
data(countriesHigh)
countriesLow.aea = st_transform(countriesLow, CRS=CRS(prj.aea))
countriesHigh.aea = st_transform(countriesHigh, CRS=CRS(prj.aea))

# North American polygon
northAmerica.aea=subset(countriesLow.aea,continent=="North America"|GEO3=="Meso-America")

## state polygons
#sta = readOGR("data/raw_data/gz_2010_us_040_00_500k/gz_2010_us_040_00_500k.shp")
#projection(sta) = CRS(prj.wgs)
#sta.aea = spTransform(sta, CRS=CRS(prj.aea))

## gridlines
# create unprojected gridlines
grd.wgs = gridlines(ecoreg.crop.wgs, ndiscr=200)

# project gridlines
grd.aea = spTransform(grd.wgs, CRS=CRS(prj.aea))

# Make base plotting frame with gridlines
frame.aea <- crop(countriesHigh.aea, bbox.aea); 
frame.wgs = spTransform(frame.aea, CRS=CRS(prj.wgs))
frame.grd = gridlines(frame.wgs, ndiscr=100)
frame.grd.aea = spTransform(frame.grd, CRS=CRS(prj.aea))

# prepare labels for gridlines in unprojected space
gridat <- gridat(frame.grd, side="EN")
#project labels for gridlines
gridat.aea=spTransform(gridat, CRS=CRS(prj.aea), side="EN")

# slim down labels to fit on plots
lab.all = parse(text=as.character(gridat.aea$labels))
lab.cull = lab.all[c(2,4,6,7,8,9)]
coord.all = coordinates(gridat.aea)
coord.cull = coord.all[c(2,4,6,7,8,9),]

# plot North America map that will be inserted as inset in Fig. 2a
# print to pdf
pdf(file="Figures/Fig2a_map_inset.pdf",width=11,height=8.5)
plot(northAmerica.aea,col="white")
plot(frame.aea,add=TRUE,col="lightgrey")
dev.off()

# print to pdf; Fig. 2a
pdf(file="Figures/Fig2a.pdf",width=11,height=8.5)

par(mfrow=c(1,3),las=1,bty="l",xpd=NA,cex.lab=2.3,mar=c(10,5,4,4.5)+.1)
plot(frame.aea, col="lightgrey",xlim=c(-2004864,-1852484))
sta.aea2 <- crop(sta.aea,frame.aea); plot(sta.aea2, col="lightgrey",add=TRUE)
grd.aea2 <- crop(frame.grd.aea,frame.aea); plot(grd.aea2, add=T, lty="dashed", col="darkgrey", lwd=1)
points(localities.aea,pch=21,col="black",bg="white",cex=2)
plot(sta.aea2,add=T,border="black",lwd=0.2)
points(dat.aea, pch=21, col="black", bg=adjustcolor(data$Col,alpha=0.75),cex=4) # this weights color by latitude
text(coord.cull, labels=lab.cull, offset=0.5, pos=c(3,3,4,4,4,4),col="black",cex=2) #this works for single panels, but not multi-panel
mtext("A",side=3,cex=1.5,adj=-0.18)
dev.off()















