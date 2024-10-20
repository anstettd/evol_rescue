#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Define a fixed color for each population, so that color assignment does not vary across plots depending on which populations are included 
#### AUTHOR: Amy Angert
#### DATE LAST MODIFIED: 20240627

#*******************************************************************************
#### 1. Clean workspace and load required packages
#*******************************************************************************

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

# Make vector of packages needed
packages_needed <- c("tidyverse", "RColorBrewer")

# Install packages needed (if not already installed)
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
}

# Load packages needed
for (i in 1:length(packages_needed)){
  library(packages_needed[i], character.only = TRUE)
}

#*******************************************************************************
#### 2. Read in metadata file (sorted high to low latitude)
#*******************************************************************************
pop_meta <- read_csv("data/genomic_data/pop_meta_data.csv") %>% #filter(Paper_ID!=12) 
  arrange(Latitude)

#*******************************************************************************
#### 3. Define color ramp
#*******************************************************************************

# N-S color gradient
lat_cols=colorRampPalette(brewer.pal(11,"Spectral"))
n.sites <- length(unique(pop_meta$Site))
color.list <- as.data.frame(lat_cols(n.sites))

#*******************************************************************************
#### 4. Save color assignments
#*******************************************************************************

pop_meta <- cbind(pop_meta, color.list)
names(pop_meta)[7] = "Lat.Color"

write_csv(pop_meta, "data/genomic_data/pop_meta_data.csv")
