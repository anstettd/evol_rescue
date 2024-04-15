##################################################################################
## Climate associated SNP distribution across chromosomes
## Author Daniel Anstett
## 
## 
## Last Modified Marc 22, 2021
###################################################################################


##### Annual #####
# env 1 is MAT = Mean annual temperature (°C)
# env 2 is MAP = Mean annual precipitation (mm)
# env 3 is PAS = Precipitation as snow (mm) between August in previous year and July in current year
# env 4 is EXT = Extreme temperature over 30 years
# env 5 is CMD = Hargreaves climatic moisture deficit (mm)

##### Seasonal #####
# env 6 is Tave_wt = Winter mean temperature (°C)
# env 7 is Tave_sm = Summer mean temperature (°C)
# env 8 is PPT_wt = Winter precipitation (mm)
# env 9 is PPT_sm = Summer precipitation (mm)



###################################################################################
#Import libraries
library(tidyverse)

#Import final SNP set per ENV variable
snp1_filter <- read_csv("data/genomic_data/snp_set_env1.csv")
snp2_filter <- read_csv("data/genomic_data/snp_set_env2.csv")
snp3_filter <- read_csv("data/genomic_data/snp_set_env3.csv")
snp4_filter <- read_csv("data/genomic_data/snp_set_env4.csv")
snp5_filter <- read_csv("data/genomic_data/snp_set_env5.csv")
snp6_filter <- read_csv("data/genomic_data/snp_set_env6.csv")
snp7_filter <- read_csv("data/genomic_data/snp_set_env7.csv")
snp8_filter <- read_csv("data/genomic_data/snp_set_env8.csv")
snp9_filter <- read_csv("data/genomic_data/snp_set_env9.csv")

#Add env variable identifier 
chr_obs_env1 <- snp1_filter %>% select(chr_snp) %>% separate(chr_snp,c("CE","chr","Position")) %>% select(-CE)
chr_obs_env2 <- snp2_filter %>% select(chr_snp) %>% separate(chr_snp,c("CE","chr","Position")) %>% select(-CE)
chr_obs_env3 <- snp3_filter %>% select(chr_snp) %>% separate(chr_snp,c("CE","chr","Position")) %>% select(-CE)
chr_obs_env4 <- snp4_filter %>% select(chr_snp) %>% separate(chr_snp,c("CE","chr","Position")) %>% select(-CE)
chr_obs_env5 <- snp5_filter %>% select(chr_snp) %>% separate(chr_snp,c("CE","chr","Position")) %>% select(-CE)
chr_obs_env6 <- snp6_filter %>% select(chr_snp) %>% separate(chr_snp,c("CE","chr","Position")) %>% select(-CE)
chr_obs_env7 <- snp7_filter %>% select(chr_snp) %>% separate(chr_snp,c("CE","chr","Position")) %>% select(-CE)
chr_obs_env8 <- snp8_filter %>% select(chr_snp) %>% separate(chr_snp,c("CE","chr","Position")) %>% select(-CE)
chr_obs_env9 <- snp9_filter %>% select(chr_snp) %>% separate(chr_snp,c("CE","chr","Position")) %>% select(-CE)


#Add env identfier
chr_obs_env1$env <- "MAT"
chr_obs_env2$env <- "MAP"
chr_obs_env3$env <- "PAS"
chr_obs_env4$env <- "EXT"
chr_obs_env5$env <- "CMD"
chr_obs_env6$env <- "Tave_wt"
chr_obs_env7$env <- "Tave_sm"
chr_obs_env8$env <- "PPT_wt"
chr_obs_env9$env <- "PPT_sm"

#Join env variables
chr_all <- rbind(chr_obs_env1,
                 chr_obs_env2,
                 chr_obs_env3,
                 chr_obs_env4,
                 chr_obs_env5,
                 chr_obs_env6,
                 chr_obs_env7,
                 chr_obs_env8,
                 chr_obs_env9)
chr_all$Position <- as.integer(chr_all$Position)


#############################################################################################
#Plot all 
# Calculate evenly spaced breaks for the y-axis
env_breaks <- seq_along(unique(chr_all$env))
env_labels <- rev(unique(chr_all$env))

chr1_graph_1 <- ggplot(chr_all, aes(x = Position, y = env)) +
  geom_point(size = 2, shape = 4, color = "black") +
  scale_x_continuous(name = "Position", breaks = c(0, 1e+07, 2e+07, 3e+07, 4e+07, 5e+07, 6e+07)) +
  theme_classic() +
  facet_wrap(~ chr, ncol = 2) +
  labs(y = "") +  # Clear default y-axis label
  theme(
    strip.text.x = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(color = "black", size = 14, vjust = 0.5, face = "bold"),
    axis.title.y = element_blank()
  )
chr1_graph_1 + scale_y_discrete(labels = env_labels)  # Add custom y-axis labels

ggsave("Graphs/chr_map/chr_map_env.pdf",width=8, height = 7, units = "in")
