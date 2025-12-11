##################################################################################
#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Correlate lambda decline vs recovery as a quick check for density dependence
#### AUTHOR: Daniel Anstett
#### DATE LAST MODIFIED: 20251211
###################################################################################

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

#Library install and import
library(tidyverse) 
library(car)
library(RColorBrewer)
library(cowplot)

#Import population metadata
#pop_meta <- read_csv("data/genomic_data/pop_meta_data.csv") %>% dplyr::select(Site, Paper_ID) #just to get translation of pop names <--> numbers
#pop_meta[20,1] <- "Mill Creek" #Fill in missing information
#pop_meta[20,2] <- 12

#Import Demography Data
demog_recovery <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv") %>% 
  filter(Paper_ID!=12) 


###########################################################################################################

hist(demog_recovery$geomean.l.drought)
hist(demog_recovery$mean.r.drought)
hist(demog_recovery$geomean.l.recovery)
hist(demog_recovery$mean.r.recovery)

#cor.test(demog_recovery$geomean.l.drought, demog_recovery$geomean.l.recovery, method = 'spearman')

cor.test(demog_recovery$mean.r.drought, demog_recovery$mean.r.recovery, method = 'spearman')

demog_recovery_cull <- demog_recovery %>% filter(mean.r.drought>-2)
cor.test(demog_recovery_cull$mean.r.drought, demog_recovery_cull$mean.r.recovery, method = 'spearman')
# still not correlated after removing three populations that were extirpated

demog_recovery$Lat.Color<-as.factor(demog_recovery$Lat.Color)
demog_recovery$Lat.Color<-factor(demog_recovery$Lat.Color,levels=demog_recovery$Lat.Color)

#Plot lambda decline vs lambda recovery (Fig. S9)
ggplot(demog_recovery, aes(x=mean.r.drought, y=mean.r.recovery)) + 
  geom_point(aes(fill=Lat.Color), shape=21, size =6)+
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.7) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = 0.7) +
  scale_y_continuous(name="Mean Pop. Growth after Drought")+
  scale_x_continuous(name="Mean Pop. Growth during Drought")+
  #scale_fill_manual(name = "Latitude (°N)",labels=round(demog_recovery$Latitude,1), values=as.character(demog_recovery$Lat.Color)) +
  scale_fill_manual(values=as.character(demog_recovery$Lat.Color), 
                    labels = unique(demog_recovery$Paper_ID) ) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    #legend.title = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce height
  ) +
  guides(
    color = guide_legend(reverse = TRUE, override.aes = list(linetype = 0)),
    fill  = guide_legend(reverse = TRUE)
  )

ggsave("Graphs/Demography/decline_vs_recovery.pdf",width=8, height = 8, units = "in")
 
