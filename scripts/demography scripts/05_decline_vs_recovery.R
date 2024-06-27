##################################################################################
#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Correlate lambda decline vs recovery as a quick check for density dependence
#### AUTHOR: Daniel Anstett
#### DATE LAST MODIFIED: 20230815
###################################################################################

#Library install and import
library(tidyverse) 
library(car)
library(RColorBrewer)
library(cowplot)

#Import population metadata
pop_meta <- read_csv("data/genomic_data/pop_meta_data.csv") %>% dplyr::select(Site, Paper_ID) #just to get translation of pop names <--> numbers
pop_meta[20,1] <- "Mill Creek" #Fill in missing information
pop_meta[20,2] <- 12

#Import Demography Data
demog_recovery <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv") #%>% 
  #filter(Site!="Mill Creek") %>% filter(Site!= "Buck Meadows")


###########################################################################################################

# N-S color gradient
lat_cols=colorRampPalette(brewer.pal(11,"Spectral"))
n.sites <- 20
color.list <- lat_cols(n.sites)


cor.test(demog_recovery$lambda.mean.drought, demog_recovery$lambda.mean.recovery, method = 'pearson')

#Plot lambda decline vs lamda recovery
ggplot(demog_recovery, aes(x=lambda.mean.drought, y=lambda.mean.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =6)+
  geom_smooth(method=lm,color="black", lty="dashed", se=FALSE)+
  scale_y_continuous(name="Mean Lambda after Drought")+#, limits=c(-0.3,2.5),breaks=seq(0,2.5,0.5))+
  scale_x_continuous(name="Mean Lambda during Drought")+
  #, limits=c(0.2,0.35),breaks=seq(0.1,0.35,0.05))+
  scale_fill_manual(values=color.list) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce hight
  )

ggsave("Graphs/Demography/decline_vs_recovery.pdf",width=8, height = 6, units = "in")

