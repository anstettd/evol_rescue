##################################################################################
#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Test whether lambda is predicted by genetic diversity (pi) 
#### for pi calalcualted from select SNP set with only three enviornmental variables
#### AUTHOR: Daniel Anstett and Amy Angert
#### DATE LAST MODIFIED: 20251211
###################################################################################

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

#Library install and import
library(tidyverse) 
library(car)
library(ggrepel)
library(cowplot)
library(MASS)


#Import demography data (includes metadata)
demog_recovery <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv")

#Import pi for demography populations
pi_raw <- read_csv("data/genomic_data/raw_pi_clump.csv") 

#Join data sets
pi_pop <- left_join(demog_recovery, pi_raw, by=c("Paper_ID"="Site")) 

#Visualize scatter plot
ggplot(pi_pop, aes(x=pi_env578, y=mean.r.recovery)) + geom_point()
pi_pop_graph <- pi_pop %>% dplyr::select(Site,Paper_ID,Latitude,Lat.Color,mean.r.recovery,pi_env578)
pi_pop_graph <- drop_na(pi_pop_graph)
lm_snp_3 <- lm(mean.r.recovery~pi_env578,data=pi_pop_graph)
summary(lm_snp_3)
Anova(lm_snp_3,type="III") #P-value = 0.004, R2= 0.68

#Robust regression instead
rob.mod_pi_all <- rlm(mean.r.recovery~pi_env578,data=pi_pop_graph)
summary(rob.mod_pi_all)
# this yields a coefficient estimate & std error that is similar to the culled model
# but no p-value or R2
library(sfsmisc)
f.robftest(rob.mod_pi_all, var="pi_env578") #p-value = 0.009




###########################################################################################################

# Graphs of lambda recovery vs genetic diversity

pi_pop_graph$Lat.Color<-as.factor(pi_pop_graph$Lat.Color)
pi_pop_graph$Lat.Color<-factor(pi_pop_graph$Lat.Color,levels=pi_pop_graph$Lat.Color)

#Depicting ordinary regression & robust regression
ggplot(pi_pop_graph, aes(x=pi_env578, y=mean.r.recovery)) +
  geom_smooth(method=lm,color="black", size=1.8, linetype="dashed", fill="gray75")+
  geom_smooth(method=MASS::rlm, color="black", size=1.8, fill="grey50") +
  #geom_smooth(data=pi_pop_graph_cull1, aes(x=pi_snp_set, y=log(mean.lambda.recovery+0.5)), method=lm, color="black", size=1.8, linetype="longdash", fill="grey35") +
  geom_point(aes(fill=Lat.Color), shape=21, size =6)+
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.7) +
  scale_y_continuous(name="Mean Pop. Growth after Drought")+
  #, limits=c(-0.3,2.5),breaks=seq(0,2.5,0.5))+
  scale_x_continuous(name="Pi (CMD + Tave_sm + PPT_wt SNP)")+
  #, limits=c(0.2,0.35), breaks=seq(0.1,0.35,0.05)) +  
  #scale_fill_manual(name = "Latitude (°N)",labels=round(pi_pop_graph$Latitude,1), values=as.character(pi_pop_graph$Lat.Color)) +
  scale_fill_manual(values=as.character(pi_pop_graph$Lat.Color), 
                    labels = unique(pi_pop_graph$Paper_ID) ) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines"))+ # Reduce height
  guides(color = guide_legend(reverse = TRUE, override.aes = list(linetype = 0)),
         fill  = guide_legend(reverse = TRUE))
ggsave("Graphs/Demography_Genomics/02_pi_meanr_env578_S13B.pdf",width=8, height = 6, units = "in")


