##################################################################################
#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Test whether lambda is predicted by genetic diversity (pi) 
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
library(sfsmisc)

#Import demography data (includes metadata)
demog_recovery <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv")

#Import pi for demography populations
pi_raw <- read_csv("data/genomic_data/raw_pi_clump.csv") 

#Join data sets
pi_pop <- left_join(demog_recovery, pi_raw, by=c("Paper_ID"="Site")) 

###################################################################################
###### PI OF CLIMATE SNP

#Visualize scatter plot
ggplot(pi_pop, aes(x=pi_snp_set, y=mean.r.recovery)) + geom_point()
pi_pop_graph <- pi_pop %>% dplyr::select(Site,Paper_ID,Latitude,Lat.Color,mean.r.recovery,pi_snp_set,pi_all_snps)
pi_pop_graph <- drop_na(pi_pop_graph)

#Ordinary least squares regression
mod_pi_snp <- lm(mean.r.recovery~pi_snp_set,data=pi_pop)
summary(mod_pi_snp) #b=33.25, p=0.00475
Anova(mod_pi_snp,type="III")
qqnorm(resid(mod_pi_snp))
qqline(resid(mod_pi_snp))

#Robust regression instead
rob.mod_pi_snp <- rlm(pi_pop$mean.r.recovery~pi_snp_set,data=pi_pop)
summary(rob.mod_pi_snp)
f.robftest(rob.mod_pi_snp, var="pi_snp_set") #p=0.003963


###### GENOME-WIDE PI 

#Visualize scatter plot
ggplot(pi_pop, aes(x=pi_all_snps, y=mean.r.recovery)) + geom_point()

#Ordinary least squares regression
mod_pi_all <- lm(mean.r.recovery~pi_all_snps,data=pi_pop)
summary(mod_pi_all) #b=-8.337, p=0.554
Anova(mod_pi_all,type="III")
qqnorm(resid(mod_pi_all))
qqline(resid(mod_pi_all))

#Robust regression instead
rob.mod_pi_all <- rlm(pi_pop$mean.r.recovery~pi_all_snps,data=pi_pop)
summary(rob.mod_pi_all)
f.robftest(rob.mod_pi_all, var="pi_all_snps") #p=0.8606



#####################################################################################

### Graphs of lambda recovery vs genetic diversity

## CLIMATE SNP (FIG 3D)
pi_pop_graph$Lat.Color<-as.factor(pi_pop_graph$Lat.Color)
pi_pop_graph$Lat.Color<-factor(pi_pop_graph$Lat.Color,levels=pi_pop_graph$Lat.Color)

#Depicting ordinary regression with and without robust regression
ggplot(pi_pop_graph, aes(x=pi_snp_set, y=mean.r.recovery)) +
  geom_smooth(method=lm,color="black", size=1.8, linetype="dashed", fill="gray75")+
  geom_smooth(method=MASS::rlm, color="black", size=1.8, fill="grey50") +
  geom_point(aes(fill=Lat.Color), shape=21, size =6)+
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.7) +
  scale_y_continuous(name="Mean Pop. Growth after Drought")+
  #, limits=c(-0.3,2.5),breaks=seq(0,2.5,0.5))+
  scale_x_continuous(name="Pi (Climate SNP)")+
  #, limits=c(0.2,0.35), breaks=seq(0.1,0.35,0.05)) +  
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

ggsave("Graphs/Demography_2/06_pi_demography_snpset_all.pdf",width=8, height = 6, units = "in")


## GENOME SNP (FIG 3E)

#Depicting ordinary regression with and without robust regression
ggplot(pi_pop_graph, aes(x=pi_all_snps, y=mean.r.recovery)) +
  geom_smooth(method=lm,color="black", size=1.8, linetype="dashed", fill="gray75")+
  geom_smooth(method=MASS::rlm, color="black", size=1.8, fill="grey50") +
  geom_point(aes(fill=Lat.Color), shape=21, size =6)+
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.7) +
  scale_y_continuous(name="Mean Pop. Growth after Drought")+
  #, limits=c(-0.3,2.5),breaks=seq(0,2.5,0.5))+
  scale_x_continuous(name="Pi (All SNP)")+
  #, limits=c(0.2,0.35), breaks=seq(0.1,0.35,0.05)) +  
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
ggsave("Graphs/Demography_2/07_pi_demography_global_all.pdf",width=8, height = 6, units = "in")

