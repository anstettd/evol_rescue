##################################################################################
#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Test whether lambda is predicted by genetic diversity (pi) 
#### AUTHOR: Daniel Anstett and Amy Angert
#### DATE LAST MODIFIED: 20250605
###################################################################################

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

#Library install and import
library(tidyverse) 
library(car)
library(ggrepel)
library(cowplot)


#Import demography data (includes metadata)
demog_recovery <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv")

#Import pi for demography populations
pi_raw <- read_csv("data/genomic_data/raw_pi_clump.csv") 

#Join data sets
pi_pop <- left_join(demog_recovery, pi_raw, by=c("Paper_ID"="Site")) 

#Visualize scatter plot
ggplot(pi_pop, aes(x=pi_snp_set, y=mean.lambda.recovery)) + geom_point()
#looks like some big outliers - needs transformation &/or test for outliers

#Distribution of response variable
hist(pi_pop$mean.lambda.recovery) #long tail
hist(log(pi_pop$mean.lambda.recovery+0.5)) #much better

#Visualize scatter plot with log-transformed response
ggplot(pi_pop, aes(x=pi_snp_set, y=log(mean.lambda.recovery+0.5))) + geom_point()

mod <- lm(log(mean.lambda.recovery+0.5)~pi_snp_set,data=pi_pop)
summary(mod)
qqnorm(resid(mod))
qqline(resid(mod))

#Check for influential outliers
inflm <- influence.measures(mod)
summary(inflm) #observation 9 (Buck Meadows) is potentially an outlier according to cov.r

#Cull the outlier population
pi_pop_cull1 <- pi_pop %>% 
  filter(row_number()!=row.names(summary(inflm)))  # remove the 9th row of data frame

# sensitivity test without influential outlier Buck Meadows
mod.cull <- lm(log(mean.lambda.recovery+0.5)~pi_snp_set, data=pi_pop_cull1)
summary(mod.cull)
qqnorm(resid(mod.cull))
qqline(resid(mod.cull))

#Cook's distance check for influential outliers after removing Buck
inflm_cull1 <- influence.measures(mod.cull)
summary(inflm_cull1) #2 more observations flagged by cov.r, 1 more observation flagged by df metrics
#Cull more outlier populations...?? This seems highly questionable.

#Robust regression instead
library(MASS)
rob.mod <- rlm(log(pi_pop$mean.lambda.recovery+0.5)~pi_snp_set,data=pi_pop)
summary(rob.mod)
# this yields a coefficient estimate & std error that is similar to the culled model
# but no p-value or R2
library(sfsmisc)
f.robftest(rob.mod, var="pi_snp_set")

###########################################################################################################

# Graphs of lambda recovery vs genetic diversity

#### All Data #### 

pi_pop_graph <- drop_na(pi_pop)
lm_snp_3 <- lm(mean.lambda.recovery~pi_snp_set,data=pi_pop_graph)
summary(lm_snp_3)
Anova(lm_snp_3,type="III")

pi_pop_graph$Lat.Color<-as.factor(pi_pop_graph$Lat.Color)
pi_pop_graph$Lat.Color<-factor(pi_pop_graph$Lat.Color,levels=pi_pop_graph$Lat.Color)

# pi snp set
ggplot(pi_pop_graph, aes(x=pi_snp_set, y=mean.lambda.recovery)) +
  geom_smooth(method=lm,color="black",size=1.8,fill="gray75")+
  geom_point(aes(fill=pi_pop_graph$Lat.Color), shape=21, size =6)+
  geom_hline(yintercept = 1, linetype = "dotted", color = "black", size = 0.7) +
  scale_y_continuous(name="Mean Lambda after Drought")+
  #, limits=c(-0.3,2.5),breaks=seq(0,2.5,0.5))+
  scale_x_continuous(name="Pi (Climate SNP)")+
  #, limits=c(0.2,0.35), breaks=seq(0.1,0.35,0.05)) +  
  scale_fill_manual(name = "Latitude (°N)",labels=round(pi_pop_graph$Latitude,1), values=as.character(pi_pop_graph$Lat.Color)) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce height
  )
ggsave("Graphs/Demography_2/07_pi_demography_snpset.pdf",width=8, height = 6, units = "in")


#global pi
lm_all_3 <- lm(mean.lambda.recovery~pi_all_snps,data=pi_pop_graph)
summary(lm_all_3)
Anova(lm_all_3,type="III")

pi_pop_graph$Lat.Color<-as.factor(pi_pop_graph$Lat.Color)
pi_pop_graph$Lat.Color<-factor(pi_pop_graph$Lat.Color,levels=pi_pop_graph$Lat.Color)

ggplot(pi_pop_graph, aes(x=pi_all_snps, y=mean.lambda.recovery)) +
  geom_smooth(method=lm,color="black", lty="dashed", size=1.8, se=FALSE)+
  geom_point(aes(fill=pi_pop_graph$Lat.Color), shape=21, size =6)+
  geom_hline(yintercept = 1, linetype = "dotted", color = "black", size = 0.7) +
  scale_y_continuous(name="Mean Lambda after Drought")+
  #, limits=c(-0.3,2.5), breaks=seq(0,2.5,0.5))+
  scale_x_continuous(name="Pi (Genome-Wide)")+#, breaks=c(0.04,0.045,0.05,0.055,0.06))+
  scale_fill_manual(name = "Latitude (°N)",labels=round(pi_pop_graph$Latitude,1), values=as.character(pi_pop_graph$Lat.Color)) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce hight
  )

ggsave("Graphs/Demography_2/09_pi_demography_global.pdf",width=8, height = 6, units = "in")




#### All outlier removal #### 

pi_pop_graph_cull9 <- drop_na(pi_pop_cull9)
lm_snp_1 <- lm(mean.lambda.recovery~pi_snp_set,data=pi_pop_cull9)
summary(lm_snp_1)
Anova(lm_snp_1,type="III")

pi_pop_graph_cull9$Lat.Color<-as.factor(pi_pop_graph_cull9$Lat.Color)
pi_pop_graph_cull9$Lat.Color<-factor(pi_pop_graph_cull9$Lat.Color,levels=pi_pop_graph_cull9$Lat.Color)

# pi snp set
ggplot(pi_pop_graph_cull9, aes(x=pi_snp_set, y=mean.lambda.recovery)) +
  geom_smooth(method=lm,color="black",size=1.8,fill="gray75")+
  geom_point(aes(fill=pi_pop_graph_cull9$Lat.Color), shape=21, size =6)+
  scale_y_continuous(name="Mean Lambda after Drought", limits=c(-0.3,2.5), breaks=seq(0,2.5,0.5))+
  scale_x_continuous(name="Pi (Climate SNP)", limits=c(0.2,0.35), breaks=seq(0.1,0.35,0.05)) +  
  scale_fill_manual(labels=round(pi_pop_graph_cull9$Latitude,1), 
                    values=as.character(pi_pop_graph_cull9$Lat.Color)) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce height
  )
ggsave("Graphs/Demography_2/10_pi_demography_snpset_cull9.pdf",width=8, height = 6, units = "in")




#global pi
lm_all_1 <- lm(mean.lambda.recovery~pi_all_snps,data=pi_pop_graph_cull9)
summary(lm_all_1)
Anova(lm_all_1,type="III")

pi_pop_graph_cull9$Lat.Color<-as.factor(pi_pop_graph_cull9$Lat.Color)
pi_pop_graph_cull9$Lat.Color<-factor(pi_pop_graph_cull9$Lat.Color,levels=pi_pop_graph_cull9$Lat.Color)

ggplot(pi_pop_graph_cull9, aes(x=pi_all_snps, y=mean.lambda.recovery)) +
  geom_smooth(method=lm,color="black", lty="dashed", size=1.8, se=FALSE)+
  geom_point(aes(fill=pi_pop_graph_cull9$Lat.Color), shape=21, size =6)+
  scale_y_continuous(name="Mean Lambda after Drought", limits=c(-0.3,2.5), breaks=seq(0,2.5,0.5))+
  scale_x_continuous(name="Pi (Genome-Wide)")+#, breaks=c(0.04,0.045,0.05,0.055,0.06))+
  scale_fill_manual(labels=round(pi_pop_graph_cull9$Latitude,1), 
                    values=as.character(pi_pop_graph_cull9$Lat.Color)) +
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

`ggsave("Graphs/Demography_2/12_pi_demography_global_cull9.pdf",width=8, height = 6, units = "in")



# graphs with slopes with and without outliers overlain

#climate pi, with and without outliers
ggplot(pi_pop_graph, aes(x=pi_snp_set, y=mean.lambda.recovery)) +
  geom_smooth(method=lm,color="black",size=1.8,fill="gray75")+
  geom_smooth(data=pi_pop_graph_cull9, aes(x=pi_snp_set, y=mean.lambda.recovery), method=lm, color="black", linetype="dashed", fill="grey50") +
  geom_point(aes(fill=pi_pop_graph$Lat.Color), shape=21, size =6)+
  geom_hline(yintercept = 1, linetype = "dotted", color = "black", size = 0.7) +
  scale_y_continuous(name="Mean Lambda after Drought")+
  #, limits=c(-0.3,2.5),breaks=seq(0,2.5,0.5))+
  scale_x_continuous(name="Pi (Climate SNP)")+
  #, limits=c(0.2,0.35), breaks=seq(0.1,0.35,0.05)) +  
  scale_fill_manual(name = "Latitude (°N)",labels=round(pi_pop_graph$Latitude,1), values=as.character(pi_pop_graph$Lat.Color)) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines")) #Reduce height

ggsave("Graphs/Demography_2/13_pi_demography_snpset_all.pdf",width=8, height = 6, units = "in")

#global pi, with and without outliers
ggplot(pi_pop_graph, aes(x=pi_all_snps, y=mean.lambda.recovery)) +
  geom_smooth(method=lm,color="black",size=1.8,fill="gray75")+
  geom_smooth(data=pi_pop_graph_cull4, aes(x=pi_all_snps, y=mean.lambda.recovery), method=lm, color="black", linetype="dashed", fill="grey50") +
  geom_point(aes(fill=pi_pop_graph$Lat.Color), shape=21, size =6)+
  geom_hline(yintercept = 1, linetype = "dotted", color = "black", size = 0.7) +
  scale_y_continuous(name="Mean Lambda after Drought")+
  #, limits=c(-0.3,2.5),breaks=seq(0,2.5,0.5))+
  scale_x_continuous(name="Pi (All SNPs)")+
  #, limits=c(0.2,0.35), breaks=seq(0.1,0.35,0.05)) +  
  scale_fill_manual(name = "Latitude (°N)",labels=round(pi_pop_graph$Latitude,1), values=as.character(pi_pop_graph$Lat.Color)) +
  theme_classic() + theme(
    axis.text.x = element_text(size=20, face="bold"),
    axis.text.y = element_text(size=20,face="bold"),
    axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=24,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 14),  # Increase the size of the legend text
    legend.key.size = unit(2, "lines"),  # Increase the size of the legend dots
    legend.key.height = unit(1.6, "lines") #Reduce height
  )
ggsave("Graphs/Demography_2/15_pi_demography_global_all.pdf",width=8, height = 6, units = "in")

