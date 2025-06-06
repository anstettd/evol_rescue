##################################################################################
#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Test whether lambda is predicted by genetic diversity (pi) 
#### AUTHOR: Daniel Anstett and Amy Angert
#### DATE LAST MODIFIED: 20250310
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
ggplot(pi_pop, aes(x=pi_snp_set, y=log(mean.lambda.recovery))) + geom_point()
# looks like some extreme outliers - needs test

mod <- lm(log(mean.lambda.recovery+0.1)~pi_snp_set,data=pi_pop)

#Cook's distance check for influential outliers
mod.check <- lm(mean.lambda.recovery~pi_snp_set,data=pi_pop)
cooksd_cull1 <- cooks.distance(mod.check)
plot(cooksd_cull1, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd_cull1, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd_cull1)+1, y=cooksd_cull1, labels=ifelse(cooksd_cull1>4*mean(cooksd_cull1, na.rm=T),names(cooksd_cull1),""), col="red")  # add labels

#Cull another outlier population
pi_pop_cull2 <- pi_pop %>% 
  filter(row_number()!=9) # remove influential outlier Buck Meadows (9th row of data frame, as labeled above)

#Cook's distance check for influential outliers after removing Buck
mod.check_cull2 <- lm(mean.lambda.recovery~pi_snp_set,data=pi_pop_cull2)
cooksd_cull2 <- cooks.distance(mod.check_cull2)
plot(cooksd_cull2, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd_cull2, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd_cull2)+1, y=cooksd_cull2, labels=ifelse(cooksd_cull2>4*mean(cooksd_cull2, na.rm=T),names(cooksd_cull2),""), col="red")  # add labels

#Cull another outlier population
pi_pop_cull3 <- pi_pop_cull2 %>% 
  filter(row_number()!=7) # remove influential outlier Redwood Creek (7th row of data frame, as labeled above)

#Cook's distance check for influential outliers after removing Buck + Redwood
mod.check_cull3 <- lm(mean.lambda.recovery~pi_snp_set,data=pi_pop_cull3)
cooksd_cull3 <- cooks.distance(mod.check_cull3)
plot(cooksd_cull3, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd_cull3, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd_cull3)+1, y=cooksd_cull3, labels=ifelse(cooksd_cull3>4*mean(cooksd_cull3, na.rm=T),names(cooksd_cull3),""), col="red")  # add labels

#Cull two more outlier populations
pi_pop_cull4 <- pi_pop_cull3 %>% 
  filter(row_number()!=3, row_number()!=8) # remove influential outliers Sweetwater River and Carlon

#Cook's distance check for influential outliers after removing Buck + Redwood + Sweetwater + Carlon
mod.check_cull4 <- lm(mean.lambda.recovery~pi_snp_set,data=pi_pop_cull4)
cooksd_cull4 <- cooks.distance(mod.check_cull4)
plot(cooksd_cull4, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd_cull4, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd_cull4)+1, y=cooksd_cull4, labels=ifelse(cooksd_cull4>4*mean(cooksd_cull4, na.rm=T),names(cooksd_cull4),""), col="red")  # add labels

#Cull one more outlier population
pi_pop_cull5 <- pi_pop_cull4 %>% 
  filter(row_number()!=5) # remove influential outlier NF MF Tule

#Cook's distance check for influential outliers after removing Buck + Redwood + Sweetwater + Carlon + NFTule + Whitewater
mod.check_cull5 <- lm(mean.lambda.recovery~pi_snp_set,data=pi_pop_cull5)
cooksd_cull5 <- cooks.distance(mod.check_cull5)
plot(cooksd_cull5, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd_cull5, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd_cull5)+1, y=cooksd_cull5, labels=ifelse(cooksd_cull5>4*mean(cooksd_cull5, na.rm=T),names(cooksd_cull5),""), col="red")  # add labels

#Cull one more outlier population
pi_pop_cull6 <- pi_pop_cull5 %>% 
  filter(row_number()!=3) # remove influential outlier Whitewater

#Cook's distance check for influential outliers after removing Buck + Redwood + Sweetwater + Carlon + NFTule + Whitewater
mod.check_cull6 <- lm(mean.lambda.recovery~pi_snp_set,data=pi_pop_cull6)
cooksd_cull6 <- cooks.distance(mod.check_cull6)
plot(cooksd_cull6, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd_cull6, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd_cull6)+1, y=cooksd_cull6, labels=ifelse(cooksd_cull6>4*mean(cooksd_cull6, na.rm=T),names(cooksd_cull6),""), col="red")  # add labels

#Cull one more outlier population
pi_pop_cull7 <- pi_pop_cull6 %>% 
  filter(row_number()!=5) # remove influential outlier Rainbow Pool

#Cook's distance check for influential outliers after removing Buck + Redwood + Sweetwater + Carlon + NFTule + Whitewater + Rainbow
mod.check_cull7 <- lm(mean.lambda.recovery~pi_snp_set,data=pi_pop_cull7)
cooksd_cull7 <- cooks.distance(mod.check_cull7)
plot(cooksd_cull7, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd_cull7, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd_cull7)+1, y=cooksd_cull7, labels=ifelse(cooksd_cull7>4*mean(cooksd_cull7, na.rm=T),names(cooksd_cull7),""), col="red")  # add labels

#Cull one more outlier population
pi_pop_cull8 <- pi_pop_cull7 %>% 
  filter(row_number()!=9) # remove influential outlier Rock Creek

#Cook's distance check for influential outliers after removing Buck + Redwood + Sweetwater + Carlon + NFTule + Whitewater + Rainbow + Rock
mod.check_cull8 <- lm(mean.lambda.recovery~pi_snp_set,data=pi_pop_cull8)
cooksd_cull8 <- cooks.distance(mod.check_cull8)
plot(cooksd_cull8, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd_cull8, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd_cull8)+1, y=cooksd_cull8, labels=ifelse(cooksd_cull8>4*mean(cooksd_cull8, na.rm=T),names(cooksd_cull8),""), col="red")  # add labels

#Cull one more outlier population
pi_pop_cull9 <- pi_pop_cull8 %>% 
  filter(row_number()!=2) # remove influential outlier Kitchen Creek

#Cook's distance check for influential outliers after removing Buck + Redwood + Sweetwater + Carlon + NFTule + Whitewater + Rainbow + Rock + Kitchen
mod.check_cull9 <- lm(mean.lambda.recovery~pi_snp_set,data=pi_pop_cull9)
cooksd_cull9 <- cooks.distance(mod.check_cull9)
plot(cooksd_cull9, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd_cull9, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd_cull9)+1, y=cooksd_cull9, labels=ifelse(cooksd_cull9>4*mean(cooksd_cull9, na.rm=T),names(cooksd_cull9),""), col="red")  # add labels
#no more influential outliers


###########################################################################################################
# Graphs of lambda recovery vs genetic diversity
# Not overlain
# Stats are included
###########################################################################################################

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
#ggsave("Graphs/Demography_2/pi_demography_snpset.pdf",width=8, height = 6, units = "in")


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

#ggsave("Graphs/Demography_2/pi_demography_global.pdf",width=8, height = 6, units = "in")




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
<<<<<<< HEAD
#ggsave("Graphs/Demography_2/pi_demography_snpset_cull4.pdf",width=8, height = 6, units = "in")
=======
ggsave("Graphs/Demography_2/10_pi_demography_snpset_cull9.pdf",width=8, height = 6, units = "in")
>>>>>>> 396d2e6f55d6fef1a1b63c9c1937b575c3668000




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

<<<<<<< HEAD
#ggsave("Graphs/Demography_2/pi_demography_global_cull4.pdf",width=8, height = 6, units = "in")
=======
`ggsave("Graphs/Demography_2/12_pi_demography_global_cull9.pdf",width=8, height = 6, units = "in")

>>>>>>> 396d2e6f55d6fef1a1b63c9c1937b575c3668000

###########################################################################################################
###########################################################################################################
###########################################################################################################
# Graphs with and without outliers overlain
# Use for paper and presentations
###########################################################################################################
###########################################################################################################
###########################################################################################################


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
<<<<<<< HEAD
    legend.key.height = unit(1.6, "lines") #Reduce height
  )
ggsave("Graphs/Demography_2/03_pi_demography_snpset_all.pdf",width=8, height = 6, units = "in")
=======
    legend.key.height = unit(1.6, "lines")) #Reduce height

ggsave("Graphs/Demography_2/13_pi_demography_snpset_all.pdf",width=8, height = 6, units = "in")
>>>>>>> 396d2e6f55d6fef1a1b63c9c1937b575c3668000

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
ggsave("Graphs/Demography_2/04_pi_demography_global_all.pdf",width=8, height = 6, units = "in")

