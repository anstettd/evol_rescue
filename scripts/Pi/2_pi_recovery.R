##################################################################################
#### PROJECT: Evolutionary rescue of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Test whether lambda is predicted by genetic diversity (pi) 
#### AUTHOR: Daniel Anstett and Amy Angert
#### DATE LAST MODIFIED: 20240628
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
pi_raw <- read_csv("data/genomic_data/raw_pi.csv")

#Join data sets
pi_pop <- left_join(demog_recovery, pi_raw, by=c("Paper_ID"="Site")) 

#Visualize scatter plot
ggplot(pi_pop, aes(x=pi_snp_set, y=mean.lambda.recovery)) + geom_point()
# looks like some extreme outliers - needs test

#Cook's distance check for influential outliers
mod.check <- lm(mean.lambda.recovery~pi_snp_set,data=pi_pop)
cooksd <- cooks.distance(mod.check)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels

#Cull outlier population
pi_pop_cull1 <- pi_pop %>% 
  filter(Paper_ID!=12) #remove influential outlier Mill Creek

#Cook's distance check for influential outliers after removing Mill
mod.check_cull1 <- lm(mean.lambda.recovery~pi_snp_set,data=pi_pop_cull1)
cooksd_cull1 <- cooks.distance(mod.check_cull1)
plot(cooksd_cull1, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd_cull1, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd_cull1)+1, y=cooksd_cull1, labels=ifelse(cooksd_cull1>4*mean(cooksd_cull1, na.rm=T),names(cooksd_cull1),""), col="red")  # add labels

#Cull another outlier population
pi_pop_cull2 <- pi_pop_cull1 %>% 
  filter(Paper_ID!=27) # remove influential outlier Buck Meadows

#Cook's distance check for influential outliers after removing Mill + Buck
mod.check_cull2 <- lm(mean.lambda.recovery~pi_snp_set,data=pi_pop_cull2)
cooksd_cull2 <- cooks.distance(mod.check_cull2)
plot(cooksd_cull2, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd_cull2, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd_cull2)+1, y=cooksd_cull2, labels=ifelse(cooksd_cull2>4*mean(cooksd_cull2, na.rm=T),names(cooksd_cull2),""), col="red")  # add labels

#Cull another outlier population
pi_pop_cull3 <- pi_pop_cull2 %>% 
  filter(Paper_ID!=4) # remove influential outlier Redwood

#Cook's distance check for influential outliers after removing Mill + Buck + Redwood
mod.check_cull3 <- lm(mean.lambda.recovery~pi_snp_set,data=pi_pop_cull3)
cooksd_cull3 <- cooks.distance(mod.check_cull3)
plot(cooksd_cull3, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd_cull3, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd_cull3)+1, y=cooksd_cull3, labels=ifelse(cooksd_cull3>4*mean(cooksd_cull3, na.rm=T),names(cooksd_cull3),""), col="red")  # add labels

#Cull another outlier population
pi_pop_cull4 <- pi_pop_cull3 %>% 
  filter(Paper_ID!=1) # remove influential outlier Sweetwater R

#Cook's distance check for influential outliers after removing Mill + Buck + Redwood + Sweetwater
mod.check_cull4 <- lm(mean.lambda.recovery~pi_snp_set,data=pi_pop_cull4)
cooksd_cull4 <- cooks.distance(mod.check_cull4)
plot(cooksd_cull4, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd_cull4, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd_cull4)+1, y=cooksd_cull4, labels=ifelse(cooksd_cull4>4*mean(cooksd_cull4, na.rm=T),names(cooksd_cull4),""), col="red")  # add labels
#no more influential outliers

#Run linear models
lm_snp <- lm(mean.lambda.recovery~pi_snp_set,data=pi_pop_cull4)
lm_all <- lm(mean.lambda.recovery~pi_all_snps,data=pi_pop_cull4)

#Get summary and Anova for each model
summary(lm_snp)
Anova(lm_snp,type="III")
summary(lm_all)
Anova(lm_all,type="III")


###########################################################################################################

# Graphs of lambda recovery vs genetic diversity

pi_pop_graph <- drop_na(pi_pop_cull4)

# pi snp set
ggplot(pi_pop_graph, aes(x=pi_snp_set, y=mean.lambda.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))), shape=21, size =6)+
  geom_smooth(method=lm,color="black")+
  scale_y_continuous(name="Mean Lambda after Drought", limits=c(-0.3,2.5), breaks=seq(0,2.5,0.5))+
  scale_x_continuous(name="Pi (Climate SNP)", limits=c(0.2,0.35), breaks=seq(0.1,0.35,0.05)) +  
  scale_fill_manual(values=pi_pop_graph$Lat.Color) +
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
ggsave("Graphs/Demography_pi/3_pi_demography_snpset.pdf",width=8, height = 6, units = "in")


#global pi
ggplot(pi_pop_graph, aes(x=pi_all_snps, y=mean.lambda.recovery)) + 
  geom_point(aes(fill=as.factor(round(Latitude, 1))),shape=21,size =6)+
  geom_smooth(method=lm,color="black", lty="dashed", se=FALSE)+
  scale_y_continuous(name="Mean Lambda after Drought", limits=c(-0.3,2.5), breaks=seq(0,2.5,0.5))+
  scale_x_continuous(name="Pi (Genome-Wide)")+#, breaks=c(0.04,0.045,0.05,0.055,0.06))+
  scale_fill_manual(values=pi_pop_graph$Lat.Color) +
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

ggsave("Graphs/Demography_pi/4_pi_demography_global.pdf",width=8, height = 6, units = "in")






