##################################################################################
## Calculate median S
## Plot along chromosomes
## 
## Author Daniel Anstett
## 
## 
## Last Modified 07082025
###################################################################################
# Remove objects and clear workspace
rm(list = ls(all=TRUE))

#Function
theme_ci <- function(){ 
  theme_classic() %+replace%    #replace elements we want to change
    theme(axis.text.x = element_text(size = 14, face = "bold", angle = 0,hjust = 0.4, vjust = 0.7), 
          axis.title = element_text(size = 18, face = "bold"), 
          axis.text.y = element_text(size = 16, face = "bold"),
          strip.background = element_blank(),strip.text.x = element_text(size = 16, face = "bold"))
}
###################################################################################
#Import libraries
library(tidyverse)
library(car)

#Import
snp_all <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_all.csv") %>% filter(SE<5) %>% filter(Site!=12)

#Get top 212 SNPs per pop
p1_212 <- snp_all %>% filter(Site==1) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 212)
p2_212 <- snp_all %>% filter(Site == 2) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 212)
p3_212 <- snp_all %>% filter(Site == 3) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 212)
p4_212 <- snp_all %>% filter(Site == 4) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 212)
p5_212 <- snp_all %>% filter(Site == 5) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 212)
p6_212 <- snp_all %>% filter(Site == 6) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 212)
p7_212 <- snp_all %>% filter(Site == 7) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 212)
p8_212 <- snp_all %>% filter(Site == 8) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 212)
p9_212 <- snp_all %>% filter(Site == 9) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 212)
p10_212 <- snp_all %>% filter(Site == 10) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 212)
p11_212 <- snp_all %>% filter(Site == 11) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 212)


all_212 <- rbind(p1_212,
                 p2_212,
                 p3_212,
                 p4_212,
                 p5_212,
                 p6_212,
                 p7_212,
                 p8_212,
                 p9_212,
                 p10_212,
                 p11_212)



##Get median S for each pop
snp_all_median <- all_212 %>% group_by(Site) %>% summarise(median = median(Slope_abs, na.rm = TRUE))

snp_all_median$Site <- as.factor(snp_all_median$Site)
snp_all_median$Site <- factor(snp_all_median$Site, levels = c(1,2,3,4,5,6,7,8,9,10,11))

#Median 
histPop <- ggplot(all_212,aes(x=Slope_abs))+
  geom_histogram(color="black",fill = "grey70")+
  labs(x = "Response to Selection", y = "Count") +
  #geom_vline(xintercept=0) +
  theme_ci() + facet_wrap(.~Site) +
  geom_vline(data = snp_all_median, aes(xintercept = median), size=0.5, linetype="dashed",color="red")
histPop 

ggsave("Graphs/Regression_all/median_high_S.pdf",width=12, height = 8, units = "in")

###################################################################################
#High S Median vs demograph

#Import demography data + metadata
demo_pop <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv")

#Join frames
snp_all_median <- snp_all_median %>%
  mutate(Site = as.numeric(Site))
slope_pop <- left_join(demo_pop, snp_all_median, by=c("Paper_ID"= "Site"))
slope_pop <- slope_pop %>% dplyr::select(Latitude,Site,Paper_ID,Lat.Color,mean.r.recovery,median)
slope_pop <- drop_na(slope_pop)

#Visualize scatter plot
ggplot(slope_pop, aes(x=median, y=mean.r.recovery)) + geom_point()

mod_S <- lm(mean.r.recovery~median,data=slope_pop)
summary(mod_S)
Anova(mod_S,type="III")
qqnorm(resid(mod_S))
qqline(resid(mod_S))

#Robust regression instead
library(MASS)
rob.mod_S <- rlm(slope_pop$mean.r.recovery~median,data=slope_pop)
summary(rob.mod_S)
# this yields a coefficient estimate & std error but no p-value or R2
library(sfsmisc)
f.robftest(rob.mod_S, var="median")


###################################################################################
slope_pop$Lat.Color<-as.factor(slope_pop$Lat.Color)
slope_pop$Lat.Color<-factor(slope_pop$Lat.Color,levels=slope_pop$Lat.Color)

ggplot(slope_pop, aes(x=median, y=mean.r.recovery)) + 
  geom_smooth(method=lm,color="black", size=1.8, linetype="dashed", fill="gray75")+
  geom_smooth(method=MASS::rlm, color="black", size=1.8, fill="grey50") +
  #geom_smooth(data=slope_pop_cull1, aes(x=Median, y=log(mean.lambda.recovery+0.5)), method=lm, color="black", size=1.8, linetype="longdash", fill="grey35") +
  geom_point(aes(fill=Lat.Color), shape=21, size =6)+
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.7) +
  scale_y_continuous(name="Mean Pop. Growth after Drought")+#,breaks=c(0.5,1,1.5,2,2.5))+
  scale_x_continuous(name="Median Response to Selection")+#,breaks=c(-0.1,0,0.1,0.2))+
  #,breaks=c(0.025,0.03,0.035,0.04,0.045))+
  scale_fill_manual(name = "Latitude (°N)",labels=round(slope_pop$Latitude,1),
                    values=as.character(slope_pop$Lat.Color)) +
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
ggsave("Graphs/Regression_all/highS_demo.pdf",width=8, height = 6, units = "in")



