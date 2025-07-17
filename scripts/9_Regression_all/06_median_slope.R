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

#Import top plink clumped 215 SNPs per pop
p1_215 <- read_csv("data/genomic_data/p01_215.csv")
p2_215 <- read_csv("data/genomic_data/p02_215.csv")
p3_215 <- read_csv("data/genomic_data/p03_215.csv")
p4_215 <- read_csv("data/genomic_data/p04_215.csv")
p5_215 <- read_csv("data/genomic_data/p05_215.csv")
p6_215 <- read_csv("data/genomic_data/p06_215.csv")
p7_215 <- read_csv("data/genomic_data/p07_215.csv")
p8_215 <- read_csv("data/genomic_data/p08_215.csv")
p9_215 <- read_csv("data/genomic_data/p09_215.csv")
p10_215 <- read_csv("data/genomic_data/p10_215.csv")
p11_215 <- read_csv("data/genomic_data/p11_215.csv")


all_215 <- rbind(p1_215,
                 p2_215,
                 p3_215,
                 p4_215,
                 p5_215,
                 p6_215,
                 p7_215,
                 p8_215,
                 p9_215,
                 p10_215,
                 p11_215)



##Get median S for each pop
snp_all_median <- all_215 %>% group_by(Site) %>% summarise(median = median(Slope_abs, na.rm = TRUE))

snp_all_median$Site <- as.factor(snp_all_median$Site)
snp_all_median$Site <- factor(snp_all_median$Site, levels = c(1,2,3,4,5,6,7,8,9,10,11))

#Import demography data + metadata
demo_pop <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv")

#Join frames
snp_all_median <- snp_all_median %>%
  mutate(Site = as.numeric(Site))
slope_pop <- left_join(demo_pop, snp_all_median, by=c("Paper_ID"= "Site"))
slope_pop <- slope_pop %>% dplyr::select(Latitude,Site,Paper_ID,Lat.Color,mean.r.recovery,median)
slope_pop <- drop_na(slope_pop)

slope_pop$Lat.Color<-as.factor(slope_pop$Lat.Color)
slope_pop$Lat.Color<-factor(slope_pop$Lat.Color,levels=slope_pop$Lat.Color)

###################################################################################
#Median and S histogram
histPop <- ggplot(all_215,aes(x=Slope_abs))+
  geom_histogram(color="black",fill = "grey70")+
  labs(x = "Response to Selection", y = "Count") +
  #geom_vline(xintercept=0) +
  theme_ci() + facet_wrap(.~Site) +
  geom_vline(data = snp_all_median, aes(xintercept = median), size=0.5, linetype="dashed",color="red")
histPop 

ggsave("Graphs/Regression_all/median_high_S.pdf",width=12, height = 8, units = "in")


#Medians
ggplot(slope_pop, aes(x = Paper_ID, y = median)) +
  geom_point(aes(fill=Lat.Color), shape=21, size =3)+
  coord_flip() +
  scale_y_continuous(expand = c(0, 0), limits = c(1.5, 3.5)) +
  scale_fill_manual(name = "Latitude (°N)",labels=round(slope_pop$Latitude,1),
                    values=as.character(slope_pop$Lat.Color)) +
  theme_classic() +
  labs(x = "Site", y = "Median value", title = "Median per Site")


###################################################################################
#High S Median vs demograph
#Visualize scatter plot
ggplot(slope_pop, aes(x=median, y=mean.r.recovery)) + geom_point()

mod_S <- lm(mean.r.recovery~median,data=slope_pop)
summary(mod_S) #P= 0.88, R2= -0.12
Anova(mod_S,type="III")
qqnorm(resid(mod_S))
qqline(resid(mod_S))

#Robust regression instead
library(MASS)
rob.mod_S <- rlm(slope_pop$mean.r.recovery~median,data=slope_pop)
summary(rob.mod_S)
# this yields a coefficient estimate & std error but no p-value or R2
library(sfsmisc)
f.robftest(rob.mod_S, var="median") #p= 0.9757


# High S Median vs demography plot
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



