##################################################
# Trait change over time

# This script recalculates trait evolution from Anstett et al. 2021 Evolution Letters
# and then relates it to population dynamics and genomic metrics

# Last updated 2025 03 12
##################################################


###################################################################################
#Import libraries
library(tidyverse)
library(lme4)
library(lmtest)
library(car)
library(visreg)
library(RColorBrewer)
###################################################################################


###################################################################################
#Import main trait dataset from resurrection common garden in greenhouse
y5 <- read_csv("Data/traits_anstett2021.csv")

#Set factors
y5$Block <- as.factor(y5$Block) ; y5$Family <- as.factor(y5$Family) # prep factors

# Set up vectors with treatment and regional information
treatment.v<-c("W", "D")
region.v<-c("1.North", "2.Center", "3.South")
order.row<-1


###################################################################################
###################################################################################
#Get slopes of SLA Vs Year
fullmod.SLA <- lmer(SLA ~ Site.Lat*Year*Drought + (1|Family) + (1|Block),
                    control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)), data=y5)
effectsize::standardize_parameters(fullmod.SLA)

for (i in 1:2){
  vis_SLA<-visreg(fullmod.SLA, xvar="Year", by="Site.Lat", cond=list(Drought=treatment.v[i]))
  Res_SLA<-vis_SLA$res
  for (j in 1:3){
    Ref_SLA_filter<- Res_SLA %>% filter(Region==region.v[j])
    Ref_SLA_filter<- Ref_SLA_filter %>% mutate(Res.scale=scale(visregRes))
    lm_SLA<-lm(Res.scale~Year, data=Ref_SLA_filter)
    summary_SLA<-summary(lm_SLA)
    slope.reg[order.row,4]<-summary_SLA$coefficients[2,1]
    slope.reg[order.row,5]<-summary_SLA$coefficients[2,2]
    order.row<-order.row+1
  }
}

#Get slopes of Date of FLowering Vs Year
fullmod.FT <- lmer(Experiment_Date ~ Region*Year*Drought + (1|Family) + (1|Block) + (1|Site.Lat),
                   control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)), data=y5)
for (i in 1:2){
  vis_FT<-visreg(fullmod.FT, xvar="Year", by="Region", cond=list(Drought=treatment.v[i]))
  Res_FT<-vis_FT$res
  for (j in 1:3){
    Ref_FT_filter<-Res_FT %>% filter(Region==region.v[j])
    Ref_FT_filter<- Ref_FT_filter %>% mutate(Res.scale=scale(visregRes))
    lm_FT<-lm(Res.scale~Year, data=Ref_FT_filter)
    summary_FT<-summary(lm_FT)
    slope.reg[order.row,4]<-summary_FT$coefficients[2,1]
    slope.reg[order.row,5]<-summary_FT$coefficients[2,2]
    order.row<-order.row+1
  }
}


#Get slopes of Water Content Vs Year
fullmod.WC <- lmer(Water_Content ~ Region + Year + Drought + (1|Family) + (1|Block) + (1|Site.Lat),
                   control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)), data=y5)

for (i in 1:2){
  vis_WC<-visreg(fullmod.WC, xvar="Year", by="Region", cond=list(Drought=treatment.v[i]))
  Res_WC<-vis_WC$res
  for (j in 1:3){
    Ref_WC_filter<-Res_WC %>% filter(Region==region.v[j])
    Ref_WC_filter<- Ref_WC_filter %>% mutate(Res.scale=scale(visregRes))
    lm_WC<-lm(Res.scale~Year, data=Ref_WC_filter)
    summary_WC<-summary(lm_WC)
    slope.reg[order.row,4]<-summary_WC$coefficients[2,1]
    slope.reg[order.row,5]<-summary_WC$coefficients[2,2]
    order.row<-order.row+1
  }
}


#Get slopes of Assimilation Vs Year
fullmod.A <- lmer(Assimilation ~ Region*Year*Drought + (1|Family) + (1|Block) + (1|Site.Lat),
                  control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)), data=y5)
for (i in 1:2){
  vis_A<-visreg(fullmod.A, xvar="Year", by="Region", cond=list(Drought=treatment.v[i]))
  Res_A<-vis_A$res
  for (j in 1:3){
    Ref_A_filter<-Res_A %>% filter(Region==region.v[j])
    Ref_A_filter<- Ref_A_filter %>% mutate(Res.scale=scale(visregRes))
    lm_A<-lm(Res.scale~Year, data=Ref_A_filter)
    summary_A<-summary(lm_A)
    slope.reg[order.row,4]<-summary_A$coefficients[2,1]
    slope.reg[order.row,5]<-summary_A$coefficients[2,2]
    order.row<-order.row+1
  }
}

#Get slopes of Stomatal Conductance Vs Year
fullmod.gs <- lmer(Stomatal_Conductance ~ Region*Year*Drought  + (1|Family) + (1|Block) + (1|Site.Lat),
                   control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)), data=y5)
for (i in 1:2){
  vis_gs<-visreg(fullmod.gs, xvar="Year", by="Region", cond=list(Drought=treatment.v[i]))
  Res_gs<-vis_gs$res
  for (j in 1:3){
    Ref_gs_filter<-Res_gs %>% filter(Region==region.v[j])
    Ref_gs_filter<- Ref_gs_filter %>% mutate(Res.scale=scale(visregRes))
    lm_gs<-lm(Res.scale~Year, data=Ref_gs_filter)
    summary_gs<-summary(lm_gs)
    slope.reg[order.row,4]<-summary_gs$coefficients[2,1]
    slope.reg[order.row,5]<-summary_gs$coefficients[2,2]
    order.row<-order.row+1
  }
}

###################################################################################
#Add in dummy varaible to organize X-axis into showing Wet and Dry seperately per region
slope.reg<- slope.reg %>% mutate(x_order=paste(Region,Drought,sep="_"))

#Create dummy variable to sort variables in non-aphabtical order
slope.reg<- slope.reg %>% mutate(order.var=ifelse(Variable == "SLA", 1, 
                                                  ifelse(Variable =="FT", 2,
                                                         ifelse(Variable =="WC", 3,
                                                                ifelse(Variable =="A", 4,5)))))

#Set correct order for all levels
slope.reg$Variable<-as.factor(slope.reg$Variable)
slope.reg$Variable<-factor(slope.reg$Variable,levels=c("SLA","FT","WC","A","gs"))
slope.reg$Drought<-as.factor(slope.reg$Drought)
slope.reg$Drought<-factor(slope.reg$Drought,levels=c("Wet","Dry"))
slope.reg$x_order<-as.factor(slope.reg$x_order)
slope.reg$x_order<-factor(slope.reg$x_order,levels=c("1.North_Wet", "1.North_Dry","2.Center_Wet","2.Center_Dry",
                                                     "3.South_Wet","3.South_Dry"))


###################################################################################
###################################################################################
# Plot slopes across traits and treatments
Trait_Labs<-c("SLA"="(A) SLA", "FT"=" (B) Date of Flowering", "WC"="(C) Water Content", 
              "A"="(D) Assimilation","gs"="(E) Stomatal Conductance")
ggplot(slope.reg, aes(x=x_order, y=Slopes, fill=Drought)) + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values= c("Wet"="#006600","Dry"="#FF7700")) +
  scale_y_continuous(name="Slope")+
  geom_errorbar(mapping=aes(ymin=Slopes-Slopes_STDER, ymax=Slopes+Slopes_STDER), width=0.2, size=1)+
  geom_hline(yintercept=0)+
  theme( axis.title.x=element_blank(),
         #axis.ticks.x = element_blank(),
         axis.text.x = element_text(size=16, face="bold", angle=0,hjust=0.08,vjust = 0.5),
         axis.text.y = element_text(size=12,face="bold"),
         axis.title.y = element_text(color="black", size=16,vjust = 2, face="bold",hjust=0.5))+
  facet_wrap(.~Variable,nrow=3,ncol=2,labeller = labeller(Variable=Trait_Labs)) +
  theme(legend.title = element_blank(),legend.text = element_text(size=14,face="bold"),
        strip.background = element_blank(), strip.text.x=element_text(size=14,face="bold",hjust=0.05,vjust=0))+
  scale_x_discrete(labels=c("1.North_Wet" = "North", "1.North_Dry" = "", 
                            "2.Center_Wet" = "Center","2.Center_Dry" = "",
                            "3.South_Wet" = "South","3.South_Dry" = ""))

#Save Fig S3
#ggsave("Slopes_all_traits.pdf", width = 7, height = 7, units = "in") 
##################################################################################





### Import and merge datasets

# Read in trait change from Evol Lett paper
slopes.traits <- read_csv("data/trait_slopes_anstett2021.csv") %>% 
  #select(-X.1, -X) %>% 
  mutate(tot.change.dry = SLA_Dry - Flowering_Dry + Stomatal_Conductance_Dry + Assimilation_Dry,
         tot.change.wet = SLA_Wet - Flowering_Wet + Stomatal_Conductance_Wet + Assimilation_Wet,
         tot.change = tot.change.dry + tot.change.wet)
# values are slopes of trait values over time
# for all traits except flowering time, positive slopes indicate evolution towards drought escape
# for flowering time, negative slope indicates evolution towards drought escape
# so, for index of total trait change, flowering time slope is subtracted so that bigger values mean greater evolution towards drought escape across all measured traits

slopes.traits$Site <- as.numeric(gsub("S", "", slopes.traits$Site))
slopes.traits$Site.Lat <- as.numeric(gsub("_S", "", slopes.traits$Site.Lat))

# Read in demography data (includes metadata)
demog_recovery <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv") %>% mutate(Site=gsub(" ", "", Site))

# Import genetic diversity (pi)
pi_raw <- read_csv("data/genomic_data/raw_pi.csv") 

# Import genomic selection response
slope.summary <- read_csv("data/snp_change_data/mean_median_S.csv") %>% select(Site, Median, Mean)

# Import density trends 
density_trends <- read_csv("data/demography data/density_trends.csv")

# Join data sets
pi_pop <- left_join(demog_recovery, pi_raw, by=c("Paper_ID"="Site")) 
geno_pop <- left_join(pi_pop, slope.summary, by=c("Paper_ID"="Site"))
geno_pop <- left_join(geno_pop, density_trends, by="Site")
trait_geno_pop <- left_join(geno_pop, slopes.traits, by=c("SiteID"="Site"))

color.list <- trait_geno_pop$Lat.Color

# Recovery lambda ~ SLA evolution (dry treatment)
ggplot(data=trait_geno_pop, aes(x=SLA_Dry, y=mean.lambda.recovery, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  geom_hline(yintercept=1, linetype="dotted") +
  xlab("Rate of evolution in SLA") +
  ylab("Rate of population recovery (lambda)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(mean.lambda.recovery ~ SLA_Dry, data=trait_geno_pop))  

# Recovery lambda ~ flowering time evolution (dry treatment)
ggplot(data=trait_geno_pop, aes(x=Flowering_Dry, y=mean.lambda.recovery, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  geom_hline(yintercept=1, linetype="dotted") +
  xlab("Rate of evolution in flowering time") +
  ylab("Rate of population recovery (lambda)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(mean.lambda.recovery ~ SLA_Dry, data=trait_geno_pop)) 

# Recovery lambda ~ total trait evolution (dry treatment)
ggplot(data=trait_geno_pop, aes(x=tot.change.dry, y=mean.lambda.recovery, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  geom_hline(yintercept=1, linetype="dotted") +
  xlab("Rate of evolution towards drought escape (dry treatment") +
  ylab("Rate of population recovery (lambda)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(mean.lambda.recovery ~ tot.change.dry, data=trait_geno_pop)) 

# Recovery lambda ~ total trait evolution (wet treatments)
ggplot(data=trait_geno_pop, aes(x=tot.change.wet, y=mean.lambda.recovery, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  geom_hline(yintercept=1, linetype="dotted") +
  xlab("Rate of evolution towards drought escape (wet treatment)") +
  ylab("Rate of population recovery (lambda)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(mean.lambda.recovery ~ tot.change.wet, data=trait_geno_pop))  

# Recovery lambda ~ total trait evolution (both treatments)
ggplot(data=trait_geno_pop, aes(x=tot.change, y=mean.lambda.recovery, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  geom_hline(yintercept=1, linetype="dotted") +
  xlab("Rate of evolution towards drought escape") +
  ylab("Rate of population recovery (lambda)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(mean.lambda.recovery ~ tot.change, data=trait_geno_pop))  







# Density convexity ~ leaf evolution (dry treatment)
ggplot(data=trait_geno_pop, aes(x=SLA_Dry, y=curves.dens.quad, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  geom_hline(yintercept=0, linetype="dotted") +
  xlab("Rate of evolution in SLA") +
  ylab("Density convexity") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(curves.dens.quad ~ SLA_Dry, data=trait_geno_pop))  

# Density convexity ~ flowering time evolution (dry treatment)
ggplot(data=trait_geno_pop, aes(x=Flowering_Dry, y=curves.dens.quad, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  geom_hline(yintercept=0, linetype="dotted") +
  xlab("Rate of evolution in flowering time") +
  ylab("Density convexity") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(curves.dens.quad ~ Flowering_Dry, data=trait_geno_pop))  
