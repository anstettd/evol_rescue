##################################################
# Trait change over time

# This script recalculates trait evolution from Anstett et al. 2021 Evolution Letters and then relates it to genomic metrics of selection on climate-associated SNP

# Goal is to see if trait evolution towards drought escape versus drought avoidance can help clarify why some populations have selection against heat/drought-associated alleles (negative S) and some populations have selection for heat/drought-associated alleles (positive S)

# Last updated 2025 06 24
##################################################


###################################################################################
#Import libraries
library(tidyverse)
library(lme4)
library(lmtest)
library(car)
library(visreg)
library(RColorBrewer)
library(MASS)
library(sfsmisc)
###################################################################################


###################################################################################
#Import main trait dataset from resurrection common garden in greenhouse
y5 <- read_csv("data/trait_data/traits_anstett2021.csv") %>% 
  mutate(SLA_scaled = scale(SLA, center=T, scale=T),
         FT_scaled = scale(Experiment_Date, center=T, scale=T),
         WC_scaled = scale(Water_Content, center=T, scale=T),
         SC_scaled = scale(Stomatal_Conductance, center=T, scale=T),
         A_scaled = scale(Assimilation, center=T, scale=T),
         Year_scaled = scale(Year, center=T, scale=T),
         SiteName = gsub(" ", "", Site.Name))

#Set factors
y5$Block <- as.factor(y5$Block) ; y5$Family <- as.factor(y5$Family) # prep factors

# Set up vectors with treatment and regional information
treatment.v<-c("W", "D")
site.v<-c("SweetwaterRiver", "WestForkMojaveRiver", "RedwoodCreek", "NorthForkMiddleForkTule", "RockCreek", "O'NeilCreek", "DeepCreek", "LittleJamesonCreek", "OregonCreek", "Wawona", "DeerCreek")
order.row<-1
slope.reg<-matrix(nrow=length(treatment.v)*length(site.v)*5, ncol=5)

###################################################################################
###################################################################################
#Get slopes of SLA Vs Year for each population in each treatment
hist(y5$SLA_scaled)

fullmod.SLA <- lmer(SLA_scaled ~ SiteName*Year_scaled*Drought + (1|Family) + (1|Block),
                    control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)), data=y5)

for (i in 1:length(treatment.v)){
  vis_SLA<-visreg(fullmod.SLA, xvar="Year_scaled", by="SiteName", cond=list(Drought=treatment.v[i]))
  Res_SLA<-vis_SLA$res
  for (j in 1:length(site.v)){
    Ref_SLA_filter<- Res_SLA %>% filter(SiteName==site.v[j])
    lm_SLA<-lm(visregRes~Year_scaled, data=Ref_SLA_filter)
    summary_SLA<-summary(lm_SLA)
    slope.reg[order.row,1]<-"SLA"
    slope.reg[order.row,2]<-site.v[j]
    slope.reg[order.row,3]<-treatment.v[i]
    slope.reg[order.row,4]<-summary_SLA$coefficients[2,1]
    slope.reg[order.row,5]<-summary_SLA$coefficients[2,2]
    order.row<-order.row+1
  }
}

#Get slopes of Date of FLowering Vs Year for each population in each treatment
hist(y5$FT_scaled)

fullmod.FT <- lmer(FT_scaled ~ SiteName*Year_scaled*Drought + (1|Family) + (1|Block),
                   control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)), data=y5)

for (i in 1:length(treatment.v)){
  vis_FT<-visreg(fullmod.FT, xvar="Year_scaled", by="SiteName", cond=list(Drought=treatment.v[i]))
  Res_FT<-vis_FT$res
  for (j in 1:length(site.v)){
    Ref_FT_filter<-Res_FT %>% filter(SiteName==site.v[j])
    lm_FT<-lm(visregRes~Year_scaled, data=Ref_FT_filter)
    summary_FT<-summary(lm_FT)
    slope.reg[order.row,1]<-"FT"
    slope.reg[order.row,2]<-site.v[j]
    slope.reg[order.row,3]<-treatment.v[i]
    slope.reg[order.row,4]<-summary_FT$coefficients[2,1]
    slope.reg[order.row,5]<-summary_FT$coefficients[2,2]
    order.row<-order.row+1
  }
}

#Get slopes of Water Content Vs Year for each population in each treatment
hist(y5$WC_scaled)

fullmod.WC <- lmer(WC_scaled ~ SiteName*Year_scaled*Drought + (1|Family) + (1|Block),
                   control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)), data=y5)

for (i in 1:length(treatment.v)){
  vis_WC<-visreg(fullmod.WC, xvar="Year_scaled", by="SiteName", cond=list(Drought=treatment.v[i]))
  Res_WC<-vis_WC$res
  for (j in 1:length(site.v)){
    Ref_WC_filter<-Res_WC %>% filter(SiteName==site.v[j])
    lm_WC<-lm(visregRes~Year_scaled, data=Ref_WC_filter)
    summary_WC<-summary(lm_WC)
    slope.reg[order.row,1]<-"WC"
    slope.reg[order.row,2]<-site.v[j]
    slope.reg[order.row,3]<-treatment.v[i]
    slope.reg[order.row,4]<-summary_WC$coefficients[2,1]
    slope.reg[order.row,5]<-summary_WC$coefficients[2,2]
    order.row<-order.row+1
  }
}

#Get slopes of Assimilation Vs Year for each population in each treatment
hist(y5$A_scaled)

fullmod.A <- lmer(A_scaled ~ SiteName*Year_scaled*Drought + (1|Family) + (1|Block),
                  control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)), data=y5)

for (i in 1:length(treatment.v)){
  vis_A<-visreg(fullmod.A, xvar="Year_scaled", by="SiteName", cond=list(Drought=treatment.v[i]))
  Res_A<-vis_A$res
  for (j in 1:length(site.v)){
    Ref_A_filter<-Res_A %>% filter(SiteName==site.v[j])
    lm_A<-lm(visregRes~Year_scaled, data=Ref_A_filter)
    summary_A<-summary(lm_A)
    slope.reg[order.row,1]<-"A"
    slope.reg[order.row,2]<-site.v[j]
    slope.reg[order.row,3]<-treatment.v[i]
    slope.reg[order.row,4]<-summary_A$coefficients[2,1]
    slope.reg[order.row,5]<-summary_A$coefficients[2,2]
    order.row<-order.row+1
  }
}

#Get slopes of Stomatal Conductance Vs Year for each population in each treatment
hist(y5$SC_scaled)

fullmod.gs <- lmer(SC_scaled ~ SiteName*Year_scaled*Drought  + (1|Family) + (1|Block),
                   control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)), data=y5)

for (i in 1:length(treatment.v)){
  vis_gs<-visreg(fullmod.gs, xvar="Year_scaled", by="SiteName", cond=list(Drought=treatment.v[i]))
  Res_gs<-vis_gs$res
  for (j in 1:length(site.v)){
    Ref_gs_filter<-Res_gs %>% filter(SiteName==site.v[j])
    lm_gs<-lm(visregRes~Year_scaled, data=Ref_gs_filter)
    summary_gs<-summary(lm_gs)
    slope.reg[order.row,1]<-"SC"
    slope.reg[order.row,2]<-site.v[j]
    slope.reg[order.row,3]<-treatment.v[i]
    slope.reg[order.row,4]<-summary_gs$coefficients[2,1]
    slope.reg[order.row,5]<-summary_gs$coefficients[2,2]
    order.row<-order.row+1
  }
}

trait_change <- as.data.frame(slope.reg)
colnames(trait_change)=c("Trait", "Site", "Treatment", "Slope", "SE_slope")
trait_change$Slope=as.numeric(trait_change$Slope)
trait_change$SE_slope=as.numeric(trait_change$SE_slope)


### Make Cumulative Indices of total trait change towards drought avoidance vs drought escape
# for all traits except flowering time and water content, positive slopes indicate evolution towards drought escape
# for flowering time and water content, negative slope indicates evolution towards drought escape
# so, for index of total trait change, flowering time and water content slopes are subtracted so that bigger values mean greater evolution towards drought escape across all measured traits

trait_change <- trait_change %>% select(-SE_slope) %>% 
  pivot_wider(names_from = c(Trait, Treatment), values_from=Slope) %>% 
  mutate(trait.change.all.dry = SLA_D - FT_D + SC_D + A_D - WC_D,
         trait.change.all.wet = SLA_W - FT_W + SC_W + A_W - WC_W,
         trait.change.best.dry = SLA_D - FT_D,
         trait.change.best.wet = SLA_W - FT_W,
         trait.change.all = trait.change.all.dry + trait.change.all.wet,
         trait.change.best = trait.change.best.dry + trait.change.best.wet)
write_csv(trait_change, "data/trait_data/trait_slopes.csv")

# Read in demography data (includes metadata)
demog_recovery <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv") %>% mutate(Site=gsub(" ", "", Site))

# Import genomic selection response
slope.summary <- read_csv("data/snp_change_2/mean_median_S_all.csv") %>% select(Site, Median, Mean)

# Join data sets
geno_pop <- left_join(demog_recovery, slope.summary, by=c("Paper_ID"="Site"))

trait_geno_pop <- left_join(geno_pop, trait_change, by=c("Site"="Site"))

color.list <- trait_geno_pop$Lat.Color


### GENOMIC SELECTION - TRAIT EVOLUTION

# Median S ~ all trait evolution (both treatments)
ggplot(data=trait_geno_pop, aes(x=Median, y=trait.change.all, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dashed", fill="grey75") +
  geom_smooth(method=MASS::rlm, color="black", fill="grey50") +
  ylab("Rate of evolution towards drought escape") +
  xlab("Median S (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

mod<-lm(trait.change.all ~ Median, data=trait_geno_pop)
summary(mod)
rob.mod <- rlm(trait.change.all ~ Median, data=trait_geno_pop)
f.robftest(rob.mod, var="Median")

# Median S ~ best trait evolution (both treatments)
ggplot(data=trait_geno_pop, aes(x=Median, y=trait.change.best, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dashed", fill="grey75") +
  geom_smooth(method=MASS::rlm, color="black", fill="grey50") +
  ylab("Rate of evolution towards drought escape") +
  xlab("Median S (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

mod<-lm(trait.change.best ~ Median, data=trait_geno_pop)
summary(mod)
rob.mod <- rlm(trait.change.best ~ Median, data=trait_geno_pop)
f.robftest(rob.mod, var="Median")

# Median S ~ all trait evolution (dry treatment)
ggplot(data=trait_geno_pop, aes(x=Median, y=trait.change.all.dry, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dashed", fill="grey75") +
  geom_smooth(method=MASS::rlm, color="black", fill="grey50") +
  ylab("Rate of evolution towards drought escape (dry treatment") +
  xlab("Median S (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(trait.change.all.dry ~ Median, data=trait_geno_pop)) 
rob.mod <- rlm(trait.change.all.dry ~ Median, data=trait_geno_pop)
f.robftest(rob.mod, var="Median")

# Median S ~ all trait evolution (wet treatment)
ggplot(data=trait_geno_pop, aes(x=Median, y=trait.change.all.wet, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dashed", fill="grey75") +
  geom_smooth(method=MASS::rlm, color="black", fill="grey50") +
  ylab("Rate of evolution towards drought escape (wet treatment") +
  xlab("Median S (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())


# Median S ~ best trait evolution (dry treatment)
ggplot(data=trait_geno_pop, aes(x=Median, y=trait.change.best.dry, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dashed", fill="grey75") +
  geom_smooth(method=MASS::rlm, color="black", fill="grey50") +
  ylab("Rate of evolution towards drought escape (dry treatment") +
  xlab("Median S (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(trait.change.best.dry ~ Median, data=trait_geno_pop)) 
rob.mod <- rlm(trait.change.best.dry ~ Median, data=trait_geno_pop)
f.robftest(rob.mod, var="Median")

# Median S ~ best trait evolution (wet treatment)
ggplot(data=trait_geno_pop, aes(x=Median, y=trait.change.best.wet, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dashed", fill="grey75") +
  geom_smooth(method=MASS::rlm, color="black", fill="grey50") +
  ylab("Rate of evolution towards drought escape (wet treatment") +
  xlab("Median S (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())





#pops whose slopes went significantly positive: 3, 4, 11
#pops whose slopes went significantly negative: 6, 7, 8
#pops in the middle whose distributions are not distinguishable from 0: 1, 2, 5, 9, 10
#instead of regression, look for mean differences in traits between pops(3,4,11) vs pops(6,7,8)

trait_geno_pop_cull <- trait_geno_pop %>% 
  dplyr::filter(Paper_ID==3|Paper_ID==4|Paper_ID==11|Paper_ID==6|Paper_ID==7|Paper_ID==8) %>% 
  droplevels() %>% 
  mutate(Sgroup = ifelse(Median>0, "positive S", "negative S"))

# All traits, both environments
ggplot(data=trait_geno_pop_cull, aes(x=Sgroup, y=trait.change.all)) +
  geom_boxplot() +
  geom_point(col=trait_geno_pop_cull$Lat.Color, size=3)

# Best traits, both environments
ggplot(data=trait_geno_pop_cull, aes(x=Sgroup, y=trait.change.best)) +
  geom_boxplot() +
  geom_point(col=trait_geno_pop_cull$Lat.Color, size=3)
t.test(trait_geno_pop_cull$trait.change.best[trait_geno_pop_cull$Sgroup=="negative S"],trait_geno_pop_cull$trait.change.best[trait_geno_pop_cull$Sgroup=="positive S"])

# All traits, wet treatment
ggplot(data=trait_geno_pop_cull, aes(x=Sgroup, y=trait.change.all.wet)) +
  geom_boxplot() +
  geom_point(col=trait_geno_pop_cull$Lat.Color, size=3)
t.test(trait_geno_pop_cull$trait.change.all.wet[trait_geno_pop_cull$Sgroup=="negative S"],trait_geno_pop_cull$trait.change.all.wet[trait_geno_pop_cull$Sgroup=="positive S"])

# Best traits, wet treatment
ggplot(data=trait_geno_pop_cull, aes(x=Sgroup, y=trait.change.best.wet)) +
  geom_boxplot() +
  geom_point(col=trait_geno_pop_cull$Lat.Color, size=3)

# All traits, dry treatment
ggplot(data=trait_geno_pop_cull, aes(x=Sgroup, y=trait.change.all.dry)) +
  geom_boxplot() +
  geom_point(col=trait_geno_pop_cull$Lat.Color, size=3)
t.test(trait_geno_pop_cull$trait.change.all.dry[trait_geno_pop_cull$Sgroup=="negative S"],trait_geno_pop_cull$trait.change.all.dry[trait_geno_pop_cull$Sgroup=="positive S"])

# Best traits, dry treatment
ggplot(data=trait_geno_pop_cull, aes(x=Sgroup, y=trait.change.best.dry)) +
  geom_boxplot() +
  geom_point(col=trait_geno_pop_cull$Lat.Color, size=3) 
  

summary(manova(cbind(SLA_D, FT_D, A_D, SC_D, WC_D) ~ Median, data=trait_geno_pop))
summary(manova(cbind(SLA_W, FT_W, A_W, SC_W, WC_W) ~ Median, data=trait_geno_pop))

