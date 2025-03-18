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
y5 <- read_csv("data/trait_data/traits_anstett2021.csv") %>% 
  mutate(SLA_scaled = scale(SLA, center=T, scale=T),
         FT_scaled = scale(Experiment_Date, center=T, scale=T),
         WC_scaled = scale(Water_Content, center=T, scale=T),
         SC_scaled = scale(Stomatal_Conductance, center=T, scale=T),
         A_scaled = scale(Assimilation, center=T, scale=T),
         Year_scaled = scale(Year, center=T, scale=T))

#Set factors
y5$Block <- as.factor(y5$Block) ; y5$Family <- as.factor(y5$Family) # prep factors

# Set up vectors with treatment and regional information
treatment.v<-c("W", "D")
site.v<-c("S02", "S07", "S08", "S10", "S15", "S16", "S17", "S18", "S29", "S32", "S36")
order.row<-1
slope.reg<-matrix(nrow=length(treatment.v)*length(site.v)*5, ncol=5)

###################################################################################
###################################################################################
#Get slopes of SLA Vs Year
hist(y5$SLA_scaled)

fullmod.SLA <- lmer(SLA_scaled ~ Site*Year_scaled*Drought + (1|Family) + (1|Block),
                    control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)), data=y5)

for (i in 1:length(treatment.v)){
  vis_SLA<-visreg(fullmod.SLA, xvar="Year_scaled", by="Site", cond=list(Drought=treatment.v[i]))
  Res_SLA<-vis_SLA$res
  for (j in 1:length(site.v)){
    Ref_SLA_filter<- Res_SLA %>% filter(Site==site.v[j])
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

#Get slopes of Date of FLowering Vs Year
hist(y5$FT_scaled)

fullmod.FT <- lmer(FT_scaled ~ Site*Year_scaled*Drought + (1|Family) + (1|Block),
                   control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)), data=y5)

for (i in 1:length(treatment.v)){
  vis_FT<-visreg(fullmod.FT, xvar="Year_scaled", by="Site", cond=list(Drought=treatment.v[i]))
  Res_FT<-vis_FT$res
  for (j in 1:length(site.v)){
    Ref_FT_filter<-Res_FT %>% filter(Site==site.v[j])
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


#Get slopes of Water Content Vs Year
hist(y5$WC_scaled)

fullmod.WC <- lmer(WC_scaled ~ Site*Year_scaled*Drought + (1|Family) + (1|Block),
                   control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)), data=y5)

for (i in 1:length(treatment.v)){
  vis_WC<-visreg(fullmod.WC, xvar="Year_scaled", by="Site", cond=list(Drought=treatment.v[i]))
  Res_WC<-vis_WC$res
  for (j in 1:length(site.v)){
    Ref_WC_filter<-Res_WC %>% filter(Site==site.v[j])
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


#Get slopes of Assimilation Vs Year
hist(y5$A_scaled)

fullmod.A <- lmer(A_scaled ~ Site*Year_scaled*Drought + (1|Family) + (1|Block),
                  control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)), data=y5)

for (i in 1:length(treatment.v)){
  vis_A<-visreg(fullmod.A, xvar="Year_scaled", by="Site", cond=list(Drought=treatment.v[i]))
  Res_A<-vis_A$res
  for (j in 1:length(site.v)){
    Ref_A_filter<-Res_A %>% filter(Site==site.v[j])
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

#Get slopes of Stomatal Conductance Vs Year
hist(y5$SC_scaled)

fullmod.gs <- lmer(SC_scaled ~ Site*Year_scaled*Drought  + (1|Family) + (1|Block),
                   control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)), data=y5)

for (i in 1:length(treatment.v)){
  vis_gs<-visreg(fullmod.gs, xvar="Year_scaled", by="Site", cond=list(Drought=treatment.v[i]))
  Res_gs<-vis_gs$res
  for (j in 1:length(site.v)){
    Ref_gs_filter<-Res_gs %>% filter(Site==site.v[j])
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
trait_change$Site <- as.numeric(gsub("S", "", trait_change$Site))
trait_change$Slope=as.numeric(trait_change$Slope)
trait_change$SE_slope=as.numeric(trait_change$SE_slope)

# make cumulative indexes of total trait change towards drought avoidance vs drought escape
# for all traits except flowering time and water content, positive slopes indicate evolution towards drought escape
# for flowering time and water content, negative slope indicates evolution towards drought escape
# so, for index of total trait change, flowering time and water content slopes are subtracted so that bigger values mean greater evolution towards drought escape across all measured traits

trait_change <- trait_change %>% select(-SE_slope) %>% 
  pivot_wider(names_from = c(Trait, Treatment), values_from=Slope) %>% 
  mutate(trait.change.all.dry = SLA_D - FT_D + SC_D + A_D - WC_D,
         trait.change.all.wet = SLA_W - FT_W - WC_W + SC_W + A_W,
         trait.change.best.dry = SLA_D - FT_D,
         trait.change.best.wet = SLA_W - FT_W,
         trait.change.all = trait.change.all.dry + trait.change.all.wet,
         trait.change.best = trait.change.best.dry + trait.change.best.wet)
write_csv(trait_change, "data/trait_data/trait_slopes.csv")

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

trait_geno_pop <- left_join(geno_pop, trait_change, by=c("SiteID"="Site"))

color.list <- trait_geno_pop$Lat.Color


### GENOMIC SELECTION - TRAIT EVOLUTION
# Median S ~ SLA evolution (dry treatment)
ggplot(data=trait_geno_pop, aes(x=Median, y=SLA_D, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards higher SLA (dry treatment") +
  xlab("Median S (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(SLA_D ~ Median, data=trait_geno_pop)) 

# Median S ~ SLA evolution (wet treatment)
ggplot(data=trait_geno_pop, aes(x=Median, y=SLA_W, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards higher SLA (wet treatment") +
  xlab("Median S (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(SLA_W ~ Median, data=trait_geno_pop)) 

# Median S ~ FT evolution (dry treatment)
ggplot(data=trait_geno_pop, aes(x=Median, y=FT_D, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards later flowering (dry treatment") +
  xlab("Median S (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(FT_D ~ Median, data=trait_geno_pop)) 

# Median S ~ FT evolution (wet treatment)
ggplot(data=trait_geno_pop, aes(x=Median, y=FT_W, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  geom_hline(yintercept=1, linetype="dotted") +
  ylab("Rate of evolution towards later flowering (wet treatment") +
  xlab("Median S (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(FT_W ~ Median, data=trait_geno_pop)) 

# Median S ~ all trait evolution (dry treatment)
ggplot(data=trait_geno_pop, aes(x=Median, y=trait.change.all.dry, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards drought avoidance (dry treatment") +
  xlab("Median S (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(trait.change.all.dry ~ Median, data=trait_geno_pop)) 

# Median S ~ all trait evolution (wet treatment)
ggplot(data=trait_geno_pop, aes(x=Median, y=trait.change.all.wet, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards drought avoidance (wet treatment") +
  xlab("Median S (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(trait.change.all.wet ~ Median, data=trait_geno_pop)) 

# Median S ~ best trait evolution (dry treatment)
ggplot(data=trait_geno_pop, aes(x=Median, y=trait.change.best.dry, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards drought avoidance (dry treatment") +
  xlab("Median S (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(trait.change.best.dry ~ Median, data=trait_geno_pop)) 

# Median S ~ best trait evolution (wet treatment)
ggplot(data=trait_geno_pop, aes(x=Median, y=trait.change.best.wet, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards drought avoidance (wet treatment") +
  xlab("Median S (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(trait.change.best.wet ~ Median, data=trait_geno_pop)) 

# Median S ~ all trait evolution
ggplot(data=trait_geno_pop, aes(x=Median, y=trait.change.all, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards drought avoidance") +
  xlab("Median S (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(trait.change.all ~ Median, data=trait_geno_pop)) 

# Median S ~ best trait evolution (both treatments)
ggplot(data=trait_geno_pop, aes(x=Median, y=trait.change.best, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards drought avoidance") +
  xlab("Median S (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(trait.change.best ~ Median, data=trait_geno_pop)) 


### NUCLEOTIDE DIVERSITY - TRAIT EVOLUTION
# Climate Pi ~ SLA evolution (dry treatment)
ggplot(data=trait_geno_pop, aes(x=pi_snp_set, y=SLA_D, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards higher SLA (dry treatment") +
  xlab("Pi (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(SLA_D ~ pi_snp_set, data=trait_geno_pop)) 

# Climate Pi ~ SLA evolution (wet treatment)
ggplot(data=trait_geno_pop, aes(x=pi_snp_set, y=SLA_W, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards higher SLA (wet treatment") +
  xlab("Pi (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(SLA_W ~ pi_snp_set, data=trait_geno_pop)) 

# Climate Pi ~ FT evolution (dry treatment)
ggplot(data=trait_geno_pop, aes(x=pi_snp_set, y=FT_D, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards later flowering (dry treatment") +
  xlab("Pi (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(FT_D ~ pi_snp_set, data=trait_geno_pop)) 

# Climate Pi ~ FT evolution (wet treatment)
ggplot(data=trait_geno_pop, aes(x=pi_snp_set, y=FT_W, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  geom_hline(yintercept=1, linetype="dotted") +
  ylab("Rate of evolution towards later flowering (wet treatment") +
  xlab("Pi (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(FT_W ~ pi_snp_set, data=trait_geno_pop)) 

# Climate Pi ~ all trait evolution (dry treatment)
ggplot(data=trait_geno_pop, aes(x=pi_snp_set, y=trait.change.all.dry, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards drought avoidance (dry treatment") +
  xlab("Pi (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(trait.change.all.dry ~ pi_snp_set, data=trait_geno_pop)) 

# Climate Pi ~ all trait evolution (wet treatment)
ggplot(data=trait_geno_pop, aes(x=pi_snp_set, y=trait.change.all.wet, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  geom_hline(yintercept=1, linetype="dotted") +
  ylab("Rate of evolution towards drought avoidance (wet treatment") +
  xlab("Pi (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(trait.change.all.wet ~ pi_snp_set, data=trait_geno_pop)) 

# Climate Pi ~ best trait evolution (dry treatment)
ggplot(data=trait_geno_pop, aes(x=pi_snp_set, y=trait.change.best.dry, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards drought avoidance (dry treatment") +
  xlab("Pi (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(trait.change.best.dry ~ pi_snp_set, data=trait_geno_pop)) 

# Climate Pi ~ best trait evolution (wet treatment)
ggplot(data=trait_geno_pop, aes(x=pi_snp_set, y=trait.change.best.wet, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards drought avoidance (wet treatment") +
  xlab("Pi (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(trait.change.best.wet ~ pi_snp_set, data=trait_geno_pop)) 

# Climate Pi ~ all trait evolution
ggplot(data=trait_geno_pop, aes(x=pi_snp_set, y=trait.change.all, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards drought avoidance") +
  xlab("Pi (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(trait.change.all ~ pi_snp_set, data=trait_geno_pop)) 

# Climate Pi ~ best trait evolution (both treatments)
ggplot(data=trait_geno_pop, aes(x=pi_snp_set, y=trait.change.best, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards drought avoidance") +
  xlab("Pi (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(trait.change.best ~ pi_snp_set, data=trait_geno_pop)) 

# All Pi ~ SLA evolution (dry treatment)
ggplot(data=trait_geno_pop, aes(x=pi_all_snps, y=SLA_D, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards higher SLA (dry treatment") +
  xlab("Pi (all loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(SLA_D ~ pi_all_snps, data=trait_geno_pop)) 

# All Pi ~ SLA evolution (wet treatment)
ggplot(data=trait_geno_pop, aes(x=pi_all_snps, y=SLA_W, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards higher SLA (wet treatment") +
  xlab("Pi (all loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(SLA_W ~ pi_all_snps, data=trait_geno_pop)) 

# All Pi ~ FT evolution (dry treatment)
ggplot(data=trait_geno_pop, aes(x=pi_all_snps, y=FT_D, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards later flowering (dry treatment") +
  xlab("Pi (all loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(FT_D ~ pi_all_snps, data=trait_geno_pop)) 

# All Pi ~ FT evolution (wet treatment)
ggplot(data=trait_geno_pop, aes(x=pi_all_snps, y=FT_W, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  geom_hline(yintercept=1, linetype="dotted") +
  ylab("Rate of evolution towards later flowering (wet treatment") +
  xlab("Pi (all loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(FT_W ~ pi_all_snps, data=trait_geno_pop)) 

# All Pi ~ all trait evolution (dry treatment)
ggplot(data=trait_geno_pop, aes(x=pi_all_snps, y=trait.change.all.dry, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards drought avoidance (dry treatment") +
  xlab("Pi (all loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(trait.change.all.dry ~ pi_all_snps, data=trait_geno_pop)) 

# All Pi ~ all trait evolution (wet treatment)
ggplot(data=trait_geno_pop, aes(x=pi_all_snps, y=trait.change.all.wet, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards drought avoidance (wet treatment") +
  xlab("Pi (all loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(trait.change.all.wet ~ pi_all_snps, data=trait_geno_pop)) 

# All Pi ~ best trait evolution (dry treatment)
ggplot(data=trait_geno_pop, aes(x=pi_all_snps, y=trait.change.best.dry, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards drought avoidance (dry treatment") +
  xlab("Pi (all loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(trait.change.best.dry ~ pi_all_snps, data=trait_geno_pop)) 

# All Pi ~ best trait evolution (wet treatment)
ggplot(data=trait_geno_pop, aes(x=pi_all_snps, y=trait.change.best.wet, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards drought avoidance (wet treatment") +
  xlab("Pi (all loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(trait.change.best.wet ~ pi_all_snps, data=trait_geno_pop)) 

# All Pi ~ all trait evolution
ggplot(data=trait_geno_pop, aes(x=pi_all_snps, y=trait.change.all, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards drought avoidance") +
  xlab("Pi (all loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(trait.change.all ~ pi_all_snps, data=trait_geno_pop)) 

# All Pi ~ best trait evolution (both treatments)
ggplot(data=trait_geno_pop, aes(x=pi_all_snps, y=trait.change.best, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  ylab("Rate of evolution towards drought avoidance") +
  xlab("Pi (all loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(trait.change.best ~ pi_all_snps, data=trait_geno_pop)) 



### LAMBDA - TRAIT EVOLUTION
# Recovery lambda ~ SLA evolution (dry treatment)
ggplot(data=trait_geno_pop, aes(x=SLA_D, y=mean.lambda.recovery, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  geom_hline(yintercept=1, linetype="dotted") +
  xlab("Rate of evolution towards higher SLA (dry treatment") +
  ylab("Rate of population recovery (lambda)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(mean.lambda.recovery ~ SLA_D, data=trait_geno_pop)) 

# Recovery lambda ~ SLA evolution (wet treatment)
ggplot(data=trait_geno_pop, aes(x=SLA_W, y=mean.lambda.recovery, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  geom_hline(yintercept=1, linetype="dotted") +
  xlab("Rate of evolution towards higher SLA (dry treatment") +
  ylab("Rate of population recovery (lambda)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(mean.lambda.recovery ~ SLA_W, data=trait_geno_pop)) 

# Recovery lambda ~ FT evolution (dry treatment)
ggplot(data=trait_geno_pop, aes(x=FT_D, y=mean.lambda.recovery, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  geom_hline(yintercept=1, linetype="dotted") +
  xlab("Rate of evolution towards later flowering (dry treatment") +
  ylab("Rate of population recovery (lambda)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(mean.lambda.recovery ~ FT_D, data=trait_geno_pop)) 

# Recovery lambda ~ SLA evolution (wet treatment)
ggplot(data=trait_geno_pop, aes(x=SLA_W, y=mean.lambda.recovery, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  geom_hline(yintercept=1, linetype="dotted") +
  xlab("Rate of evolution towards later flowering (wet treatment") +
  ylab("Rate of population recovery (lambda)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(mean.lambda.recovery ~ SLA_W, data=trait_geno_pop)) 

# Recovery lambda ~ all trait evolution (dry treatment)
ggplot(data=trait_geno_pop, aes(x=trait.change.all.dry, y=mean.lambda.recovery, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  geom_hline(yintercept=1, linetype="dotted") +
  xlab("Rate of evolution towards drought escape (dry treatment") +
  ylab("Rate of population recovery (lambda)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(mean.lambda.recovery ~ trait.change.all.dry, data=trait_geno_pop)) 

# Recovery lambda ~ total trait evolution, SLA-FT (dry treatment)
ggplot(data=trait_geno_pop, aes(x=trait.change.best.dry, y=mean.lambda.recovery, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  geom_hline(yintercept=1, linetype="dotted") +
  xlab("Rate of evolution towards drought escape (dry treatment") +
  ylab("Rate of population recovery (lambda)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(mean.lambda.recovery ~ trait.change.best.dry, data=trait_geno_pop)) 

# Recovery lambda ~ total trait evolution (wet treatments)
ggplot(data=trait_geno_pop, aes(x=trait.change.all.wet, y=mean.lambda.recovery, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  geom_hline(yintercept=1, linetype="dotted") +
  xlab("Rate of evolution towards drought escape (wet treatment)") +
  ylab("Rate of population recovery (lambda)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(mean.lambda.recovery ~ trait.change.all.wet, data=trait_geno_pop))  

# Recovery lambda ~ total trait evolution, SLA-FT (wet treatments)
ggplot(data=trait_geno_pop, aes(x=trait.change.best.wet, y=mean.lambda.recovery, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  geom_hline(yintercept=1, linetype="dotted") +
  xlab("Rate of evolution towards drought escape (wet treatment)") +
  ylab("Rate of population recovery (lambda)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(mean.lambda.recovery ~ trait.change.best.wet, data=trait_geno_pop))  

# Recovery lambda ~ total trait evolution (both treatments)
ggplot(data=trait_geno_pop, aes(x=trait.change.all, y=mean.lambda.recovery, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  geom_hline(yintercept=1, linetype="dotted") +
  xlab("Rate of evolution towards drought escape") +
  ylab("Rate of population recovery (lambda)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(mean.lambda.recovery ~ trait.change.all, data=trait_geno_pop))  

# Recovery lambda ~ total trait evolution, SLA-FT (both treatments)
ggplot(data=trait_geno_pop, aes(x=trait.change.best, y=mean.lambda.recovery, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  geom_hline(yintercept=1, linetype="dotted") +
  xlab("Rate of evolution towards drought escape") +
  ylab("Rate of population recovery (lambda)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(mean.lambda.recovery ~ trait.change.best, data=trait_geno_pop))  




### DENSITY - TRAIT EVOLUTON
# Density convexity ~ total trait evolution (dry treatment)
ggplot(data=trait_geno_pop, aes(x=trait.change.all.dry, y=curves.dens.quad, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  geom_hline(yintercept=0, linetype="dotted") +
  xlab("Rate of evolution towards drought escape") +
  ylab("Density convexity") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(curves.dens.quad ~ trait.change.all.dry, data=trait_geno_pop))  

# Density convexity ~ total trait evolution (wet treatment)
ggplot(data=trait_geno_pop, aes(x=trait.change.all.wet, y=curves.dens.quad, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dotted") +
  geom_hline(yintercept=0, linetype="dotted") +
  xlab("Rate of evolution towards drought escape") +
  ylab("Density convexity") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(curves.dens.quad ~ trait.change.all.wet, data=trait_geno_pop))  
