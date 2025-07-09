##################################################
### Phenotypic trait change over time

# This script recalculates trait evolution from Anstett et al. 2021 Evolution Letters and then relates it to genomic metrics of selection on climate-associated SNP

# Anstett et al. 2021 (a) focused on regional groups of populations rather than population-specific trends and (b) did not centre and scale variables. Here we need population-specific metrics of trait change, in standardized units that allow us to compare and combine across traits.

# Goal is to see if trait evolution towards drought escape versus drought avoidance can help clarify why some populations have selection against heat/drought-associated alleles (negative S) while some populations have selection for heat/drought-associated alleles (positive S)

# Last updated 2025 06 29
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
library(ggeffects)
###################################################################################


###################################################################################
# Import main trait dataset from resurrection common garden in greenhouse
y5 <- read_csv("data/trait_data/traits_anstett2021.csv") %>% 
  mutate(SLA_scaled = scale(SLA, center=T, scale=T),
         FT_scaled = scale(Experiment_Date, center=T, scale=T),
         WC_scaled = scale(Water_Content, center=T, scale=T),
         SC_scaled = scale(Stomatal_Conductance, center=T, scale=T),
         A_scaled = scale(Assimilation, center=T, scale=T),
         WUE_scaled = scale(Assimilation/Stomatal_Conductance, center=T, scale=T),
         Year_scaled = scale(Year, center=T, scale=T),
         SiteName = gsub(" ", "", Site.Name))

# Set factors
y5$Block <- as.factor(y5$Block) ; y5$Family <- as.factor(y5$Family) # prep factors

# Set up vectors with treatment and regional information
treatment.v<-c("W", "D")
site.v<-c("SweetwaterRiver", "WestForkMojaveRiver", "RedwoodCreek", "NorthForkMiddleForkTule", "RockCreek", "O'NeilCreek", "DeepCreek", "LittleJamesonCreek", "OregonCreek", "Wawona", "DeerCreek")
order.row<-1
slope.reg<-matrix(nrow=length(treatment.v)*length(site.v)*6, ncol=5)

###################################################################################


###################################################################################
### Fit models to get slopes of change over time for each population in each treatment

# Get slopes of SLA Vs Year for each population in each treatment
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

# Get slopes of Date of FLowering Vs Year for each population in each treatment
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

#Get slopes of Water Use Efficiency Vs Year for each population in each treatment
hist(y5$WUE_scaled)

fullmod.wue <- lmer(WUE_scaled ~ SiteName*Year_scaled*Drought  + (1|Family) + (1|Block),
                   control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)), data=y5)

for (i in 1:length(treatment.v)){
  vis_wue<-visreg(fullmod.wue, xvar="Year_scaled", by="SiteName", cond=list(Drought=treatment.v[i]))
  Res_wue<-vis_wue$res
  for (j in 1:length(site.v)){
    Ref_wue_filter<-Res_wue %>% filter(SiteName==site.v[j])
    lm_wue<-lm(visregRes~Year_scaled, data=Ref_wue_filter)
    summary_wue<-summary(lm_wue)
    slope.reg[order.row,1]<-"WUE"
    slope.reg[order.row,2]<-site.v[j]
    slope.reg[order.row,3]<-treatment.v[i]
    slope.reg[order.row,4]<-summary_wue$coefficients[2,1]
    slope.reg[order.row,5]<-summary_wue$coefficients[2,2]
    order.row<-order.row+1
  }
}


trait_change <- as.data.frame(slope.reg)
colnames(trait_change)=c("Trait", "Site", "Treatment", "Slope", "SE_slope")
trait_change$Slope=as.numeric(trait_change$Slope)
trait_change$SE_slope=as.numeric(trait_change$SE_slope)

###################################################################################

###################################################################################
### Make Cumulative Indices of total trait change towards drought avoidance vs drought escape

# For all traits except flowering time and water content, positive slopes indicate evolution towards drought escape
# For flowering time and water content, negative slope indicates evolution towards drought escape
# So, for index of total trait change, flowering time and water content slopes are subtracted so that bigger values mean greater evolution towards drought escape across all measured traits

trait_change <- trait_change %>% dplyr::select(-SE_slope) %>% 
  pivot_wider(names_from = c(Trait, Treatment), values_from=Slope) %>% 
  mutate(trait.change.all.dry = SLA_D - FT_D + SC_D + A_D - WC_D,
         trait.change.all.wet = SLA_W - FT_W + SC_W + A_W - WC_W,
         trait.change.all = trait.change.all.dry + trait.change.all.wet)
write_csv(trait_change, "data/trait_data/trait_slopes.csv")

###################################################################################

###################################################################################
### Merge trait slopes with population metadata and genomic selection data

# Read in demography data (because includes population metadata)
demog_recovery <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv") %>% mutate(Site=gsub(" ", "", Site))

# Import genomic selection response
slope.summary <- read_csv("data/snp_change_2/mean_median_S_all.csv") %>% dplyr::select(Site, Median, Mean)

# Join data sets
geno_pop <- left_join(demog_recovery, slope.summary, by=c("Paper_ID"="Site"))

trait_geno_pop <- left_join(geno_pop, trait_change, by=c("Site"="Site"))

color.list <- trait_geno_pop$Lat.Color

###################################################################################

###################################################################################
### GENOMIC SELECTION - CUMULATIVE TRAIT EVOLUTION

# Median S ~ cumulative trait evolution (across both treatments)
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
summary(mod) #NS
rob.mod <- rlm(trait.change.all ~ Median, data=trait_geno_pop)
f.robftest(rob.mod, var="Median") #NS

# Median S ~ cumulative trait evolution (within dry treatment)
ggplot(data=trait_geno_pop, aes(x=Median, y=trait.change.all.dry, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dashed", fill="grey75") +
  geom_smooth(method=MASS::rlm, color="black", fill="grey50") +
  ylab("Rate of evolution towards drought escape (dry treatment") +
  xlab("Median S (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(trait.change.all.dry ~ Median, data=trait_geno_pop)) #NS
rob.mod <- rlm(trait.change.all.dry ~ Median, data=trait_geno_pop)
f.robftest(rob.mod, var="Median") #+

# Median S ~ cumulative trait evolution (within wet treatment)
ggplot(data=trait_geno_pop, aes(x=Median, y=trait.change.all.wet, color=as.character(Latitude))) +
  geom_point(size=3) +
  geom_smooth(method="lm", color="black", linetype="dashed", fill="grey75") +
  geom_smooth(method=MASS::rlm, color="black", fill="grey50") +
  ylab("Rate of evolution towards drought escape (wet treatment") +
  xlab("Median S (climate loci)") + 
  scale_color_manual(values=color.list) +
  theme_classic() +
  theme(legend.title = element_blank())

summary(lm(trait.change.all.wet ~ Median, data=trait_geno_pop)) #NS
rob.mod <- rlm(trait.change.all.wet ~ Median, data=trait_geno_pop)
f.robftest(rob.mod, var="Median") #NS

### Conclusion from trait index analyses: weak trend for negative S to be more consistent with escape and positive S to be more consistent with avoidance within dry treatment, but only if influence of Oregon Creek is minimized (and Oregon Creek is one of the 3 populations with significantly negative S)

###################################################################################


###################################################################################
### Focus on populations whose genomic selection was significantly different from null: explore bins of populations with negative vs positive slopes (ignoring the mushy middle)

# Pops whose slopes went significantly positive: 3, 4, 11
# Pops whose slopes went significantly negative: 6, 7, 8
# Pops in the middle whose distributions are not distinguishable from 0: 1, 2, 5, 9, 10
# Instead of regression, look for mean differences in traits between pops(3,4,11) vs pops(6,7,8)

trait_geno_pop_cull <- trait_geno_pop %>% 
  dplyr::filter(Paper_ID==3|Paper_ID==4|Paper_ID==11|Paper_ID==6|Paper_ID==7|Paper_ID==8) %>% 
  droplevels() %>% 
  mutate(Sgroup = ifelse(Median>0, "positive S", "negative S"))

# Cumulative rate of change in traits, both environments
ggplot(data=trait_geno_pop_cull, aes(x=Sgroup, y=trait.change.all)) +
  geom_boxplot() +
  geom_point(col=trait_geno_pop_cull$Lat.Color, size=3) #NS

# Cumulative rate of change in traits, wet treatment
ggplot(data=trait_geno_pop_cull, aes(x=Sgroup, y=trait.change.all.wet)) +
  geom_boxplot() +
  geom_point(col=trait_geno_pop_cull$Lat.Color, size=3)
t.test(trait_geno_pop_cull$trait.change.all.wet[trait_geno_pop_cull$Sgroup=="negative S"],trait_geno_pop_cull$trait.change.all.wet[trait_geno_pop_cull$Sgroup=="positive S"]) #NS

# Cumulative rate of change in traits, dry treatment
ggplot(data=trait_geno_pop_cull, aes(x=Sgroup, y=trait.change.all.dry)) +
  geom_boxplot() +
  geom_point(col=trait_geno_pop_cull$Lat.Color, size=3)
t.test(trait_geno_pop_cull$trait.change.all.dry[trait_geno_pop_cull$Sgroup=="negative S"],trait_geno_pop_cull$trait.change.all.dry[trait_geno_pop_cull$Sgroup=="positive S"])

### Conclusion from binned analyses: no support for relationship between S bin and cumulative trait change, but this test has VERY low power

###################################################################################

###################################################################################
### Cumulative indices assume that we know how traits should cluster into syndromes, but Haley's work shows that traits can change in complicated ways that don't map neatly onto simple escape-avoid dichotomy

# Use multiple regression to test whether trait change can collectively predict selection response, allowing each trait to have its own direction of contribution

trait_geno_pop_na <- trait_geno_pop %>% drop_na(SLA_D) %>% dplyr::select(Median, SLA_D, FT_D, A_D, SC_D, WC_D)

#Export file
write_csv(trait_geno_pop_na,"data/trait_data/trait_pop.csv")

mod.dry <- lm(Median ~ SLA_D + FT_D + A_D + SC_D + WC_D, dat=trait_geno_pop_na, na.action=na.fail)
summary(mod.dry) #positive effect of increased carbon assimilation, negative effect of increased stomatal conductance (traits that we were forcing to go in same direction on simple escape-avoid continuum)
plot(mod.dry)
rob.mod.dry <- rlm(Median ~ SLA_D + FT_D + A_D + SC_D + WC_D, dat=trait_geno_pop, maxit=100)
f.robftest(rob.mod.dry, var="SLA_D")
f.robftest(rob.mod.dry, var="FT_D")
f.robftest(rob.mod.dry, var="A_D") #*
f.robftest(rob.mod.dry, var="SC_D") #**
f.robftest(rob.mod.dry, var="WC_D")

# This is way too many variables for such a small dataset, though. Try coherent subsets?

mod.dry.anat <- lm(Median ~ SLA_D + WC_D, dat=trait_geno_pop)
summary(mod.dry.anat) #NS  
rob.mod.dry.anat <- rlm(Median ~ SLA_D + WC_D, dat=trait_geno_pop)
f.robftest(rob.mod.dry.anat, var="SLA_D") #NS
f.robftest(rob.mod.dry.anat, var="WC_D") #NS

mod.dry.gasx <- lm(Median ~ A_D + SC_D, dat=trait_geno_pop)
summary(mod.dry.gasx) #positive effects of evolution towards greater carbon assimilation, decreased stomatal conductance
rob.mod.dry.gasx <- rlm(Median ~ A_D + SC_D, dat=trait_geno_pop)
f.robftest(rob.mod.dry.gasx, var="A_D") #*
f.robftest(rob.mod.dry.gasx, var="SC_D") #***

mod.dry.phen <- lm(Median ~ FT_D, dat=trait_geno_pop)
summary(mod.dry.phen) #NS
rob.mod.dry.phen <- rlm(Median ~ FT_D, dat=trait_geno_pop)
f.robftest(rob.mod.dry.phen, var="FT_D") #NS

mod.wet <- lm(Median ~ SLA_W + FT_W + A_W + SC_W + WC_W, dat=trait_geno_pop)
summary(mod.wet) #positive effects of change towards increased SLA, earlier FT, greater carbon assimilation, and higher leaf water content on direction of selection
plot(mod.wet)
rob.mod.wet <- rlm(Median ~ SLA_W + FT_W + A_W + SC_W + WC_W, dat=trait_geno_pop, maxit=100)
f.robftest(rob.mod.wet, var="SLA_W") #****
f.robftest(rob.mod.wet, var="FT_W") #****
f.robftest(rob.mod.wet, var="A_W") #****
f.robftest(rob.mod.wet, var="SC_W") #****
f.robftest(rob.mod.wet, var="WC_W") #****

# Again, though, too many variables for the data set; try subsets
mod.wet.anat <- lm(Median ~ SLA_W + WC_W, dat=trait_geno_pop)
summary(mod.wet.anat) #NS  
rob.mod.wet.anat <- rlm(Median ~ SLA_W + WC_W, dat=trait_geno_pop)
f.robftest(rob.mod.wet.anat, var="SLA_W") #NS
f.robftest(rob.mod.wet.anat, var="WC_W") #NS

mod.wet.gasx <- lm(Median ~ A_W + SC_W, dat=trait_geno_pop)
summary(mod.wet.gasx) #NS
rob.mod.wet.gasx <- rlm(Median ~ A_W + SC_W, dat=trait_geno_pop)
f.robftest(rob.mod.wet.gasx, var="A_W") #NS
f.robftest(rob.mod.wet.gasx, var="SC_W") #NS

mod.wet.phen <- lm(Median ~ FT_W, dat=trait_geno_pop)
summary(mod.wet.phen) #NS
rob.mod.wet.phen <- rlm(Median ~ FT_W, dat=trait_geno_pop)
f.robftest(rob.mod.wet.phen, var="FT_W") #NS

### Conclusion from multiple regression analyses: within dry treatment, evolution of gas exchange can explain some differences in genomic selection. Specifically, positive S is associated with evolution towards increased photosynthesis and decreased stomatal conductance (sounds adaptive) and negative S is associated with evolution towards decreased photosynthesis and increased conductance (sounds maladaptive). Results for wet treatment are spurious in an overdetermined model.

# Alternative route for exploring multiple regression: dredge on all subsets?
library(MuMIn)
dd <- dredge (mod.dry, evaluate=TRUE, m.lim=c(1, 2))

# Visualize partial effects in gas exchange model
pred_df_A <- predict_response(mod.dry.gasx, terms=c("A_D","SC_D"), margin="mean_reference")
pred_A_plot <- ggplot(filter(pred_df_A, group==-0.07), aes(x=x, y=predicted)) + #, colour=group
  geom_ribbon(aes(x=x, ymin=conf.low, ymax=conf.high), fill="grey75") + #, group=group, colour=group
  stat_smooth(method="lm", color="black") +
  xlab("Evolution of photosynthetic rate") + 
  ylab("Predicted Median S") +
  theme_classic()

pred_df_B <- predict_response(mod.dry.gasx, terms=c("SC_D","A_D"), margin="mean_reference")
pred_B_plot <- ggplot(filter(pred_df_B, group==0.01), aes(x=x, y=predicted)) + #, colour=group
  geom_ribbon(aes(x=x, ymin=conf.low, ymax=conf.high), fill="grey75") + #, group=group, colour=group
  stat_smooth(method="lm", color="black") +
  xlab("Evolution of stomatal conductance") + 
  ylab("Predicted Median S") +
  theme_classic()

# Supplemental Figure Sx
plot_grid(pred_A_plot, pred_B_plot)
ggsave("Graphs/Traits_Selection.pdf")

###################################################################################


###################################################################################
### MANOVA because arguably this is a more sensible direction to test things (genotype to phenotype)...but cannot support this many variables with such a small dataset
man.dry <- manova(cbind(SLA_D, FT_D, A_D, SC_D, WC_D) ~ Median, data=trait_geno_pop)
summary(man.dry) #*

man.wet <- manova(cbind(SLA_W, FT_W, A_W, SC_W, WC_W) ~ Median, data=trait_geno_pop)
summary(man.wet) #NS

### Conclusion from MANOVA analyses: not sure this is robust but still indicating support for effect of evolution in dry treatment.

###################################################################################


###################################################################################
### Collapse rates of trait evolution into multivariate PC axes?

# Both treatments PCA
pc.traits <- princomp(~SLA_D+SLA_W+FT_D+FT_W+A_D+A_W+SC_D+SC_W+WC_D+WC_W, data=trait_geno_pop, cor=TRUE, na.action=na.exclude)
summary(pc.traits) #PC1&2 only get to 65% variance explained
pc.traits$loadings
biplot(pc.traits) #dry and wet plasticity in rates of change and relationships among different traits
pc.traits$scores
trait_geno_pop <- cbind(trait_geno_pop, pc.traits$scores)

mod.pc <- lm(Median ~ Comp.1+Comp.2, data=trait_geno_pop)
summary(mod.pc) #NS

# Dry treatment PCA
pc.traits.dry <- princomp(~SLA_D+FT_D+A_D+SC_D+WC_D, data=trait_geno_pop, cor=TRUE, na.action=na.exclude)
summary(pc.traits.dry) #PC1&2 get to 83% variance explained
pc.traits.dry$loadings
biplot(pc.traits.dry) #SLA and WC basically give the same info (in bivariate space of first 2 PC axes), other traits not coupled quite as fully as simple avoid-escape dichotomy would suggest
pc.traits.dry$scores
colnames(pc.traits.dry$scores) = paste(colnames(pc.traits.dry$scores), ".dry", sep="")
trait_geno_pop <- cbind(trait_geno_pop, pc.traits.dry$scores)

mod.pc.dry <- lm(Median ~ Comp.1.dry+Comp.2.dry, data=trait_geno_pop)
summary(mod.pc.dry) #NS

# Wet treatment PCA
pc.traits.wet <- princomp(~SLA_W+FT_W+A_W+SC_W+WC_W, data=trait_geno_pop, cor=TRUE, na.action=na.exclude)
summary(pc.traits.wet) #PC1&2 get to 80% variance explained
pc.traits.wet$loadings
biplot(pc.traits.wet) #SLA and WC basically give the same info (in bivariate space of first 2 PC axes), other traits not coupled quite as fully as simple avoid-escape dichotomy would suggest
pc.traits.wet$scores
colnames(pc.traits.wet$scores) = paste(colnames(pc.traits.wet$scores), ".wet", sep="")
trait_geno_pop <- cbind(trait_geno_pop, pc.traits.wet$scores)

mod.pc.wet <- lm(Median ~ Comp.1.wet+Comp.2.wet, data=trait_geno_pop)
summary(mod.pc.wet) #NS

### Conclusion from PCA analyses: This attempt at distillation is not helping to produce easily interpretable biological axes. Abandon.

###################################################################################

