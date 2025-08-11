##################################################
### Phenotypic trait change over time

# This script recalculates trait evolution from Anstett et al. 2021 Evolution Letters and then relates it to genomic metrics of selection on climate-associated SNP

# Anstett et al. 2021 (a) focused on regional groups of populations rather than population-specific trends and (b) did not centre and scale variables. Here we need population-specific metrics of trait change, in standardized units that allow us to compare and combine across traits.

# Goal is to test if different directions of trait evolution (e.g. towards drought escape versus drought avoidance) can help clarify why some populations have selection against heat/drought-associated alleles (negative S) while some populations have selection for heat/drought-associated alleles (positive S)

# Last updated 2025 07 29
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
### Merge trait slopes with population metadata and genomic selection data

# Read in demography data (because includes population metadata)
demog_recovery <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv") %>% mutate(Site=gsub(" ", "", Site))

# Import genomic selection response
slope.summary <- read_csv("data/snp_change_2/mean_median_S_all.csv") %>% dplyr::select(Site, Median, Mean)

# Join data sets
geno_pop <- left_join(demog_recovery, slope.summary, by=c("Paper_ID"="Site"))

trait_geno_pop <- left_join(geno_pop, trait_change, by=c("Site"="Site")) %>% 
  dplyr::select(-SE_slope) %>% 
  pivot_wider(names_from = c(Trait, Treatment), values_from=Slope)

trait_geno_pop_na <- trait_geno_pop %>% drop_na(SLA_D) %>% dplyr::select(Median, SLA_D, FT_D, A_D, SC_D, WC_D)

#Export file
write_csv(trait_geno_pop_na,"data/trait_data/trait_pop.csv")

color.list <- trait_geno_pop$Lat.Color

###################################################################################


###################################################################################
### Use multiple regression to test whether trait change can  predict selection response

# We don't have a large enough dataset to include all 5 variables at once. Instead do subsets grouped by function.

# Focus on trait expression in dry treatment (the relevant selective environment)

# Anatomical traits: specific leaf area and leaf water content
mod.dry.anat <- lm(Median ~ SLA_D + WC_D, dat=trait_geno_pop)
summary(mod.dry.anat) #NS  
rob.mod.dry.anat <- rlm(Median ~ SLA_D + WC_D, dat=trait_geno_pop)
f.robftest(rob.mod.dry.anat, var="SLA_D") #NS
f.robftest(rob.mod.dry.anat, var="WC_D") #NS

# Gas exchange traits: photosynthetic carbon assimilation and stomatal conductance
mod.dry.gasx <- lm(Median ~ A_D + SC_D, dat=trait_geno_pop)
summary(mod.dry.gasx) #positive effects of evolution towards greater carbon assimilation, decreased stomatal conductance
rob.mod.dry.gasx <- rlm(Median ~ A_D + SC_D, dat=trait_geno_pop)
f.robftest(rob.mod.dry.gasx, var="A_D") #*
f.robftest(rob.mod.dry.gasx, var="SC_D") #***

# Phenological trait: flowering time
mod.dry.phen <- lm(Median ~ FT_D, dat=trait_geno_pop)
summary(mod.dry.phen) #NS
rob.mod.dry.phen <- rlm(Median ~ FT_D, dat=trait_geno_pop)
f.robftest(rob.mod.dry.phen, var="FT_D") #NS

### Conclusion from multiple regression analyses: within dry treatment, evolution of gas exchange can explain some differences in genomic selection. Specifically, positive S is associated with evolution towards increased photosynthesis and decreased stomatal conductance (sounds adaptive) and negative S is associated with evolution towards decreased photosynthesis and increased conductance (seems maladaptive). 

###################################################################################
# Visualize partial effects in gas exchange model
###################################################################################
#Evolution of Photosynthetic Rate

pred_df_A <- predict_response(mod.dry.gasx, terms=c("A_D","SC_D"), margin="mean_reference")
plot(pred_df_A, show_data=TRUE)

raw_points <- attr(pred_df_A, "rawdata")
#Raw plot
pred_df_A$group <- factor(pred_df_A$group)
raw_points$group <- factor(raw_points$group, levels = levels(pred_df_A$group))

#Plot different levels of other variable
pred_A_gg <-ggplot(pred_df_A, aes(x = x, y = predicted, color = group,fill = group)) + 
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, group = group), alpha = 0.1, color = NA) +
  geom_point(data = raw_points, aes(x = x, y = response), shape=21, color = "black", fill="black", size=1.5,inherit.aes = FALSE) +
  xlab("Evolution of Photosynthetic Rate") + 
  ylab("Predicted Median S") +
  theme_classic()+ theme(
    axis.text.x = element_text(size=16, face="bold"),
    axis.text.y = element_text(size=18,face="bold"),
    axis.title.x = element_text(color="black", size=16, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=18,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.position = "none"
  )


#Single regression plot
pred_A_plot <- ggplot(filter(pred_df_A, group==-0.07), aes(x=x, y=predicted)) + #, colour=group
  geom_ribbon(aes(x=x, ymin=conf.low, ymax=conf.high), fill="grey75") + #, group=group, colour=group
  stat_smooth(method="lm", color="black") +
  xlab("Evolution of Photosynthetic Rate") + 
  ylab("Predicted Median S") +
  theme_classic()
###################################################################################
##Evolution of Stomatal Conductance

#Raw plot
pred_df_B <- predict_response(mod.dry.gasx, terms=c("SC_D","A_D"), margin="mean_reference")
plot(pred_df_B, show_data=TRUE)

raw_points <- attr(pred_df_B, "rawdata")
#Raw plot
pred_df_B$group <- factor(pred_df_B$group)
raw_points$group <- factor(raw_points$group, levels = levels(pred_df_B$group))

#Plot different levels of other variable
pred_B_gg <-ggplot(pred_df_B, aes(x = x, y = predicted, color = group,fill = group)) + 
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, group = group), alpha = 0.1, color = NA) +
  geom_point(data = raw_points, aes(x = x, y = response), shape=21, color = "black", fill="black", size=1.5,inherit.aes = FALSE) +
  xlab("Evolution of Stomatal Conductance") + 
  ylab("Predicted Median S") +
  theme_classic()+ theme(
    axis.text.x = element_text(size=16, face="bold"),
    axis.text.y = element_text(size=18,face="bold"),
    axis.title.x = element_text(color="black", size=16, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=18,vjust = 1.7, face="bold",hjust=0.5),
    legend.title = element_blank(),
    legend.position = "none"
  )
pred_B_gg

#Single regression plot
pred_B_plot <- ggplot(filter(pred_df_B, group==0.01), aes(x=x, y=predicted)) + #, colour=group
  geom_ribbon(aes(x=x, ymin=conf.low, ymax=conf.high), fill="grey75") + #, group=group, colour=group
  stat_smooth(method="lm", color="black") +
  xlab("Evolution of Stomatal Conductance") + 
  ylab("Predicted Median S") +
  theme_classic()

# Supplemental Figure Sx
plot_grid(pred_A_gg, pred_B_gg)
ggsave("Graphs/Traits/gas_selection_gg.pdf", width=10, height = 5, units = "in")

#Combined Single regression plot
plot_grid(pred_A_plot, pred_B_plot)
ggsave("Graphs/Traits/Traits_Selection.pdf")

###################################################################################












