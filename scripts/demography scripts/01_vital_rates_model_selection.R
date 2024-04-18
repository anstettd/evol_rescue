#### PROJECT: Genomic offsets and demographic trajectories of Mimulus cardinalis populations during extreme drought
#### PURPOSE OF THIS SCRIPT: Perform model selection for each vital rate for subsequent use in IPMs
#### AUTHOR: Seema Sheth and Amy Angert
#### DATE LAST MODIFIED: 20230619

#*******************************************************************************
#### 0. Clean workspace and load required packages
#*******************************************************************************

# Remove objects and clear workspace
rm(list = ls(all=TRUE))

# Make vector of packages needed
packages_needed <- c("lme4", "MuMIn", "MASS", "pscl", "tidyverse", "glmmTMB", "R2admb") 

# Install packages needed (if not already installed)
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
}

# if there are errors related to glmmTMB, try: install.packages("glmmTMB", type="source")

# Load packages needed
for (i in 1:length(packages_needed)){
  library(packages_needed[i], character.only = TRUE)
}


#*******************************************************************************
#### 2. Read in vital rate data ###
#*******************************************************************************

data <- read.csv("data/demography data/Mcard_demog_data_2010-2019_cleanindivs.csv")
data$Site = factor(data$Site)
data$Year = factor(data$Year)

#*******************************************************************************
#### 2. Survival models ###
#*******************************************************************************

# Fixed effects model w/ and w/out size
s1 <- glm(Surv ~ logSize, data=data, family=binomial)
s2 <- glm(Surv ~ 1, data=data, family=binomial)
model.sel(s1, s2) # model w/ size is preferred

# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (nested within Year)
s3 <- glmer(Surv  ~logSize + (logSize|Year/Site), data=data, family=binomial, control=glmerControl(optimizer = "bobyqa")) 
# NOTE: this model has a singularity warning

# Random intercepts & constant slopes for Year
# Random intercepts & constant slopes for Site (nested within Year)
s4 <- glmer(Surv ~ logSize  +(1|Year/Site), data=data, family=binomial, control=glmerControl(optimizer = "bobyqa")) 

# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
s5 <- glmer(Surv ~ logSize + (logSize|Year) + (logSize|Site), data=data, family=binomial, control=glmerControl(optimizer = "bobyqa")) 

# Random intercepts & random slopes for Year
# Random intercepts & constant slopes for Site (not nested within Year)
s6 <- glmer(Surv ~ logSize + (logSize|Year) + (1|Site), data=data, family=binomial, control=glmerControl(optimizer = "bobyqa")) 

# Random intercepts & constant slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
s7 <- glmer(Surv ~ logSize + (1|Year) + (logSize|Site), data=data, family=binomial, control=glmerControl(optimizer = "bobyqa")) 

# Compare models
anova(s3, s4, s5, s6, s7)
model.sel(s1, s2, s3, s4, s5, s6, s7) 

# PREFERRED MODEL IS s3, but due to singularity issues, we are going with the next best model that doesn't have singularity issues, s4
r.squaredGLMM(s4) 

# Save top survival model to .rda file
save(s4, file='data/demography data/surv.reg.rda')   

#*******************************************************************************
#### 4. Growth ###
#*******************************************************************************

# Fixed effects model w/ and w/out size
g1 <- glm(logSizeNext ~ logSize, data=data[!is.na(data$logSize),])
g2 <- glm(logSizeNext ~ 1, data=data[!is.na(data$logSize),])
model.sel(g1,g2) # model w/ size is preferred

# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (nested within Year)
g3 <- lmer(logSizeNext ~ logSize + (logSize|Year/Site), data=data, control=lmerControl(optimizer = "bobyqa")) 

# Random intercepts & constant slopes for Year
# Random intercepts & constant slopes for Site (nested within Year)
g4 <- lmer(logSizeNext ~ logSize + (1|Year/Site), data=data, control=lmerControl(optimizer = "bobyqa")) 

# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
g5 <- lmer(logSizeNext ~ logSize + (logSize|Year) + (logSize|Site), data=data, control=lmerControl(optimizer = "bobyqa")) 

# Random intercepts & random slopes for Year
# Random intercepts & constant slopes for Site (not nested within Year)
g6 <- lmer(logSizeNext ~ logSize + (logSize|Year) + (1|Site), data=data, control=lmerControl(optimizer = "bobyqa")) 

# Random intercepts & constant slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
g7 <- lmer(logSizeNext ~ logSize + (1|Year) + (logSize|Site), data=data, control=lmerControl(optimizer = "bobyqa")) 

# Compare models
anova(g3, g4, g5, g6, g7)
model.sel(g1, g2, g3, g4, g5, g6, g7)

# # PREFERRED MODEL IS g3 and this has no singularity issues after addition of 2017-19 data
r.squaredGLMM(g3) 

# Save top growth model to .rda file
save(g3, file='data/demography data/growth.reg.rda')   

#*******************************************************************************
#### 5. Flowering ###
#*******************************************************************************

# fixed effects model w/ and w/out size
fl1 <- glm(Fec0 ~ logSize, data=data, family=binomial)
fl2 <- glm(Fec0 ~ 1, data=data, family=binomial)
model.sel(fl1, fl2) # model w/ size is preferred

# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (nested within Year)
fl3 <- glmer(Fec0 ~ logSize + (logSize|Year/Site), data=data, family=binomial, control=glmerControl(optimizer = "bobyqa")) 
# NOTE: singularity warning

# Random intercepts & constant slopes for Year
# Random intercepts & constant slopes for Site (nested within Year)
fl4 <- glmer(Fec0 ~ logSize + (1|Year/Site), data=data, family=binomial, control=glmerControl(optimizer = "bobyqa")) 

# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
fl5 <- glmer(Fec0 ~ logSize + (logSize|Year) + (logSize|Site), data=data, family=binomial, control=glmerControl(optimizer = "bobyqa")) 

# Random intercepts & random slopes for Year
# Random intercepts & constant slopes for Site (not nested within Year)
fl6 <- glmer(Fec0 ~ logSize + (logSize|Year) + (1|Site), data=data, family=binomial, control=glmerControl(optimizer = "bobyqa")) 

# Random intercepts & constant slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
fl7 <- glmer(Fec0 ~ logSize + (1|Year) + (logSize|Site), data=data, family=binomial, control=glmerControl(optimizer = "bobyqa")) 

# Compare models
anova(fl3, fl4, fl5, fl6, fl7)
model.sel(fl1, fl2, fl3, fl4, fl5, fl6, fl7) 

# # PREFERRED MODEL IS fl3, but due to singularity issues, we are going with the next best model, fl4
r.squaredGLMM(fl4) 

# Save top flowering model to .rda file
save(fl4, file='data/demography data/flowering.reg.rda')   

#*******************************************************************************
#### 6. Fruit number ###
#*******************************************************************************
  
##6A. Fit fixed effects models (GLMs) only for initial exploratory model selection of variance structure (poisson vs. negative binomial vs. 0- inflation)
fr1_glm <- glm(Fec1 ~ logSize, data=data, na.action=na.omit, family=poisson) # poisson without 0-inflation 
fr2_glm <- glm.nb(Fec1 ~ logSize, data=data, na.action=na.omit) # negative binomial without 0-inflation
fr3_glm <- zeroinfl(Fec1 ~ logSize, data=data, na.action=na.omit, dist="poisson") # poisson with 0-inflation
fr4_glm <- zeroinfl(Fec1 ~ logSize, data=data, na.action=na.omit, dist="negbin") # negative binomial with 0-inflation

model.sel(fr1_glm, fr2_glm, fr3_glm, fr4_glm) 
	   
# PREFERRED MODEL IS fr2 (negative binomial w/out 0-inflation)

# Fixed effects model w/ and w/out size
fr5_glm <- glm.nb(Fec1 ~ 1, data=data, na.action=na.omit)
model.sel(fr2_glm, fr5_glm) # model w/ size is preferred
	  	
##6B. Model selection of random effects structure

# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (nested within Year)
fr3 <- glmmTMB(Fec1 ~ logSize + (logSize|Year/Site), data=data[!is.na(data$Fec1),], family=nbinom1()) 
# NOTE: singularity warning

# Random intercepts & constant slopes for Year
# Random intercepts & constant slopes for Site (nested within Year)
fr4 <- glmmTMB(Fec1 ~ logSize + (1|Year/Site), data=data[!is.na(data$Fec1),], family=nbinom1()) 
      
# Random intercepts & random slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
fr5 <- glmmTMB(Fec1 ~ logSize + (logSize|Year) + (logSize|Site), data=data[!is.na(data$Fec1),], family=nbinom1) 

# Random intercepts & random slopes for Year
# Random intercepts & constant slopes for Site (not nested within Year)
fr6 <- glmmTMB(Fec1 ~ logSize + (logSize|Year) + (1|Site), data=data[!is.na(data$Fec1),], family=nbinom1()) 

# Random intercepts & constant slopes for Year
# Random intercepts & random slopes for Site (not nested within Year)
fr7 <- glmmTMB(Fec1 ~ logSize + (1|Year) + (logSize|Site), data=data[!is.na(data$Fec1),], family=nbinom1) 

# No random effects
fr8 <- glmmTMB(Fec1 ~ logSize, data=data[!is.na(data$Fec1),], family=nbinom1) 

# Intercept only
fr9 <- glmmTMB(Fec1 ~ 1, data=data[!is.na(data$Fec1),], family=nbinom1) 

# Compare models
anova(fr3, fr4, fr5, fr6, fr7)
model.sel(fr3, fr4, fr5, fr6, fr7, fr8, fr9) 

# PREFERRED MODEL IS fr3, but due to singularity issues, we are going with the next best model, fr5
r.squaredGLMM(fr5)

# Save top fruit # model to .rda file 
save(fr5, file='data/demography data/fruit.reg.rda')   
