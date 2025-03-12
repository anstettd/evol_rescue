##################################################
# Population density over time

# This script visualizes density trends over time
# Last updated 2025 03 11
##################################################

library(tidyverse)
library(RColorBrewer)

### Read in tidy data & remove Mill Creek and Deer Creek
counts <- read_csv("data/demography data/tidy_counts.csv")

pop_meta <- read_csv("data/genomic_data/pop_meta_data.csv") %>% mutate(Site=gsub(" ", "", Site)) %>% select(Site, Lat.Color)

counts <- left_join(counts, pop_meta) %>% drop_na(Lat.Color) %>% arrange(SiteLatitude)
unique(counts$Site)

### Graphs of density trends over time
color.list = unique(counts$Lat.Color)

# log density
ggplot(counts, aes(x=Year, y=logdens, color=as.factor(SiteLatitude))) +
  geom_point(alpha=0.3) +
  geom_smooth(method="lm", formula=y~poly(x,2)) +
  ylab("Log(density)") +
  scale_color_manual(values=color.list) +
  facet_wrap(~SiteLatitude, scale="free") +
  theme_classic()


# relative density
ggplot(counts, aes(x=Year, y=reldens, color=as.factor(SiteLatitude))) +
  geom_point(alpha=0.3) +
  geom_smooth(method="lm", formula=y~poly(x,2)) +
  ggpubr::stat_regline_equation(formula=y~poly(x,2), label.x=2010, label.y=0.95, size=2.5) + 
  ylim(0,1) +
  ylab("Relative density") +
  scale_color_manual(values=color.list) +
  facet_wrap(~SiteLatitude, nrow=3, scale="free") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        legend.title = element_blank())


### Summarizing indices of density trajectories
site.vec <- unique(counts$SiteID)
slopes.dens.drought <- c()
slopes.dens.recovery <- c()
curves.dens.linear <- c()
curves.dens.quad <- c()

for (i in 1:length(site.vec)) {
  dat <- counts %>% filter(SiteID==site.vec[i])
  mod1 <- lm(reldens ~ Year, data=subset(dat, Year<2016))
  mod2 <- lm(reldens ~ Year, data=subset(dat, Year>2015))
  mod3 <- lm(reldens ~ poly(Year,2), data=dat)
  slopes.dens.drought[i] <- coefficients(mod1)[[2]] #negative coeff = decline
  slopes.dens.recovery[i] <- coefficients(mod2)[[2]] #positive coeff = recovery
  curves.dens.quad[i] <- coefficients(mod3)[[3]] #positive coeff=convex (consistent with rescue) 
  }

slopes.pop <- as.data.frame(cbind(site.vec, slopes.dens.drought, slopes.dens.recovery, curves.dens.quad))

density_trends <- left_join(counts, slopes.pop, by=c("SiteID"="site.vec")) %>% 
  distinct(SiteID, .keep_all=TRUE) %>% 
  select(Site, SiteID, slopes.dens.drought, slopes.dens.recovery, curves.dens.quad)

write_csv(density_trends, "data/demography data/density_trends.csv")

### Relationships of density with lambda, selection, and pi

#Import demography data (includes metadata)
demog_recovery <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv") %>% mutate(Site=gsub(" ", "", Site))

# Import pi for demography populations
pi_raw <- read_csv("data/genomic_data/raw_pi.csv") 

#Import selection data
slope.summary <- read_csv("data/snp_change_data/mean_median_S.csv") %>% select(Site, Median, Mean)

# Join data sets
pi_pop <- left_join(demog_recovery, pi_raw, by=c("Paper_ID"="Site")) 
geno_pop <- left_join(pi_pop, slope.summary, by=c("Paper_ID"="Site"))
geno_pop <- left_join(geno_pop, density_trends, by="Site")

## Visualize scatterplots

color.list = unique(geno_pop$Lat.Color)

# density ~ latitude
ggplot(geno_pop, aes(x=Latitude, y=slopes.dens.drought)) + 
  geom_point() + 
  geom_smooth(method="lm", se=TRUE, color="black") + 
  #ggpubr::stat_regline_equation(formula=y~x) + 
  ggpubr::stat_cor(method = "pearson") + 
  theme_classic()
ggplot(geno_pop, aes(x=Latitude, y=slopes.dens.recovery)) + 
  geom_point() + 
  scale_color_manual(values=color.list) + 
  geom_smooth(method="lm", se=TRUE, color="black") + 
  #ggpubr::stat_regline_equation(formula=y~x) + 
  ggpubr::stat_cor(method = "pearson") + 
  theme_classic()
ggplot(geno_pop, aes(x=Latitude, y=curves.dens.quad)) + 
  geom_point() + 
  geom_smooth(method="lm", se=TRUE, color="black") + 
  #ggpubr::stat_regline_equation(formula=y~x) + 
  ggpubr::stat_cor(method = "pearson") + 
  theme_classic()

# density ~ climate snp diversity
ggplot(geno_pop, aes(x=pi_snp_set, y=slopes.dens.drought)) + 
  geom_point() + 
  geom_smooth(method="lm", se=TRUE, color="black") + 
  #ggpubr::stat_regline_equation(formula=y~x) + 
  ggpubr::stat_cor(method = "pearson") + 
  theme_classic()
ggplot(geno_pop, aes(x=pi_snp_set, y=slopes.dens.recovery)) + 
  geom_point() + 
  geom_smooth(method="lm", se=TRUE, color="black") + 
  #ggpubr::stat_regline_equation(formula=y~x) + 
  ggpubr::stat_cor(method = "pearson") + 
  theme_classic()
ggplot(geno_pop, aes(x=pi_snp_set, y=curves.dens.quad)) + 
  geom_point() + 
  geom_smooth(method="lm", se=TRUE, color="black") + 
  #ggpubr::stat_regline_equation(formula=y~x) + 
  ggpubr::stat_cor(method = "pearson") + 
  theme_classic()

# density ~ all snp diversity
ggplot(geno_pop, aes(x=pi_all_snps, y=slopes.dens.drought)) + 
  geom_point() + 
  geom_smooth(method="lm", se=TRUE, color="black") + 
  #ggpubr::stat_regline_equation(formula=y~x) + 
  ggpubr::stat_cor(method = "pearson") + 
  theme_classic()
ggplot(geno_pop, aes(x=pi_all_snps, y=slopes.dens.recovery)) + 
  geom_point() + 
  geom_smooth(method="lm", se=TRUE, color="black") + 
  #ggpubr::stat_regline_equation(formula=y~x) + 
  ggpubr::stat_cor(method = "pearson") + 
  theme_classic()
ggplot(geno_pop, aes(x=pi_all_snps, y=curves.dens.quad)) + 
  geom_point() + 
  geom_smooth(method="lm", se=TRUE, color="black") + 
  #ggpubr::stat_regline_equation(formula=y~x) + 
  ggpubr::stat_cor(method = "pearson") + 
  theme_classic()

# density ~ selection
ggplot(geno_pop, aes(x=Median, y=slopes.dens.drought)) + 
  geom_point() + 
  geom_smooth(method="lm", se=TRUE, color="black") + 
  #ggpubr::stat_regline_equation(formula=y~x) + 
  ggpubr::stat_cor(method = "pearson") + 
  theme_classic()
ggplot(geno_pop, aes(x=Median, y=slopes.dens.recovery)) + 
  geom_point() + 
  geom_smooth(method="lm", se=TRUE, color="black") + 
  #ggpubr::stat_regline_equation(formula=y~x) + 
  ggpubr::stat_cor(method = "pearson") + 
  theme_classic()
ggplot(geno_pop, aes(x=Median, y=curves.dens.quad)) + 
  geom_point() + 
  geom_smooth(method="lm", se=TRUE, color="black") + 
  #ggpubr::stat_regline_equation(formula=y~x) + 
  ggpubr::stat_cor(method = "pearson") + 
  theme_classic()

# selection ~ diversity
ggplot(geno_pop, aes(x=pi_snp_set, y=Median)) + 
  geom_point() + 
  geom_smooth(method="lm", se=TRUE, color="black") + 
  #ggpubr::stat_regline_equation(formula=y~x) + 
  ggpubr::stat_cor(method = "pearson") + 
  theme_classic()

ggplot(geno_pop, aes(x=pi_all_snps, y=Median)) + 
  geom_point() + 
  geom_smooth(method="lm", se=TRUE, color="black") + 
  #ggpubr::stat_regline_equation(formula=y~x) + 
  ggpubr::stat_cor(method = "pearson") + 
  theme_classic()

