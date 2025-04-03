##################################################
# Population density over time

# This script visualizes density trends over time 
# and relates them to genetic metrics
# Last updated 2025 04 03
##################################################

library(tidyverse)
library(RColorBrewer)

### Read in raw count data 
counts <- read_csv("data/demography data/counts_by_plot.csv") %>% 
  filter(Year>2009) %>% #drop years outside the scope of this analysis
  filter(Year<2020) %>% 
  select(-EffectiveTransectLength) #drop unwanted column

### Derive sampled area across all years (J + A + D; exclude E)
plot_area <- counts %>% 
  group_by(PlotID) %>% 
  filter(Class!="E") %>% 
  mutate(min_x = min(xmin),
         max_x = max(xmax),
         min_y = min(ymin),
         max_y = max(ymax),
         max_area = (max_x - min_x)*(max_y - min_y))

density <- plot_area %>% 
  filter(Class=="J"|Class=="A") %>% 
  select(PlotID, Year, Class, Count, max_area) %>% 
  pivot_wider(names_from = Class, values_from = Count) %>% 
  mutate_at(c(4:5), ~replace(., is.na(.), 0)) %>% 
  mutate(N = J + A,
         Density_J = J/max_area,
         Density_A = A/max_area,
         Density_N = N/max_area)

density <- density %>% 
  group_by(PlotID) %>% 
  mutate(max_dj = max(Density_J),
         max_da = max(Density_A),
         max_dn = max(Density_N)) %>% 
  ungroup() %>% 
  mutate(RelDens_J = Density_J/max_dj,
         RelDens_A = Density_A/max_da,
         RelDens_N = Density_N/max_dn)

### Join with site metadata
counts_meta <- read_csv("data/demography data/tidy_counts.csv") %>% #obsolete density estimates; use file only to link PlotIDs to Sites
  select(PlotID, Site) %>% 
  distinct()
  
density <- left_join(density, counts_meta) %>% 
  filter(Site!="MillCreek") %>% #remove unusable sites
  filter(Site!="DeerCreek")


pop_meta <- read_csv("data/genomic_data/pop_meta_data.csv") %>% mutate(Site=gsub(" ", "", Site)) %>% 
  select(1:7)

density <- left_join(density, pop_meta) %>% drop_na(Lat.Color) %>% arrange(Latitude)
unique(density$Site)


### Graphs of density trends over time
color.list = unique(density$Lat.Color)

# density of all plants
ggplot(density, aes(x=Year, y=Density_N, color=as.factor(Latitude))) +
  geom_point(alpha=0.3) +
  geom_smooth(method="lm", formula=y~poly(x,2)) +
  ylab("Total plant density") +
  scale_color_manual(values=color.list) +
  facet_wrap(~Latitude, scale="free") +
  theme_classic()

# density of juveniles
ggplot(density, aes(x=Year, y=Density_J, color=as.factor(Latitude))) +
  geom_point(alpha=0.3) +
  geom_smooth(method="lm", formula=y~poly(x,2)) +
  ylab("Juvenile plant density") +
  scale_color_manual(values=color.list) +
  facet_wrap(~Latitude, scale="free") +
  theme_classic()

# density of adults
ggplot(density, aes(x=Year, y=Density_A, color=as.factor(Latitude))) +
  geom_point(alpha=0.3) +
  geom_smooth(method="lm", formula=y~poly(x,2)) +
  ylab("Adult plant density") +
  scale_color_manual(values=color.list) +
  facet_wrap(~Latitude, scale="free") +
  theme_classic()

# relative density of all plants
ggplot(density, aes(x=Year, y=RelDens_N, color=as.factor(Latitude))) +
  geom_point(alpha=0.3) +
  geom_smooth(method="lm", formula=y~poly(x,2)) +
  ggpubr::stat_regline_equation(formula=y~poly(x,2), label.x=2010, label.y=0.95, size=2.5) + 
  ylim(0,1) +
  ylab("Relative plant density") +
  scale_color_manual(values=color.list) +
  facet_wrap(~Latitude, nrow=3, scale="free") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        legend.title = element_blank())

# relative density of juveniles
ggplot(density, aes(x=Year, y=RelDens_J, color=as.factor(Latitude))) +
  geom_point(alpha=0.3) +
  geom_smooth(method="glm", method.args = list(family="binomial"), formula=y~poly(x,2)) +
  ggpubr::stat_regline_equation(formula=y~poly(x,2), label.x=2010, label.y=0.95, size=2.5) + 
  ylim(0,1) +
  ylab("Relative density of juvenile plants") +
  scale_color_manual(values=color.list) +
  facet_wrap(~Latitude, nrow=3, scale="free") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        legend.title = element_blank())

# relative density of adults
ggplot(density, aes(x=Year, y=RelDens_A, color=as.factor(Latitude))) +
  geom_point(alpha=0.3) +
  geom_smooth(method="glm", method.args = list(family="binomial"), formula=y~poly(x,2)) +
  ggpubr::stat_regline_equation(formula=y~poly(x,2), label.x=2010, label.y=0.95, size=2.5) + 
  ylim(0,1) +
  ylab("Relative density of adult plants") +
  scale_color_manual(values=color.list) +
  facet_wrap(~Latitude, nrow=3, scale="free") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        legend.title = element_blank())

# classes overlain
density_tall <- density %>% select(Year, Site, Lat.Color, Latitude, RelDens_J, RelDens_A, RelDens_N) %>%   pivot_longer(cols=5:7, names_to="Class", values_to="RelDens")

ggplot(density_tall, aes(x=Year, y=RelDens, group=Class, color=as.factor(Latitude))) +
  geom_point(alpha=0.3) +
  geom_smooth(method="glm", method.args = list(family="binomial"), formula=y~poly(x,2)) +
  ylab("Total plant density") +
  scale_color_manual(values=color.list) +
  facet_wrap(~Latitude, scale="free") +
  theme_classic()  +
  theme(strip.background = element_blank(), strip.text.x = element_blank())




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

# density ~ lambda
ggplot(geno_pop, aes(x=mean.lambda.drought, y=slopes.dens.drought)) + 
  geom_point() + 
  geom_smooth(method="lm", se=TRUE, color="black") + 
  #ggpubr::stat_regline_equation(formula=y~x) + 
  ggpubr::stat_cor(method = "pearson") + 
  theme_classic()
ggplot(geno_pop, aes(x=mean.lambda.recovery, y=slopes.dens.recovery)) + 
  geom_point() + 
  scale_color_manual(values=color.list) + 
  geom_smooth(method="lm", se=TRUE, color="black") + 
  #ggpubr::stat_regline_equation(formula=y~x) + 
  ggpubr::stat_cor(method = "pearson") + 
  theme_classic()
ggplot(geno_pop, aes(x=mean.lambda.drought, y=curves.dens.quad)) + 
  geom_point() + 
  geom_smooth(method="lm", se=TRUE, color="black") + 
  #ggpubr::stat_regline_equation(formula=y~x) + 
  ggpubr::stat_cor(method = "pearson") + 
  theme_classic()
ggplot(geno_pop, aes(x=mean.lambda.recovery, y=curves.dens.quad)) + 
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

