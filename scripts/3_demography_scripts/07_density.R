##################################################
# Population density over time

# This script visualizes density trends over time 
# and relates them to genetic metrics
# Last updated 2025 05 16
##################################################

library(tidyverse)
library(RColorBrewer)

### Read in raw count data 
counts <- read_csv("data/demography data/counts_by_plot.csv") %>% 
  filter(Year>2009) %>% #drop years outside the scope of this analysis
  filter(Year<2020) %>% 
  filter(PlotID>17) %>% #2 mistaken observations from plots at Carlon that were discarded prior to 2009
  select(-EffectiveTransectLength) #drop unwanted column

### Derive sampling area from plant coordinates
plot_area <- counts %>% 
  filter(Class!="E") %>% filter(Class!="?") %>% # get rid of areas that could not be sorted out, but keep Class=D to encompass sampled areas without living plants
  group_by(PlotID) %>% 
  mutate(min_x_all = min(xmin),
         max_x_all = max(xmax),
         min_y_all = min(ymin),
         max_y_all = max(ymax),
         max_area_all = (max_x_all - min_x_all)*(max_y_all - min_y_all)) %>% #this is the max extent over which identifiable plants were recorded across all years
  # using this as denominator assumes that a consistent plot area was scanned yearly, so that years with narrower coordinate extents reflect sparser plants and not decreases in sampling extent
  # this will underestimate density in years with exclusion zones and when plots are truncated/expanded in different years
  ungroup() %>% 
  group_by(PlotID, Year) %>% 
  mutate(min_x_year = min(xmin),
         max_x_year = max(xmax),
         min_y_year = min(ymin),
         max_y_year = max(ymax),
         max_area_year = (max_x_year - min_x_year)*(max_y_year - min_y_year)) #this is the max extent over which identifiable plants occurred within a particular year
  # using this as denominator will reflect local density where plants were definitely sampled, allowing sampled area to wax and wane with exclusion zones
  # this makes density inestimable when there is only one plant on a plot, and it will overestimate density when there is a small patch of plants in an otherwise empty plot
ggplot(data=plot_area, aes(x=log(max_area_year), y=log(max_area_all))) + geom_point()

#best solution for exclusion zones would be to calculate their areas and subtract them from sampled area (like swiss cheese)
# but with this data frame we can only calculate one outer bound of all exclusions, more like a donut than swiss cheese
#excluded_area <- counts %>% 
#  filter(Class=="E") %>% 
#  group_by(PlotID, Year) %>% 
#  mutate(min_x_E = min(xmin),
#         max_x_E = max(xmax),
#         min_y_E = min(ymin),
#         max_y_E = max(ymax),
#         max_area_E = (max_x_E - min_x_E)*(max_y_E - min_y_E))

# quick inspect: how different are these?
ggplot(data=plot_area, aes(x=log(max_area_year), y=log(max_area_all))) + geom_point()


density <- plot_area %>% 
  filter(Class=="J"|Class=="A") %>% 
  select(PlotID, Year, Class, Count, max_area_all, max_area_year) %>% 
  pivot_wider(names_from = Class, values_from = Count) %>% 
  mutate_at(c(4:6), ~replace(., is.na(.), 0)) %>% 
  mutate(N = J + A,
         Density_J_all = J/max_area_all,
         Density_A_all = A/max_area_all,
         Density_N_all = N/max_area_all,
         Density_J_year = J/max_area_year,
         Density_A_year = A/max_area_year,
         Density_N_year = N/max_area_year)

density <- density %>% 
  group_by(PlotID) %>% 
  mutate(max_dj_all = max(Density_J_all),
         max_da_all = max(Density_A_all),
         max_dn_all = max(Density_N_all),
         max_dj_year = max(Density_J_year),
         max_da_year = max(Density_A_year),
         max_dn_year = max(Density_N_year)) %>% 
  ungroup() %>% 
  mutate(RelDens_J_all = Density_J_all/max_dj_all,
         RelDens_A_all = Density_A_all/max_da_all,
         RelDens_N_all = Density_N_all/max_dn_all,
         RelDens_J_year = Density_J_year/max_dj_year,
         RelDens_A_year = Density_A_year/max_da_year,
         RelDens_N_year = Density_N_year/max_dn_year) %>% 
  ungroup() #1117 plot-year combos with non-zero density values

### Read in plot census history so that sampled years with no recorded plants can be recoded as density=0
census_hist <- read_csv("data/demography data/PlotCensusHistory.csv") %>% 
  # NA means plot did not exist that year (e.g., not yet created, discarded)
  # n means plot was skipped that year (e.g., site or plot was inaccessible, triaged for time)
  # y means site and plot were visited <-- these are the plot-year combinations that should be 0 density if values are not calculable above
  select(5:15) %>% 
  pivot_longer(-PlotCode, names_to="Year", values_to="Census") #1600 plot*year combos
census_hist$Year <- as.numeric(census_hist$Year)

hist_density <- left_join(census_hist, density, by=c("PlotCode"="PlotID", "Year"="Year")) %>% 
  mutate(
    across(.cols=Density_J_all:RelDens_N_year,
           .fns = ~ ifelse(Census=="y" & is.na(max_area_all), 0, .)
           )
  )

# quick inspect: how different are these?
ggplot(data=hist_density, aes(x=Density_J_all, y=Density_J_year)) + geom_point() + geom_abline(slope=1, intercept=0)
ggplot(data=hist_density, aes(x=Density_A_all, y=Density_A_year)) + geom_point() + geom_abline(slope=1, intercept=0)
ggplot(data=hist_density, aes(x=Density_N_all, y=Density_N_year)) + geom_point() + geom_abline(slope=1, intercept=0)
# Pretty darn different. Year-specific areas generate much higher densities than single max area.
# Though sampling extents do wax and wane a bit from year to year, this is probably less of a factor than plant extents waxing and waning --> use _all instead of _year estimates

### Join with site metadata
counts_meta <- read_csv("data/demography data/tidy_counts.csv") %>% #obsolete density estimates; use file only to link PlotIDs to Sites
  select(PlotID, Site) %>% 
  distinct()
  
density <- left_join(hist_density, counts_meta, by=c("PlotCode"="PlotID")) %>% 
  filter(Site!="MillCreek") %>% #remove unusable sites
  filter(Site!="DeerCreek")


pop_meta <- read_csv("data/genomic_data/pop_meta_data.csv") %>% mutate(Site=gsub(" ", "", Site)) %>% 
  select(1:7)

density <- left_join(density, pop_meta) %>% drop_na(Lat.Color) %>% arrange(Latitude)
unique(density$Site)


### Graphs of density trends over time
color.list = unique(density$Lat.Color)

# density of all plants, calculated with area across all sampled years
ggplot(density, aes(x=Year, y=Density_N_all, color=as.factor(Latitude))) +
  geom_point(alpha=0.3) +
  geom_smooth(method="lm", formula=y~poly(x,2)) +
  ylab("Total plant density") +
  scale_color_manual(values=color.list) +
  facet_wrap(~Latitude, scale="free") +
  theme_classic()

# density of juveniles
ggplot(density, aes(x=Year, y=Density_J_all, color=as.factor(Latitude))) +
  geom_point(alpha=0.3) +
  geom_smooth(method="lm", formula=y~poly(x,2)) +
  ylab("Juvenile plant density") +
  scale_color_manual(values=color.list) +
  facet_wrap(~Latitude, scale="free") +
  theme_classic()

# density of adults
ggplot(density, aes(x=Year, y=Density_A_all, color=as.factor(Latitude))) +
  geom_point(alpha=0.3) +
  geom_smooth(method="lm", formula=y~poly(x,2)) +
  ylab("Adult plant density") +
  scale_color_manual(values=color.list) +
  facet_wrap(~Latitude, scale="free") +
  theme_classic()

##### relative density of all plants <-- USE THIS ONE 
ggplot(density, aes(x=Year, y=RelDens_N_all, color=as.factor(Latitude))) +
  geom_point(alpha=0.3) +
  geom_smooth(method="lm", formula=y~poly(x,2)) +
  ggpubr::stat_regline_equation(formula=y~poly(x,2), label.x=2010, label.y=0.95, size=2.5) + 
  ylim(0,1) +
  ylab("Relative plant density") +
  scale_x_discrete(breaks=c(2010, 2014, 2018), limits=c(2010:2018)) +
  scale_color_manual(values=color.list) +
  facet_wrap(~Latitude, nrow=3) +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text.x = element_blank(),
        legend.title = element_blank())

# relative density of juveniles
ggplot(density, aes(x=Year, y=RelDens_J_all, color=as.factor(Latitude))) +
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
ggplot(density, aes(x=Year, y=RelDens_A_all, color=as.factor(Latitude))) +
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
density_tall <- density %>% select(Year, Site, Lat.Color, Latitude, RelDens_J_all, RelDens_A_all, RelDens_N_all) %>%   pivot_longer(cols=5:7, names_to="Class", values_to="RelDens")

ggplot(density_tall, aes(x=Year, y=RelDens, group=Class, color=as.factor(Latitude))) +
  geom_point(alpha=0.3) +
  geom_smooth(method="glm", method.args = list(family="binomial"), formula=y~poly(x,2)) +
  ylab("Total plant density") +
  scale_color_manual(values=color.list) +
  facet_wrap(~Latitude, scale="free") +
  theme_classic()  +
  theme(strip.background = element_blank(), strip.text.x = element_blank())




### Summarizing indices of density trajectories
site.vec <- unique(density$Site)
slopes.dens.drought <- c()
slopes.dens.recovery <- c()
curves.dens.linear <- c()
curves.dens.quad <- c()

for (i in 1:length(site.vec)) {
  dat <- density %>% filter(Site==site.vec[i])
  mod1 <- lm(RelDens_N_all ~ Year, data=subset(dat, Year<2016))
  mod2 <- lm(RelDens_N_all ~ Year, data=subset(dat, Year>2015))
  mod3 <- lm(RelDens_N_all ~ poly(Year,2), data=dat)
  slopes.dens.drought[i] <- coefficients(mod1)[[2]] #negative coeff = decline
  slopes.dens.recovery[i] <- coefficients(mod2)[[2]] #positive coeff = recovery
  curves.dens.quad[i] <- coefficients(mod3)[[3]] #positive coeff=convex (consistent with rescue) 
  }

slopes.pop <- as.data.frame(cbind(site.vec, slopes.dens.drought, slopes.dens.recovery, curves.dens.quad))

density_trends <- left_join(density, slopes.pop, by=c("Site"="site.vec")) %>% 
  distinct(Site, .keep_all=TRUE) %>% 
  select(Site, slopes.dens.drought, slopes.dens.recovery, curves.dens.quad)

write_csv(density_trends, "data/demography data/density_trends.csv")

### Relationships of density with lambda, selection, and pi

#Import demography data (includes metadata)
demog_recovery <- read_csv("data/demography data/siteYear.lambda_responses_2010-2019.csv") %>% mutate(Site=gsub(" ", "", Site))

# Import pi for demography populations
pi_raw <- read_csv("data/genomic_data/raw_pi_clump.csv") 

#Import selection data
slope.summary <- read_csv("data/snp_change_2/mean_median_S_all.csv") %>% select(Site, Median, Mean)

# Join data sets
pi_pop <- left_join(demog_recovery, pi_raw, by=c("Paper_ID"="Site")) 
geno_pop <- left_join(pi_pop, slope.summary, by=c("Paper_ID"="Site"))
geno_pop <- left_join(geno_pop, density_trends, by="Site")
geno_pop$slopes.dens.drought <- as.numeric(geno_pop$slopes.dens.drought)
geno_pop$slopes.dens.recovery <- as.numeric(geno_pop$slopes.dens.recovery)
geno_pop$curves.dens.quad <- as.numeric(geno_pop$curves.dens.quad)

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

