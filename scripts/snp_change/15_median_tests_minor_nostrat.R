##################################################################################
## Statistical tests for median slope
## Test climate-associated median slope difference from 0
## Test against permuted random slope median
## Done for all 12 populations
## Author Daniel Anstett
## 
## 
## Last Modified Sept 2, 2023
###################################################################################

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

###################################################################################
#HistPop unique SNP ascorss all env
#Updated for binomial data
obs_env_unique <- read_csv("data/binomial_strong/slope_obs_all_unique.csv") %>% 
  filter(SE<5.5) %>% mutate(abs_slope = abs(Slope))

setwd("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files")
rand_env_unique <- read_csv("rand_slope_histPop_strong_50_50_no_strat.csv")
setwd("~/Dropbox/AM_Workshop/snp_change")

#Get slope median
median_obs <- obs_env_unique %>% group_by(Site) %>% summarise(median = median(Slope, na.rm = TRUE))
median_rand <- rand_env_unique %>% group_by(Site,Seed_ID) %>% summarise(median = median(Slope, na.rm = TRUE))

#Get slope mean
mean_obs <- obs_env_unique %>% group_by(Site) %>% summarise(mean = mean(Slope, na.rm = TRUE))
mean_rand <- rand_env_unique %>% group_by(Site,Seed_ID) %>% summarise(mean = mean(Slope, na.rm = TRUE))


###################################################################################
#Test median and mean slope is different from zero
wilcox.out <- as.data.frame(median_obs) 
wilcox.out <- cbind(wilcox.out,mean_obs[,2])

for (i in 1:12){
  obs_test <- as.data.frame(obs_env_unique  %>% filter(Site==i))
  wilcox.out[i,4] <- wilcox.test(obs_test$Slope, mu = 0)$p.value
}

colnames(wilcox.out) <- c("Site","Median","Mean","p-value")


#write_csv(wilcox.out, "data/binomial_strong/wilcox_S_minor_nostrat.csv")


###################################################################################
#Test climate-associated median slope against permuted slope median

emp_out<-as.data.frame(median_obs %>% select(Site,median))



for(i in 1:12){
median_rand_1 <- as.data.frame(median_rand %>% filter(Site==i))
median_rand_vec <- as.vector(median_rand_1[,3])
mean_rand_1 <- as.data.frame(mean_rand %>% filter(Site==i))
mean_rand_vec <- as.vector(mean_rand_1[,3])

#percentile <- pull(median_obs[1,2]) %>% percent_rank(median_rand_vec)
emp_out[i,3] <- ecdf(median_rand_vec)(pull(median_obs[i,2]))
emp_out[i,4] <- 1 - emp_out[i,3]
emp_out[i,5] <- mean_obs[i,2]
emp_out[i,6] <- ecdf(mean_rand_vec)(pull(mean_obs[i,2]))
emp_out[i,7] <- 1 - emp_out[i,6]
#emp_out[i,3] <- sum(median_rand_vec >= pull(median_obs[i,2])) / 1000

}

colnames(emp_out) <- c("Site","Median","Median_Percentile","Median_p-value",
                       "Mean","Mean_Percentile","Mean_p-value")

#write_csv(emp_out, "data/binomial_strong/mean_median_S_minor_nostrat.csv")

###################################################################################
#Make Median and Mean histograms

median_rand$Site <- as.factor(median_rand$Site)
median_rand$Site <- factor(median_rand$Site, levels = c(1,12,2,3,4,5,6,7,8,9,10,11))

mean_rand$Site <- as.factor(mean_rand$Site)
mean_rand$Site <- factor(mean_rand$Site, levels = c(1,12,2,3,4,5,6,7,8,9,10,11))


#Median 
histPop <- ggplot(median_rand,aes(x=median))+
  geom_histogram(color="black",fill = "pink")+
  labs(x = "Non-Climate Associated S Median", y = "Number of Permutations") +
  geom_vline(xintercept=0) +
  theme_ci() + facet_wrap(.~Site) +
geom_vline(data = median_obs, aes(xintercept = median), linetype="dashed",color="red")
histPop 

ggsave("graphs/mean_median_s/50_50/major_minor_median_nostrat.pdf",width=12, height = 8, units = "in")



#Mean
histPop_mean <- ggplot(mean_rand,aes(x=mean))+
  geom_histogram(color="black",fill = "pink")+
  labs(x = "Non-Climate Associated S Mean", y = "Number of Permutations") +
  geom_vline(xintercept=0) +
  theme_ci() + facet_wrap(.~Site) +
  geom_vline(data = mean_obs, aes(xintercept = mean), linetype="dashed",color="red")
histPop_mean

ggsave("graphs/mean_median_s/50_50/major_minor_mean_nostrat.pdf",width=12, height = 8, units = "in")


###################################################################################
#Site specific histograms

median_rand_pops <- median_rand %>% filter(Site == 2 | Site == 3 | Site == 11)
median_obs_pops <- median_obs %>% filter(Site == 2 | Site == 3 | Site == 11)

median_rand_pops$Site <- droplevels(median_rand_pops$Site) %>% na.omit()

#Median 
histPops <- ggplot(median_rand_pops,aes(x=median))+
  geom_histogram(color="black",fill = "pink")+
  labs(x = "Non-Climate Associated S Median", y = "Number of Permutations") +
  geom_vline(xintercept=0) +
  theme_ci() + facet_wrap(.~Site) +
  geom_vline(data = median_obs_pops, aes(xintercept = median), linetype="dashed",color="red")+
  theme(strip.text.x = element_text(size=0))
histPops

ggsave("graphs/mean_median_s/50_50/03_median_nostrat_pop.pdf",width=11, height = 3.5, units = "in")



median_rand_pops <- median_rand %>% filter(Site == 3 | Site == 11)
median_obs_pops <- median_obs %>% filter(Site == 3 | Site == 11)

median_rand_pops$Site <- droplevels(median_rand_pops$Site) %>% na.omit()

#Median 
histPops <- ggplot(median_rand_pops,aes(x=median))+
  geom_histogram(color="black",fill = "pink")+
  labs(x = "Median SNP Change", y = "Number of Permutations") +
  geom_vline(xintercept=0) +
  theme_ci() + facet_wrap(.~Site) +
  geom_vline(data = median_obs_pops, aes(xintercept = median), linetype="dashed",color="red")+
  scale_x_continuous(breaks=c(-0.05,0,0.05))+
  theme(strip.text.x = element_text(size=0))
histPops

ggsave("graphs/mean_median_s/50_50/04_median_nostrat_pop.pdf",width=9, height = 3.5, units = "in")




