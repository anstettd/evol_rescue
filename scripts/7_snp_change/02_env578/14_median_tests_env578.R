##################################################################################
## Statistical tests for median SNP slopes from 3 climate variables. Sensitivity analysis.
## Test climate-associated median slope difference from 0
## Test against permuted random slope median
## Done for all 12 populations
## Author Daniel Anstett
## 
## 
## Last Modified Sept 2, 2023
###################################################################################
# Remove objects and clear workspace
rm(list = ls(all=TRUE))

#Function
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
obs_env_unique <- read_csv("data/snp_change_2/slope_obs_all_unique_env578.csv") %>% 
  filter(SE<5) %>% mutate(abs_slope = abs(Slope)) %>% filter(Site!=12)
rand_env_unique <- read_csv("~/Dropbox/z_Documents/aLarge_files/M_gen/rand_slope_histPop_strong_50_50_clump_578.csv") %>% filter(Site!=12)

##Get slope median
median_obs <- obs_env_unique %>% group_by(Site) %>% summarise(median = median(Slope, na.rm = TRUE))
median_rand <- rand_env_unique %>% group_by(Site,Seed_ID) %>% summarise(median = median(Slope, na.rm = TRUE))

#Get slope mean
mean_obs <- obs_env_unique %>% group_by(Site) %>% summarise(mean = mean(Slope, na.rm = TRUE))
mean_rand <- rand_env_unique %>% group_by(Site,Seed_ID) %>% summarise(mean = mean(Slope, na.rm = TRUE))


###################################################################################
##Get SD per site for permuted slopes data
# Calculate standard deviation for each dataset within each population
std_dev_table <- rand_env_unique %>%
  group_by(Site, Seed_ID) %>%
  summarise(standard_deviation = sd(Slope, na.rm = TRUE)) %>%
  ungroup()

# Calculate the range (min and max) of standard deviations for each population
std_dev_range <- std_dev_table %>%
  group_by(Site) %>%
  summarise(
    mean_std_dev = mean(standard_deviation)
  )

#Get min and max
min(std_dev_range[,2])
max(std_dev_range[,2])
###################################################################################

#Test median and mean slope is different from zero
wilcox.out <- as.data.frame(median_obs) 
wilcox.out <- cbind(wilcox.out,mean_obs[,2])

for (i in 1:11){
  obs_test <- as.data.frame(obs_env_unique  %>% filter(Site==i))
  wilcox.out[i,4] <- wilcox.test(obs_test$Slope, mu = 0)$p.value
}

colnames(wilcox.out) <- c("Site","Wilcoxon_Median","Wilcoxon_Mean","Wilcoxon_p_value")

###################################################################################
#Test climate-associated median slope against permuted slope median

emp_out<-as.data.frame(median_obs %>% select(Site,median))



for(i in 1:11){
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

colnames(emp_out) <- c("Site","Median","Median_Percentile","Median_p_value",
                       "Mean","Mean_Percentile","Mean_p_value")

tableS4_578<-left_join(emp_out,wilcox.out, by="Site") %>% select(Site,Median,Wilcoxon_p_value,Median_Percentile,Median_p_value)
write_csv(tableS4_578, "Graphs/snp_change_2/tableS4_578.csv")

###################################################################################
#Make Median S histograms

median_rand$Site <- as.factor(median_rand$Site)
median_rand$Site <- factor(median_rand$Site, levels = c(1,2,3,4,5,6,7,8,9,10,11))

#Median 
histPop <- ggplot(median_rand,aes(x=median))+
  geom_histogram(color="black",fill = "grey70")+
  labs(x = "Response to Selection", y = "Number of Permutations") +
  geom_vline(xintercept=0) +
  theme_ci() + facet_wrap(.~Site) +
geom_vline(data = median_obs, aes(xintercept = median), size=0.5, linetype="dashed",color="red")
histPop 

ggsave("Graphs/snp_change_2/03_rand_median_578_S12.pdf",width=12, height = 8, units = "in")




