##################################################################################
## Get SNP abundances for baseline and timeseries 
## Done for current snp set (New WZA BF>10 + all BF > 30)
## Author Daniel Anstett
## 
## Done for all 9 env variables
## Last Modified April 24, 2023
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)

#Import strong SNP set filtered for Baseline nonlinearity
snp_set_env1  <- read_csv("data/binomial_strong/strong_snp_set_env1.csv")
snp_set_env2  <- read_csv("data/binomial_strong/strong_snp_set_env2.csv")
snp_set_env3  <- read_csv("data/binomial_strong/strong_snp_set_env3.csv")
snp_set_env4  <- read_csv("data/binomial_strong/strong_snp_set_env4.csv")
snp_set_env5  <- read_csv("data/binomial_strong/strong_snp_set_env5.csv")
snp_set_env6  <- read_csv("data/binomial_strong/strong_snp_set_env6.csv")
snp_set_env7  <- read_csv("data/binomial_strong/strong_snp_set_env7.csv")
snp_set_env8  <- read_csv("data/binomial_strong/strong_snp_set_env8.csv")
snp_set_env9  <- read_csv("data/binomial_strong/strong_snp_set_env9.csv")


#Import full snp table for timeseries
pop_order_time<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")
snp_time<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", header=F, sep=" ")
loci_time<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/timeseries_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", header=F, sep="\t")
colnames(loci_time) <- c("Chromosome","SNP")
loci_united_time <- loci_time %>% unite(chr_snp,"Chromosome","SNP",sep="_")
loci_snp_time <-cbind(loci_united_time,snp_time) #add snp lables to rows


###################################################################################
#Import Baseline data
#import full snp table
pop_order_base<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", 
                           header=F, sep="\t")
snp_base<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", 
                     header=F, sep=" ")
loci_base<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", 
                      header=F, sep="\t")
colnames(loci_base) <- c("Chromosome","SNP")
loci_united_base <- loci_base %>% unite(chr_snp,"Chromosome","SNP",sep="_")
loci_snp_base <-cbind(loci_united_base,snp_base) #add snp lables to rows


#Filter timeseries by snp set
snp_set_time_env1 <-loci_snp_time %>% filter (chr_snp %in% as.character(snp_set_env1$chr_snp))
snp_set_time_env2 <-loci_snp_time %>% filter (chr_snp %in% as.character(snp_set_env2$chr_snp))
snp_set_time_env3 <-loci_snp_time %>% filter (chr_snp %in% as.character(snp_set_env3$chr_snp))
snp_set_time_env4 <-loci_snp_time %>% filter (chr_snp %in% as.character(snp_set_env4$chr_snp))
snp_set_time_env5 <-loci_snp_time %>% filter (chr_snp %in% as.character(snp_set_env5$chr_snp))
snp_set_time_env6 <-loci_snp_time %>% filter (chr_snp %in% as.character(snp_set_env6$chr_snp))
snp_set_time_env7 <-loci_snp_time %>% filter (chr_snp %in% as.character(snp_set_env7$chr_snp))
snp_set_time_env8 <-loci_snp_time %>% filter (chr_snp %in% as.character(snp_set_env8$chr_snp))
snp_set_time_env9 <-loci_snp_time %>% filter (chr_snp %in% as.character(snp_set_env9$chr_snp))

#Filter baseline by snp set
snp_set_base_env1 <-loci_snp_base %>% filter (chr_snp %in% as.character(snp_set_env1$chr_snp))
snp_set_base_env2 <-loci_snp_base %>% filter (chr_snp %in% as.character(snp_set_env2$chr_snp))
snp_set_base_env3 <-loci_snp_base %>% filter (chr_snp %in% as.character(snp_set_env3$chr_snp))
snp_set_base_env4 <-loci_snp_base %>% filter (chr_snp %in% as.character(snp_set_env4$chr_snp))
snp_set_base_env5 <-loci_snp_base %>% filter (chr_snp %in% as.character(snp_set_env5$chr_snp))
snp_set_base_env6 <-loci_snp_base %>% filter (chr_snp %in% as.character(snp_set_env6$chr_snp))
snp_set_base_env7 <-loci_snp_base %>% filter (chr_snp %in% as.character(snp_set_env7$chr_snp))
snp_set_base_env8 <-loci_snp_base %>% filter (chr_snp %in% as.character(snp_set_env8$chr_snp))
snp_set_base_env9 <-loci_snp_base %>% filter (chr_snp %in% as.character(snp_set_env9$chr_snp))


write_csv(snp_set_time_env1,"data/binomial_strong/snp_set_time_env1.csv")
write_csv(snp_set_time_env2,"data/binomial_strong/snp_set_time_env2.csv")
write_csv(snp_set_time_env3,"data/binomial_strong/snp_set_time_env3.csv")
write_csv(snp_set_time_env4,"data/binomial_strong/snp_set_time_env4.csv")
write_csv(snp_set_time_env5,"data/binomial_strong/snp_set_time_env5.csv")
write_csv(snp_set_time_env6,"data/binomial_strong/snp_set_time_env6.csv")
write_csv(snp_set_time_env7,"data/binomial_strong/snp_set_time_env7.csv")
write_csv(snp_set_time_env8,"data/binomial_strong/snp_set_time_env8.csv")
write_csv(snp_set_time_env9,"data/binomial_strong/snp_set_time_env9.csv")


write_csv(snp_set_base_env1,"data/binomial_strong/snp_set_base_env1.csv")
write_csv(snp_set_base_env2,"data/binomial_strong/snp_set_base_env2.csv")
write_csv(snp_set_base_env3,"data/binomial_strong/snp_set_base_env3.csv")
write_csv(snp_set_base_env4,"data/binomial_strong/snp_set_base_env4.csv")
write_csv(snp_set_base_env5,"data/binomial_strong/snp_set_base_env5.csv")
write_csv(snp_set_base_env6,"data/binomial_strong/snp_set_base_env6.csv")
write_csv(snp_set_base_env7,"data/binomial_strong/snp_set_base_env7.csv")
write_csv(snp_set_base_env8,"data/binomial_strong/snp_set_base_env8.csv")
write_csv(snp_set_base_env9,"data/binomial_strong/snp_set_base_env9.csv")









