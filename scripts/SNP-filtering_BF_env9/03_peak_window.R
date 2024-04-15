#############################################################################################################
## Get peak windows from manhattan plots with empirical p-values
## Author Daniel Anstett
## 
## 
## 
## Last Modified Sept 2, 2022
#############################################################################################################
#Import libraries

library(tidyverse)
library(Kendall)
library(plotly)

#Import files

#Windows
wza_win_env1 <- read_csv("data/genomic_data/WZA_win_env1_bf.csv")
wza_win_env2 <- read_csv("data/genomic_data/WZA_win_env2_bf.csv")
wza_win_env3 <- read_csv("data/genomic_data/WZA_win_env3_bf.csv")
wza_win_env4 <- read_csv("data/genomic_data/WZA_win_env4_bf.csv")
wza_win_env5 <- read_csv("data/genomic_data/WZA_win_env5_bf.csv")

wza_win_env6 <- read_csv("data/genomic_data/WZA_win_env6_bf.csv")
wza_win_env7 <- read_csv("data/genomic_data/WZA_win_env7_bf.csv")
wza_win_env8 <- read_csv("data/genomic_data/WZA_win_env8_bf.csv")
wza_win_env9 <- read_csv("data/genomic_data/WZA_win_env9_bf.csv")

#SNPs
snps_env1 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env1_bf.csv")
snps_env2 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env2_bf.csv")
snps_env3 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env3_bf.csv")
snps_env4 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env4_bf.csv")
snps_env5 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env5_bf.csv")

snps_env6 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env6_bf.csv")
snps_env7 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env7_bf.csv")
snps_env8 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env8_bf.csv")
snps_env9 <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/WZA_snps_env9_bf.csv")


#Filter by Bonferonnii correction alpha critical = 1.29423e-06, aka 5.887988 sigma
#mat_bon <- wza_win_mat %>% filter(Z_pVal<1.294264e-06)
env1_bon <- wza_win_env1 %>% filter(Z_pVal<1.294264e-06) #21 windows
env2_bon <- wza_win_env2 %>% filter(Z_pVal<1.294264e-06) #16 windows
env3_bon <- wza_win_env3 %>% filter(Z_pVal<1.294264e-06) #32 windows
env4_bon <- wza_win_env4 %>% filter(Z_pVal<1.294264e-06) #28 windows
env5_bon <- wza_win_env5 %>% filter(Z_pVal<1.294264e-06) #12 windows

env6_bon <- wza_win_env6 %>% filter(Z_pVal<1.294264e-06) #12 windows
env7_bon <- wza_win_env7 %>% filter(Z_pVal<1.294264e-06) #27 windows
env8_bon <- wza_win_env8 %>% filter(Z_pVal<1.294264e-06) #18 windows
env9_bon <- wza_win_env9 %>% filter(Z_pVal<1.294264e-06) #19 windows
################################################################################################################
################################################################################################################

#Get peak windows only
#Filter out windows that are not peak windows
#Keep two adjoining windows if they are less than 0.3 log P-value units away from each other.
env1_bon_peak <- env1_bon %>% filter(!win %in% c(6458,6459)) #19 windows
env2_bon_peak <- env2_bon %>% filter(!win %in% c(3881,3882)) #14 windows
env3_bon_peak <- env3_bon %>% filter(!win %in% c(10527,30337)) #30 windows
env4_bon_peak <- env4_bon %>% filter(!win %in% c(6451,6457,6458,6468)) #24 windows
env5_bon_peak <- env5_bon # 12 windows

env6_bon_peak <- env6_bon %>% filter(!win %in% c(9939)) #11 windows
env7_bon_peak <- env7_bon %>% filter(!win %in% c(6451,6458,6459,6462)) #23 windows
env8_bon_peak <- env8_bon %>% filter(!win %in% c(6457)) #18 windows
env9_bon_peak <- env9_bon %>% filter(!win %in% c(3881,3882)) #19 windows


#Filter SNPs for peak
snps_env1_peak <- snps_env1 %>% filter(win %in% env1_bon_peak$win)
snps_env2_peak <- snps_env2 %>% filter(win %in% env2_bon_peak$win)
snps_env3_peak <- snps_env3 %>% filter(win %in% env3_bon_peak$win)
snps_env4_peak <- snps_env4 %>% filter(win %in% env4_bon_peak$win)
snps_env5_peak <- snps_env5 %>% filter(win %in% env5_bon_peak$win)

snps_env6_peak <- snps_env6 %>% filter(win %in% env6_bon_peak$win)
snps_env7_peak <- snps_env7 %>% filter(win %in% env7_bon_peak$win)
snps_env8_peak <- snps_env8 %>% filter(win %in% env8_bon_peak$win)
snps_env9_peak <- snps_env9 %>% filter(win %in% env9_bon_peak$win)


#Add chr_snp column
snps_env1_peak <- snps_env1_peak %>% unite(chr_snp,"chr","snp",sep="_")
snps_env2_peak <- snps_env2_peak %>% unite(chr_snp,"chr","snp",sep="_")
snps_env3_peak <- snps_env3_peak %>% unite(chr_snp,"chr","snp",sep="_")
snps_env4_peak <- snps_env4_peak %>% unite(chr_snp,"chr","snp",sep="_")
snps_env5_peak <- snps_env5_peak %>% unite(chr_snp,"chr","snp",sep="_")
snps_env6_peak <- snps_env6_peak %>% unite(chr_snp,"chr","snp",sep="_")
snps_env7_peak <- snps_env7_peak %>% unite(chr_snp,"chr","snp",sep="_")
snps_env8_peak <- snps_env8_peak %>% unite(chr_snp,"chr","snp",sep="_")
snps_env9_peak <- snps_env9_peak %>% unite(chr_snp,"chr","snp",sep="_")




#Export peak snps
write_csv(snps_env1_peak,"data/genomic_data/snps_peak_env1_bf.csv")
write_csv(snps_env2_peak,"data/genomic_data/snps_peak_env2_bf.csv")
write_csv(snps_env3_peak,"data/genomic_data/snps_peak_env3_bf.csv")
write_csv(snps_env4_peak,"data/genomic_data/snps_peak_env4_bf.csv")
write_csv(snps_env5_peak,"data/genomic_data/snps_peak_env5_bf.csv")
write_csv(snps_env6_peak,"data/genomic_data/snps_peak_env6_bf.csv")
write_csv(snps_env7_peak,"data/genomic_data/snps_peak_env7_bf.csv")
write_csv(snps_env8_peak,"data/genomic_data/snps_peak_env8_bf.csv")
write_csv(snps_env9_peak,"data/genomic_data/snps_peak_env9_bf.csv")

#Export peak windows
write_csv(env1_bon_peak,"data/genomic_data/peak_window_env1_bf.csv")
write_csv(env2_bon_peak,"data/genomic_data/peak_window_env2_bf.csv")
write_csv(env3_bon_peak,"data/genomic_data/peak_window_env3_bf.csv")
write_csv(env4_bon_peak,"data/genomic_data/peak_window_env4_bf.csv")
write_csv(env5_bon_peak,"data/genomic_data/peak_window_env5_bf.csv")
write_csv(env6_bon_peak,"data/genomic_data/peak_window_env6_bf.csv")
write_csv(env7_bon_peak,"data/genomic_data/peak_window_env7_bf.csv")
write_csv(env8_bon_peak,"data/genomic_data/peak_window_env8_bf.csv")
write_csv(env9_bon_peak,"data/genomic_data/peak_window_env9_bf.csv")















