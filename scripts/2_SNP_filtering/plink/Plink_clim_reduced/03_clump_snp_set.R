#############################################################################################################
## Make LD thinned SNP set
## Author Daniel Anstett
## For All ENV & for CMD + Winter Precipitation
## For use in Response to Selection
## 
## Last Modified April 25, 2025
#############################################################################################################


library(tidyverse)

#Import the plink clump LD thinned SNP set per env
clump_env1<-read_csv("data/genomic_data/baseline_env1_clumped.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP) %>% rename(chr_snp = SNP) 
clump_env2<-read_csv("data/genomic_data/baseline_env2_clumped.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP) %>% rename(chr_snp = SNP) 
clump_env3<-read_csv("data/genomic_data/baseline_env3_clumped.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP) %>% rename(chr_snp = SNP) 
clump_env4<-read_csv("data/genomic_data/baseline_env4_clumped.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP) %>% rename(chr_snp = SNP) 
clump_env5<-read_csv("data/genomic_data/baseline_env5_clumped.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP) %>% rename(chr_snp = SNP) 
clump_env6<-read_csv("data/genomic_data/baseline_env6_clumped.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP) %>% rename(chr_snp = SNP) 
clump_env7<-read_csv("data/genomic_data/baseline_env7_clumped.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP) %>% rename(chr_snp = SNP) 
clump_env8<-read_csv("data/genomic_data/baseline_env8_clumped.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP) %>% rename(chr_snp = SNP) 
clump_env9<-read_csv("data/genomic_data/baseline_env9_clumped.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP) %>% rename(chr_snp = SNP) 

#Merge all 9 env
strong_snp_set_clump_raw <- rbind(clump_env1,
                              clump_env2,
                              clump_env3,
                              clump_env4,
                              clump_env5,
                              clump_env6,
                              clump_env7,
                              clump_env8,
                              clump_env9)
#Get Unique SNP
strong_snp_set_clump <- as.data.frame(unique(strong_snp_set_clump_raw$chr_snp))  %>%
  rename(chr_snp = `unique(strong_snp_set_clump_raw$chr_snp)`) %>% arrange(chr_snp)

#Merge all 9 env
env58_snp_set_clump_raw <- rbind(clump_env5,clump_env8)

#Get Unique SNP
env58_snp_set_clump <- as.data.frame(unique(env58_snp_set_clump_raw$chr_snp))  %>%
  rename(chr_snp = `unique(env58_snp_set_clump_raw$chr_snp)`) %>% arrange(chr_snp)

#Write out joined sets
write_csv(strong_snp_set_clump,"data/genomic_data/strong_snp_set_clump.csv")
write_csv(env58_snp_set_clump,"data/genomic_data/strong_snp_set_clump_env58.csv")

#Write out snp sets for individual env variables
write_csv(clump_env1,"data/genomic_data/strong_snp_set_clump_env1.csv")
write_csv(clump_env2,"data/genomic_data/strong_snp_set_clump_env2.csv")
write_csv(clump_env3,"data/genomic_data/strong_snp_set_clump_env3.csv")
write_csv(clump_env4,"data/genomic_data/strong_snp_set_clump_env4.csv")
write_csv(clump_env5,"data/genomic_data/strong_snp_set_clump_env5.csv")
write_csv(clump_env6,"data/genomic_data/strong_snp_set_clump_env6.csv")
write_csv(clump_env7,"data/genomic_data/strong_snp_set_clump_env7.csv")
write_csv(clump_env8,"data/genomic_data/strong_snp_set_clump_env8.csv")
write_csv(clump_env9,"data/genomic_data/strong_snp_set_clump_env9.csv")





