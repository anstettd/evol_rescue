#############################################################################################################
## Make LD thinned SNP set
## Author Daniel Anstett
## For All high S SNPs
##
## 
## Last Modified 07172025
#############################################################################################################


library(tidyverse)

#Note .clumped files from PLINK will need to be turned into .csv files. Use text editor and copy into excel, then save as .csv
#Import the plink clump LD thinned SNP set per env
highS_p01<-read_csv("data/genomic_data/highS_p01.clumped.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP,P) %>% rename(snp_ID = SNP) 
highS_p02<-read_csv("data/genomic_data/highS_p02.clumped.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP,P) %>% rename(snp_ID = SNP) 
highS_p03<-read_csv("data/genomic_data/highS_p03.clumped.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP,P) %>% rename(snp_ID = SNP) 
highS_p04<-read_csv("data/genomic_data/highS_p04.clumped.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP,P) %>% rename(snp_ID = SNP) 
highS_p05<-read_csv("data/genomic_data/highS_p05.clumped.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP,P) %>% rename(snp_ID = SNP) 
highS_p06<-read_csv("data/genomic_data/highS_p06.clumped.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP,P) %>% rename(snp_ID = SNP) 
highS_p07<-read_csv("data/genomic_data/highS_p07.clumped.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP,P) %>% rename(snp_ID = SNP) 
highS_p08<-read_csv("data/genomic_data/highS_p08.clumped.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP,P) %>% rename(snp_ID = SNP) 
highS_p09<-read_csv("data/genomic_data/highS_p09.clumped.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP,P) %>% rename(snp_ID = SNP) 
highS_p10<-read_csv("data/genomic_data/highS_p10.clumped.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP,P) %>% rename(snp_ID = SNP) 
highS_p11<-read_csv("data/genomic_data/highS_p11.clumped.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP,P) %>% rename(snp_ID = SNP) 

#Import data frame with highS
#Import
snp_all <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_all.csv") %>% filter(SE<5)

#Split into data for each pop
p01 <- snp_all %>% filter(Site==1) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p02 <- snp_all %>% filter(Site == 2) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p03 <- snp_all %>% filter(Site == 3) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p04 <- snp_all %>% filter(Site == 4) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p05 <- snp_all %>% filter(Site == 5) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p06 <- snp_all %>% filter(Site == 6) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p07 <- snp_all %>% filter(Site == 7) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p08 <- snp_all %>% filter(Site == 8) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p09 <- snp_all %>% filter(Site == 9) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p10 <- snp_all %>% filter(Site == 10) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p11 <- snp_all %>% filter(Site == 11) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
rm(snp_all)

#Get top 215 SNPs per pop
p01_215 <- left_join(highS_p01, p01, by = "snp_ID") %>% arrange(desc(Slope_abs)) %>% slice_head(n = 215)
p02_215 <- left_join(highS_p02, p02, by = "snp_ID") %>% arrange(desc(Slope_abs)) %>% slice_head(n = 215)
p03_215 <- left_join(highS_p03, p03, by = "snp_ID") %>% arrange(desc(Slope_abs)) %>% slice_head(n = 215)
p04_215 <- left_join(highS_p04, p04, by = "snp_ID") %>% arrange(desc(Slope_abs)) %>% slice_head(n = 215)
p05_215 <- left_join(highS_p05, p05, by = "snp_ID") %>% arrange(desc(Slope_abs)) %>% slice_head(n = 215)
p06_215 <- left_join(highS_p06, p06, by = "snp_ID") %>% arrange(desc(Slope_abs)) %>% slice_head(n = 215)
p07_215 <- left_join(highS_p07, p07, by = "snp_ID") %>% arrange(desc(Slope_abs)) %>% slice_head(n = 215)
p08_215 <- left_join(highS_p08, p08, by = "snp_ID") %>% arrange(desc(Slope_abs)) %>% slice_head(n = 215)
p09_215 <- left_join(highS_p09, p09, by = "snp_ID") %>% arrange(desc(Slope_abs)) %>% slice_head(n = 215)
p10_215 <- left_join(highS_p10, p10, by = "snp_ID") %>% arrange(desc(Slope_abs)) %>% slice_head(n = 215)
p11_215 <- left_join(highS_p11, p11, by = "snp_ID") %>% arrange(desc(Slope_abs)) %>% slice_head(n = 215)

#Write out snp sets for each population
write_csv(p01_215, "data/genomic_data/p01_215.csv")
write_csv(p02_215, "data/genomic_data/p02_215.csv")
write_csv(p03_215, "data/genomic_data/p03_215.csv")
write_csv(p04_215, "data/genomic_data/p04_215.csv")
write_csv(p05_215, "data/genomic_data/p05_215.csv")
write_csv(p06_215, "data/genomic_data/p06_215.csv")
write_csv(p07_215, "data/genomic_data/p07_215.csv")
write_csv(p08_215, "data/genomic_data/p08_215.csv")
write_csv(p09_215, "data/genomic_data/p09_215.csv")
write_csv(p10_215, "data/genomic_data/p10_215.csv")
write_csv(p11_215, "data/genomic_data/p11_215.csv")






