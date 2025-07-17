#############################################################################################################
## Make LD thinned SNP set
## Author Daniel Anstett
## For All high S SNPs
##
## 
## Last Modified 07172025
#############################################################################################################


library(tidyverse)

#Import the plink clump LD thinned SNP set per env
highS_p01<-read_csv("data/genomic_data/highS_p01.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP) %>% rename(chr_snp = SNP) 
highS_p02<-read_csv("data/genomic_data/highS_p02.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP) %>% rename(chr_snp = SNP) 
highS_p03<-read_csv("data/genomic_data/highS_p03.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP) %>% rename(chr_snp = SNP) 
highS_p04<-read_csv("data/genomic_data/highS_p04.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP) %>% rename(chr_snp = SNP) 
highS_p05<-read_csv("data/genomic_data/highS_p05.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP) %>% rename(chr_snp = SNP) 
highS_p06<-read_csv("data/genomic_data/highS_p06.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP) %>% rename(chr_snp = SNP) 
highS_p07<-read_csv("data/genomic_data/highS_p07.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP) %>% rename(chr_snp = SNP) 
highS_p08<-read_csv("data/genomic_data/highS_p08.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP) %>% rename(chr_snp = SNP) 
highS_p09<-read_csv("data/genomic_data/highS_p09.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP) %>% rename(chr_snp = SNP) 
highS_p10<-read_csv("data/genomic_data/highS_p10.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP) %>% rename(chr_snp = SNP) 
highS_p11<-read_csv("data/genomic_data/highS_p11.csv", col_names = T) %>% arrange(SNP) %>% 
  select(SNP) %>% rename(chr_snp = SNP) 

#Import data frame with highS
#Import
snp_all <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_all.csv") %>% filter(SE<5)

#Split into data for each pop
p1_215 <- snp_all %>% filter(Site==1) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p2_215 <- snp_all %>% filter(Site == 2) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p3_215 <- snp_all %>% filter(Site == 3) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p4_215 <- snp_all %>% filter(Site == 4) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p5_215 <- snp_all %>% filter(Site == 5) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p6_215 <- snp_all %>% filter(Site == 6) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p7_215 <- snp_all %>% filter(Site == 7) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p8_215 <- snp_all %>% filter(Site == 8) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p9_215 <- snp_all %>% filter(Site == 9) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p10_215 <- snp_all %>% filter(Site == 10) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
p11_215 <- snp_all %>% filter(Site == 11) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 2000)
rm(snp_all)

#%>% slice_head(n = 215)



#Get top 215 SNPs per pop
p1_215 <- p01 %>% left_join(highS_p01,p1_215, by=c(snp_ID=)) %>% arrange(desc(Slope_abs)) %>% slice_head(n = 215)
p2_215 <- p02 %>% 
p3_215 <- p03 %>% 
p4_215 <- p04 %>% 
p5_215 <- p05 %>% 
p6_215 <- p06 %>% 
p7_215 <- p07 %>% 
p8_215 <- p08 %>% 
p9_215 <- p09 %>% 
p10_215 <- p10 %>% 
p11_215 <- p11 %>% 

highS_p01 <- highS_p01 %>% filter()







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
dim(strong_snp_set_clump) # 105 unique SNPs after PLINK for climate-associated SNP set reduced to only have SNPs with a monotonic increase.

#Merge all 9 env
env58_snp_set_clump_raw <- rbind(clump_env5,clump_env7,clump_env8)
dim(env58_snp_set_clump_raw) 

#Get Unique SNP
env58_snp_set_clump <- as.data.frame(unique(env58_snp_set_clump_raw$chr_snp))  %>%
  rename(chr_snp = `unique(env58_snp_set_clump_raw$chr_snp)`) %>% arrange(chr_snp)
dim(env58_snp_set_clump) #103 unique SNPs after PLINK for SNPs highly assocaited in three climate variables (CMD, Tave_sm, PPT_wt)
#for climate-associated SNP set reduced to only have SNPs with a monotonic increase.

#Write out joined sets Not currently in use
#write_csv(strong_snp_set_clump,"data/genomic_data/strong_snp_set_clump.csv")
#write_csv(env58_snp_set_clump,"data/genomic_data/strong_snp_set_clump_env578.csv")

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





