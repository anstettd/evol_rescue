##################################################################################
## Test for SNP overlap between GEA SNP set and high S SNP set
## 
## Author Daniel Anstett
## 
## 
## 
## Last Modified 07072025
###################################################################################
#Import libraries
library(tidyverse)

#Import GEA SNP set
strong_snp_set_clump <- read_csv("data/genomic_data/strong_snp_set_clump.csv")

#Import top plink clumped 215 SNPs per pop
p1_215 <- read_csv("data/genomic_data/p01_215.csv")
p2_215 <- read_csv("data/genomic_data/p02_215.csv")
p3_215 <- read_csv("data/genomic_data/p03_215.csv")
p4_215 <- read_csv("data/genomic_data/p04_215.csv")
p5_215 <- read_csv("data/genomic_data/p05_215.csv")
p6_215 <- read_csv("data/genomic_data/p06_215.csv")
p7_215 <- read_csv("data/genomic_data/p07_215.csv")
p8_215 <- read_csv("data/genomic_data/p08_215.csv")
p9_215 <- read_csv("data/genomic_data/p09_215.csv")
p10_215 <- read_csv("data/genomic_data/p10_215.csv")
p11_215 <- read_csv("data/genomic_data/p11_215.csv")

#Calculate overlap #No SNPs in common
common_snps_p01 <- intersect(strong_snp_set_clump$chr_snp, p1_215$snp_ID)
common_snps_p02 <- intersect(strong_snp_set_clump$chr_snp, p2_215$snp_ID)
common_snps_p03 <- intersect(strong_snp_set_clump$chr_snp, p3_215$snp_ID)
common_snps_p04 <- intersect(strong_snp_set_clump$chr_snp, p4_215$snp_ID)
common_snps_p05 <- intersect(strong_snp_set_clump$chr_snp, p5_215$snp_ID)
common_snps_p06 <- intersect(strong_snp_set_clump$chr_snp, p6_215$snp_ID)
common_snps_p07 <- intersect(strong_snp_set_clump$chr_snp, p7_215$snp_ID)
common_snps_p08 <- intersect(strong_snp_set_clump$chr_snp, p8_215$snp_ID)
common_snps_p09 <- intersect(strong_snp_set_clump$chr_snp, p9_215$snp_ID)
common_snps_p10 <- intersect(strong_snp_set_clump$chr_snp, p10_215$snp_ID)
common_snps_p11 <- intersect(strong_snp_set_clump$chr_snp, p11_215$snp_ID)

#Bind all SNPs
all_215 <- rbind(p1_215,
                 p2_215,
                 p3_215,
                 p4_215,
                 p5_215,
                 p6_215,
                 p7_215,
                 p8_215,
                 p9_215,
                 p10_215,
                 p11_215)

#Get site overlap per SNP
snp_site_counts <- all_215 %>%
  group_by(snp_ID) %>%
  summarise(Site_Count = n_distinct(Site)) %>%
  arrange(desc(Site_Count)) # all but 2 SNPs are unique SNPs suggesting they are involved in local adaptation and not range-wide climate adaptation





