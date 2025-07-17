##################################################################################
## Merge 22 larger data frames
## Filter to remove high SE slopes
## Export
## Author Daniel Anstett
## 
## 
## Last Modified 07072025
###################################################################################
#Import libraries
library(tidyverse)


#Import files
snp_1 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_01.csv") %>% filter(SE<5)
#plot(snp_1$Slope,snp_1$SE)

snp_2 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_02.csv") %>% filter(SE<5)
snp_3 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_03.csv") %>% filter(SE<5)
snp_4 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_04.csv") %>% filter(SE<5)
snp_5 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_05.csv") %>% filter(SE<5)
snp_6 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_06.csv") %>% filter(SE<5)
snp_7 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_07.csv") %>% filter(SE<5)
snp_8 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_08.csv") %>% filter(SE<5)
snp_9 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_09.csv") %>% filter(SE<5)
snp_10 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_10.csv") %>% filter(SE<5)
snp_11 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_11.csv") %>% filter(SE<5)
snp_12 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_12.csv") %>% filter(SE<5)
snp_13 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_13.csv") %>% filter(SE<5)
snp_14 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_14.csv") %>% filter(SE<5)
snp_15 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_15.csv") %>% filter(SE<5)
snp_16 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_16.csv") %>% filter(SE<5)
snp_17 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_17.csv") %>% filter(SE<5)
snp_18 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_18.csv") %>% filter(SE<5)
snp_19 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_19.csv") %>% filter(SE<5)
snp_20 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_20.csv") %>% filter(SE<5)
snp_21 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_21.csv") %>% filter(SE<5)
snp_22 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_glm_22.csv") %>% filter(SE<5)

#Import Baseline data
#import full snp table
pop_order_base<-read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", 
                           header=F, sep="\t")
snp_base<-read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", 
                     header=F, sep=" ")
loci_base<-read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", 
                      header=F, sep="\t")
colnames(loci_base) <- c("Chromosome","SNP")
loci_united_base <- loci_base %>% unite(chr_snp,"Chromosome","SNP",sep="_")
loci_snp_base <-cbind(loci_united_base,snp_base) #add snp lables to rows

###################################################################################
#Marge files
snp_all <- rbind(snp_1, snp_2, snp_3, snp_4, snp_5, snp_6, snp_7, snp_8, snp_9, snp_10,
                 snp_11, snp_12, snp_13, snp_14, snp_15, snp_16, snp_17, snp_18, snp_19, snp_20,
                 snp_21, snp_22) %>% mutate (Slope_abs=abs(Slope)) %>% 
  separate(snp_ID, into = c("prefix", "chr", "Position"), sep = "_", remove = FALSE) %>% 
  select(-prefix) %>%
  mutate(Position = as.numeric(Position))

dim(snp_all)
#Filter out SNPs that are only in timeseries
snp_all <- snp_all %>% filter (snp_ID %in% as.character(loci_snp_base$chr_snp))
dim(snp_all)
write_csv(snp_all,"/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/snp_all.csv")







