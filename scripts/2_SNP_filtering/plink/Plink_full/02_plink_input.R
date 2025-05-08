#############################################################################################################
## Calculate empirical p-value for PLINK dat SNP thinning
## Done on full WZA BF filtered SNPs (no removal of non-linear clim SNPs)
## Author Daniel Anstett
## 
## 
## Last Modified March 12, 2025
#############################################################################################################


library(tidyverse)

#Import BF
env1 <- read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/ENV_1_trim.tsv",header=F, sep="\t")
env2 <- read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/ENV_2_trim.tsv",header=F, sep="\t")
env3 <- read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/ENV_3_trim.tsv",header=F, sep="\t")
env4 <- read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/ENV_4_trim.tsv",header=F, sep="\t")
env5 <- read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/ENV_5_trim.tsv",header=F, sep="\t")
env6 <- read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/ENV_6_trim.tsv",header=F, sep="\t")
env7 <- read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/ENV_7_trim.tsv",header=F, sep="\t")
env8 <- read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/ENV_8_trim.tsv",header=F, sep="\t")
env9 <- read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/ENV_9_trim.tsv",header=F, sep="\t")

#Rename Columns
colnames(env1) <- c("Chromosome","SNP","Env","BF")
colnames(env2) <- c("Chromosome","SNP","Env","BF")
colnames(env3) <- c("Chromosome","SNP","Env","BF")
colnames(env4) <- c("Chromosome","SNP","Env","BF")
colnames(env5) <- c("Chromosome","SNP","Env","BF")
colnames(env6) <- c("Chromosome","SNP","Env","BF")
colnames(env7) <- c("Chromosome","SNP","Env","BF")
colnames(env8) <- c("Chromosome","SNP","Env","BF")
colnames(env9) <- c("Chromosome","SNP","Env","BF")

#Make SNP ID have chr number
env1 <- env1 %>% unite(col="chr_snp",c("Chromosome","SNP"),sep="_",remove = FALSE)
env2 <- env2 %>% unite(col="chr_snp",c("Chromosome","SNP"),sep="_",remove = FALSE)
env3 <- env3 %>% unite(col="chr_snp",c("Chromosome","SNP"),sep="_",remove = FALSE)
env4 <- env4 %>% unite(col="chr_snp",c("Chromosome","SNP"),sep="_",remove = FALSE)
env5 <- env5 %>% unite(col="chr_snp",c("Chromosome","SNP"),sep="_",remove = FALSE)
env6 <- env6 %>% unite(col="chr_snp",c("Chromosome","SNP"),sep="_",remove = FALSE)
env7 <- env7 %>% unite(col="chr_snp",c("Chromosome","SNP"),sep="_",remove = FALSE)
env8 <- env8 %>% unite(col="chr_snp",c("Chromosome","SNP"),sep="_",remove = FALSE)
env9 <- env9 %>% unite(col="chr_snp",c("Chromosome","SNP"),sep="_",remove = FALSE)

#Calc empirical p-value based on BF rank of all SNP
env1$empirical_p <- (1-rank(env1$BF)/length(env1$BF)) 
env2$empirical_p <- (1-rank(env2$BF)/length(env2$BF)) 
env3$empirical_p <- (1-rank(env3$BF)/length(env3$BF)) 
env4$empirical_p <- (1-rank(env4$BF)/length(env4$BF)) 
env5$empirical_p <- (1-rank(env5$BF)/length(env5$BF)) 
env6$empirical_p <- (1-rank(env6$BF)/length(env6$BF)) 
env7$empirical_p <- (1-rank(env7$BF)/length(env7$BF)) 
env8$empirical_p <- (1-rank(env8$BF)/length(env8$BF)) 
env9$empirical_p <- (1-rank(env9$BF)/length(env9$BF)) 

#Import SNP set per env variable. These have been filtered for Bayes, WZA cut offs
snp_set_env1 <- read_csv("data/genomic_data/snp_set_env1.csv") %>% select(chr_snp)
snp_set_env2 <- read_csv("data/genomic_data/snp_set_env2.csv") %>% select(chr_snp)
snp_set_env3 <- read_csv("data/genomic_data/snp_set_env3.csv") %>% select(chr_snp)
snp_set_env4 <- read_csv("data/genomic_data/snp_set_env4.csv") %>% select(chr_snp)
snp_set_env5 <- read_csv("data/genomic_data/snp_set_env5.csv") %>% select(chr_snp)
snp_set_env6 <- read_csv("data/genomic_data/snp_set_env6.csv") %>% select(chr_snp)
snp_set_env7 <- read_csv("data/genomic_data/snp_set_env7.csv") %>% select(chr_snp)
snp_set_env8 <- read_csv("data/genomic_data/snp_set_env8.csv") %>% select(chr_snp)
snp_set_env9 <- read_csv("data/genomic_data/snp_set_env9.csv") %>% select(chr_snp)

#Import SNP abundance table for 55 baseline populations per env variable
snp_set_abundance_env1 <- read_csv("data/genomic_data/snp_set_abundance_full_env1.csv")
snp_set_abundance_env2 <- read_csv("data/genomic_data/snp_set_abundance_full_env2.csv") 
snp_set_abundance_env3 <- read_csv("data/genomic_data/snp_set_abundance_full_env3.csv") 
snp_set_abundance_env4 <- read_csv("data/genomic_data/snp_set_abundance_full_env4.csv") 
snp_set_abundance_env5 <- read_csv("data/genomic_data/snp_set_abundance_full_env5.csv") 
snp_set_abundance_env6 <- read_csv("data/genomic_data/snp_set_abundance_full_env6.csv") 
snp_set_abundance_env7 <- read_csv("data/genomic_data/snp_set_abundance_full_env7.csv") 
snp_set_abundance_env8 <- read_csv("data/genomic_data/snp_set_abundance_full_env8.csv") 
snp_set_abundance_env9 <- read_csv("data/genomic_data/snp_set_abundance_full_env9.csv") 

#############################################################################################################
#Filter large DF with SNP set
snp_set_env1_p <- snp_set_env1 %>% left_join(env1,by=c("chr_snp"))
snp_set_env2_p <- snp_set_env2 %>% left_join(env2,by=c("chr_snp"))
snp_set_env3_p <- snp_set_env3 %>% left_join(env3,by=c("chr_snp"))
snp_set_env4_p <- snp_set_env4 %>% left_join(env4,by=c("chr_snp"))
snp_set_env5_p <- snp_set_env5 %>% left_join(env5,by=c("chr_snp"))
snp_set_env6_p <- snp_set_env6 %>% left_join(env6,by=c("chr_snp"))
snp_set_env7_p <- snp_set_env7 %>% left_join(env7,by=c("chr_snp"))
snp_set_env8_p <- snp_set_env8 %>% left_join(env8,by=c("chr_snp"))
snp_set_env9_p <- snp_set_env9 %>% left_join(env9,by=c("chr_snp"))

#Add abundance data per env
snp_set_env1_all <- snp_set_env1_p %>% left_join(snp_set_abundance_env1,by=c("chr_snp"))
snp_set_env2_all <- snp_set_env2_p %>% left_join(snp_set_abundance_env2,by=c("chr_snp"))
snp_set_env3_all <- snp_set_env3_p %>% left_join(snp_set_abundance_env3,by=c("chr_snp"))
snp_set_env4_all <- snp_set_env4_p %>% left_join(snp_set_abundance_env4,by=c("chr_snp"))
snp_set_env5_all <- snp_set_env5_p %>% left_join(snp_set_abundance_env5,by=c("chr_snp"))
snp_set_env6_all <- snp_set_env6_p %>% left_join(snp_set_abundance_env6,by=c("chr_snp"))
snp_set_env7_all <- snp_set_env7_p %>% left_join(snp_set_abundance_env7,by=c("chr_snp"))
snp_set_env8_all <- snp_set_env8_p %>% left_join(snp_set_abundance_env8,by=c("chr_snp"))
snp_set_env9_all <- snp_set_env9_p %>% left_join(snp_set_abundance_env9,by=c("chr_snp"))

#Export as table delimited file
write_tsv(snp_set_env1_all,"data/genomic_data/plink_input_full_env1.tsv")
write_tsv(snp_set_env2_all,"data/genomic_data/plink_input_full_env2.tsv")
write_tsv(snp_set_env3_all,"data/genomic_data/plink_input_full_env3.tsv")
write_tsv(snp_set_env4_all,"data/genomic_data/plink_input_full_env4.tsv")
write_tsv(snp_set_env5_all,"data/genomic_data/plink_input_full_env5.tsv")
write_tsv(snp_set_env6_all,"data/genomic_data/plink_input_full_env6.tsv")
write_tsv(snp_set_env7_all,"data/genomic_data/plink_input_full_env7.tsv")
write_tsv(snp_set_env8_all,"data/genomic_data/plink_input_full_env8.tsv")
write_tsv(snp_set_env9_all,"data/genomic_data/plink_input_full_env9.tsv")

#Run plink clump per env variable after generating .bed files for baseline vcfs



