##################################################################################
## Filter out BF>0 SNPs for all env
## Select snp datatset that will be used for permutation test
## Author Daniel Anstett
## 
##
## Last Modified April 25, 2023
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)


#Import snp env associations
env1 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_1_trim.tsv",header=F, sep=" ")
env2 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_2_trim.tsv",header=F, sep=" ")
env3 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_3_trim.tsv",header=F, sep=" ")
env4 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_4_trim.tsv",header=F, sep=" ")
env5 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_5_trim.tsv",header=F, sep=" ")
env6 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_6_trim.tsv",header=F, sep=" ")
env7 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_7_trim.tsv",header=F, sep=" ")
env8 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_8_trim.tsv",header=F, sep=" ")
env9 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_9_trim.tsv",header=F, sep=" ")

#Name Columns
colnames(env1) <- c("Chromosome","SNP","Env","BF")
colnames(env2) <- c("Chromosome","SNP","Env","BF")
colnames(env3) <- c("Chromosome","SNP","Env","BF")
colnames(env4) <- c("Chromosome","SNP","Env","BF")
colnames(env5) <- c("Chromosome","SNP","Env","BF")
colnames(env6) <- c("Chromosome","SNP","Env","BF")
colnames(env7) <- c("Chromosome","SNP","Env","BF")
colnames(env8) <- c("Chromosome","SNP","Env","BF")
colnames(env9) <- c("Chromosome","SNP","Env","BF")

#Filter Bayes factor >= 0
env1_filter <- env1 %>% filter(BF >= 0)
env2_filter <- env2 %>% filter(BF >= 0)
env3_filter <- env3 %>% filter(BF >= 0)
env4_filter <- env4 %>% filter(BF >= 0)
env5_filter <- env5 %>% filter(BF >= 0)
env6_filter <- env6 %>% filter(BF >= 0)
env7_filter <- env7 %>% filter(BF >= 0)
env8_filter <- env8 %>% filter(BF >= 0)
env9_filter <- env9 %>% filter(BF >= 0)


#Make Chromosome SNP variable for easier left_joining later
env1_united <- env1_filter %>% unite(chr_snp,Chromosome,SNP)
env1_united <- env1_united %>% select(chr_snp)
env1_filter_united <- cbind(env1_united,env1_filter)

env2_united <- env2_filter %>% unite(chr_snp,Chromosome,SNP)
env2_united <- env2_united %>% select(chr_snp)
env2_filter_united <- cbind(env2_united,env2_filter)

env3_united <- env3_filter %>% unite(chr_snp,Chromosome,SNP)
env3_united <- env3_united %>% select(chr_snp)
env3_filter_united <- cbind(env3_united,env3_filter)

env4_united <- env4_filter %>% unite(chr_snp,Chromosome,SNP)
env4_united <- env4_united %>% select(chr_snp)
env4_filter_united <- cbind(env4_united,env4_filter)

env5_united <- env5_filter %>% unite(chr_snp,Chromosome,SNP)
env5_united <- env5_united %>% select(chr_snp)
env5_filter_united <- cbind(env5_united,env5_filter)

env6_united <- env6_filter %>% unite(chr_snp,Chromosome,SNP)
env6_united <- env6_united %>% select(chr_snp)
env6_filter_united <- cbind(env6_united,env6_filter)

env7_united <- env7_filter %>% unite(chr_snp,Chromosome,SNP)
env7_united <- env7_united %>% select(chr_snp)
env7_filter_united <- cbind(env7_united,env7_filter)

env8_united <- env8_filter %>% unite(chr_snp,Chromosome,SNP)
env8_united <- env8_united %>% select(chr_snp)
env8_filter_united <- cbind(env8_united,env8_filter)

env9_united <- env9_filter %>% unite(chr_snp,Chromosome,SNP)
env9_united <- env9_united %>% select(chr_snp)
env9_filter_united <- cbind(env9_united,env9_filter)


#Export SNPs with Bf>0
write_csv(env1_filter_united, "data/env1_BF0.csv",col_names=TRUE)
write_csv(env2_filter_united, "data/env2_BF0.csv",col_names=TRUE)
write_csv(env3_filter_united, "data/env3_BF0.csv",col_names=TRUE)
write_csv(env4_filter_united, "data/env4_BF0.csv",col_names=TRUE)
write_csv(env5_filter_united, "data/env5_BF0.csv",col_names=TRUE)
write_csv(env6_filter_united, "data/env6_BF0.csv",col_names=TRUE)
write_csv(env7_filter_united, "data/env7_BF0.csv",col_names=TRUE)
write_csv(env8_filter_united, "data/env8_BF0.csv",col_names=TRUE)
write_csv(env9_filter_united, "data/env9_BF0.csv",col_names=TRUE)



