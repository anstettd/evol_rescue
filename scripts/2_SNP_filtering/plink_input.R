#############################################################################################################
## Calculate empirical p-value for PLINK dat SNP thinning
## Author Daniel Anstett
## 
## 
## Last Modified March 12, 2025
#############################################################################################################


library(tidyverse)

#Import files
env1 <- read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/ENV_1_trim.tsv",header=F, sep="\t")
colnames(env1) <- c("Chromosome","SNP","Env","BF")

#env1_test <- env1[1:10,]
env1$empirical_p <- (1-rank(env1$BF)/length(env1$BF)) 


#snps_env1 <- read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/WZA_snps_env1_bf.csv")



