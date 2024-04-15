#############################################################################################################
## Generate Final SNP Set
## BF > 30
## Peak Boferroni window, BF >5
## 
## 
## Last Modified March 20, 2023
#############################################################################################################
#Import libraries
library(tidyverse)

#Import Large Loci Win
#loci_win <- read_csv("/Users/daniel_anstett/Dropbox/AM_Workshop/Large_files/loci_win.csv")                                                                                 
#loci_win <- loci_win %>% unite(col="chr_snp", c("chr","snp"), sep="_")

#Import files
env1_united <- read_csv("data/genomic_data/snps_peak_env1_bf.csv")
env2_united <- read_csv("data/genomic_data/snps_peak_env2_bf.csv")
env3_united <- read_csv("data/genomic_data/snps_peak_env3_bf.csv")
env4_united <- read_csv("data/genomic_data/snps_peak_env4_bf.csv")
env5_united <- read_csv("data/genomic_data/snps_peak_env5_bf.csv")
env6_united <- read_csv("data/genomic_data/snps_peak_env6_bf.csv")
env7_united <- read_csv("data/genomic_data/snps_peak_env7_bf.csv")
env8_united <- read_csv("data/genomic_data/snps_peak_env8_bf.csv")
env9_united <- read_csv("data/genomic_data/snps_peak_env9_bf.csv")

#Import all BF
#Import snp env associations (Baseline)
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

env1_allBF <- env1 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env2_allBF <- env2 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env3_allBF <- env3 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env4_allBF <- env4 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env5_allBF <- env5 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env6_allBF <- env6 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env7_allBF <- env7 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env8_allBF <- env8 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env9_allBF <- env9 %>% unite(chr_snp,"Chromosome","SNP",sep="_")


#############################################################################################################
#Investigate how many SNPs are left

#Select rows you need
snp_env1_peakbf_win <- env1_united %>% select(chr_snp,snp_c,win,p_bar,q_bar,MAF,Env,BF)
snp_env2_peakbf_win <- env2_united %>% select(chr_snp,snp_c,win,p_bar,q_bar,MAF,Env,BF)
snp_env3_peakbf_win <- env3_united %>% select(chr_snp,snp_c,win,p_bar,q_bar,MAF,Env,BF)
snp_env4_peakbf_win <- env4_united %>% select(chr_snp,snp_c,win,p_bar,q_bar,MAF,Env,BF)
snp_env5_peakbf_win <- env5_united %>% select(chr_snp,snp_c,win,p_bar,q_bar,MAF,Env,BF)
snp_env6_peakbf_win <- env6_united %>% select(chr_snp,snp_c,win,p_bar,q_bar,MAF,Env,BF)
snp_env7_peakbf_win <- env7_united %>% select(chr_snp,snp_c,win,p_bar,q_bar,MAF,Env,BF)
snp_env8_peakbf_win <- env8_united %>% select(chr_snp,snp_c,win,p_bar,q_bar,MAF,Env,BF)
snp_env9_peakbf_win <- env9_united %>% select(chr_snp,snp_c,win,p_bar,q_bar,MAF,Env,BF)

#How many windowed SNPs clear BF thresholds
snp_env1_peakbf30_win <- snp_env1_peakbf_win %>% filter(BF>30) 
snp_env2_peakbf30_win <- snp_env2_peakbf_win %>% filter(BF>30) 
snp_env3_peakbf30_win <- snp_env3_peakbf_win %>% filter(BF>30) 
snp_env4_peakbf30_win <- snp_env4_peakbf_win %>% filter(BF>30) 
snp_env5_peakbf30_win <- snp_env5_peakbf_win %>% filter(BF>30) 
snp_env6_peakbf30_win <- snp_env6_peakbf_win %>% filter(BF>30) 
snp_env7_peakbf30_win <- snp_env7_peakbf_win %>% filter(BF>30) 
snp_env8_peakbf30_win <- snp_env8_peakbf_win %>% filter(BF>30) 
snp_env9_peakbf30_win <- snp_env9_peakbf_win %>% filter(BF>30) 

snp_env1_peakbf20_win <- snp_env1_peakbf_win %>% filter(BF>20)
snp_env2_peakbf20_win <- snp_env2_peakbf_win %>% filter(BF>20) 
snp_env3_peakbf20_win <- snp_env3_peakbf_win %>% filter(BF>20) 
snp_env4_peakbf20_win <- snp_env4_peakbf_win %>% filter(BF>20) 
snp_env5_peakbf20_win <- snp_env5_peakbf_win %>% filter(BF>20) 
snp_env6_peakbf20_win <- snp_env6_peakbf_win %>% filter(BF>20) 
snp_env7_peakbf20_win <- snp_env7_peakbf_win %>% filter(BF>20) 
snp_env8_peakbf20_win <- snp_env8_peakbf_win %>% filter(BF>20) 
snp_env9_peakbf20_win <- snp_env9_peakbf_win %>% filter(BF>20) 

snp_env1_peakbf10_win <- snp_env1_peakbf_win %>% filter(BF>10)
snp_env2_peakbf10_win <- snp_env2_peakbf_win %>% filter(BF>10) 
snp_env3_peakbf10_win <- snp_env3_peakbf_win %>% filter(BF>10) 
snp_env4_peakbf10_win <- snp_env4_peakbf_win %>% filter(BF>10) 
snp_env5_peakbf10_win <- snp_env5_peakbf_win %>% filter(BF>10) 
snp_env6_peakbf10_win <- snp_env6_peakbf_win %>% filter(BF>10) 
snp_env7_peakbf10_win <- snp_env7_peakbf_win %>% filter(BF>10) 
snp_env8_peakbf10_win <- snp_env8_peakbf_win %>% filter(BF>10) 
snp_env9_peakbf10_win <- snp_env9_peakbf_win %>% filter(BF>10) 

snp_env1_peakbf5_win <- snp_env1_peakbf_win %>% filter(BF>5)
snp_env2_peakbf5_win <- snp_env2_peakbf_win %>% filter(BF>5) 
snp_env3_peakbf5_win <- snp_env3_peakbf_win %>% filter(BF>5) 
snp_env4_peakbf5_win <- snp_env4_peakbf_win %>% filter(BF>5) 
snp_env5_peakbf5_win <- snp_env5_peakbf_win %>% filter(BF>5) 
snp_env6_peakbf5_win <- snp_env6_peakbf_win %>% filter(BF>5) 
snp_env7_peakbf5_win <- snp_env7_peakbf_win %>% filter(BF>5) 
snp_env8_peakbf5_win <- snp_env8_peakbf_win %>% filter(BF>5) 
snp_env9_peakbf5_win <- snp_env9_peakbf_win %>% filter(BF>5) 

snp_env1_peakbf2_win <- snp_env1_peakbf_win %>% filter(BF>2)
snp_env2_peakbf2_win <- snp_env2_peakbf_win %>% filter(BF>2) 
snp_env3_peakbf2_win <- snp_env3_peakbf_win %>% filter(BF>2) 
snp_env4_peakbf2_win <- snp_env4_peakbf_win %>% filter(BF>2) 
snp_env5_peakbf2_win <- snp_env5_peakbf_win %>% filter(BF>2) 
snp_env6_peakbf2_win <- snp_env6_peakbf_win %>% filter(BF>2) 
snp_env7_peakbf2_win <- snp_env7_peakbf_win %>% filter(BF>2) 
snp_env8_peakbf2_win <- snp_env8_peakbf_win %>% filter(BF>2) 
snp_env9_peakbf2_win <- snp_env9_peakbf_win %>% filter(BF>2) 

#BF>30 for all windows
env1_united_bf30 <- env1_allBF %>% filter(BF>30)
env2_united_bf30 <- env2_allBF %>% filter(BF>30)
env3_united_bf30 <- env3_allBF %>% filter(BF>30)
env4_united_bf30 <- env4_allBF %>% filter(BF>30)
env5_united_bf30 <- env5_allBF %>% filter(BF>30)
env6_united_bf30 <- env6_allBF %>% filter(BF>30)
env7_united_bf30 <- env7_allBF %>% filter(BF>30)
env8_united_bf30 <- env8_allBF %>% filter(BF>30)
env9_united_bf30 <- env9_allBF %>% filter(BF>30)

#BF>20 for all windows
env1_united_bf20 <- env1_allBF %>% filter(BF>20)
env2_united_bf20 <- env2_allBF %>% filter(BF>20)
env3_united_bf20 <- env3_allBF %>% filter(BF>20)
env4_united_bf20 <- env4_allBF %>% filter(BF>20)
env5_united_bf20 <- env5_allBF %>% filter(BF>20)
env6_united_bf20 <- env6_allBF %>% filter(BF>20)
env7_united_bf20 <- env7_allBF %>% filter(BF>20)
env8_united_bf20 <- env8_allBF %>% filter(BF>20)
env9_united_bf20 <- env9_allBF %>% filter(BF>20)

#Make SNPs in Window Summary Data Frame
snps_in_win <- data.frame()

snps_in_win[1,1] <- ">30 all win"
snps_in_win[1,2] <- dim(env1_united_bf30)[1]
snps_in_win[1,3] <- dim(env2_united_bf30)[1]
snps_in_win[1,4] <- dim(env3_united_bf30)[1]
snps_in_win[1,5] <- dim(env4_united_bf30)[1]
snps_in_win[1,6] <- dim(env5_united_bf30)[1]
snps_in_win[1,7] <- dim(env6_united_bf30)[1]
snps_in_win[1,8] <- dim(env7_united_bf30)[1]
snps_in_win[1,9] <- dim(env8_united_bf30)[1]
snps_in_win[1,10] <- dim(env9_united_bf30)[1]

snps_in_win[2,1] <- ">30"
snps_in_win[2,2] <- dim(snp_env1_peakbf30_win)[1]
snps_in_win[2,3] <- dim(snp_env2_peakbf30_win)[1]
snps_in_win[2,4] <- dim(snp_env3_peakbf30_win)[1]
snps_in_win[2,5] <- dim(snp_env4_peakbf30_win)[1]
snps_in_win[2,6] <- dim(snp_env5_peakbf30_win)[1]
snps_in_win[2,7] <- dim(snp_env6_peakbf30_win)[1]
snps_in_win[2,8] <- dim(snp_env7_peakbf30_win)[1]
snps_in_win[2,9] <- dim(snp_env8_peakbf30_win)[1]
snps_in_win[2,10] <- dim(snp_env9_peakbf30_win)[1]

snps_in_win[3,1] <- ">20"
snps_in_win[3,2] <- dim(snp_env1_peakbf20_win)[1]
snps_in_win[3,3] <- dim(snp_env2_peakbf20_win)[1]
snps_in_win[3,4] <- dim(snp_env3_peakbf20_win)[1]
snps_in_win[3,5] <- dim(snp_env4_peakbf20_win)[1]
snps_in_win[3,6] <- dim(snp_env5_peakbf20_win)[1]
snps_in_win[3,7] <- dim(snp_env6_peakbf20_win)[1]
snps_in_win[3,8] <- dim(snp_env7_peakbf20_win)[1]
snps_in_win[3,9] <- dim(snp_env8_peakbf20_win)[1]
snps_in_win[3,10] <- dim(snp_env9_peakbf20_win)[1]

snps_in_win[4,1] <- ">10"
snps_in_win[4,2] <- dim(snp_env1_peakbf10_win)[1]
snps_in_win[4,3] <- dim(snp_env2_peakbf10_win)[1]
snps_in_win[4,4] <- dim(snp_env3_peakbf10_win)[1]
snps_in_win[4,5] <- dim(snp_env4_peakbf10_win)[1]
snps_in_win[4,6] <- dim(snp_env5_peakbf10_win)[1]
snps_in_win[4,7] <- dim(snp_env6_peakbf10_win)[1]
snps_in_win[4,8] <- dim(snp_env7_peakbf10_win)[1]
snps_in_win[4,9] <- dim(snp_env8_peakbf10_win)[1]
snps_in_win[4,10] <- dim(snp_env9_peakbf10_win)[1]

snps_in_win[5,1] <- ">5"
snps_in_win[5,2] <- dim(snp_env1_peakbf5_win)[1]
snps_in_win[5,3] <- dim(snp_env2_peakbf5_win)[1]
snps_in_win[5,4] <- dim(snp_env3_peakbf5_win)[1]
snps_in_win[5,5] <- dim(snp_env4_peakbf5_win)[1]
snps_in_win[5,6] <- dim(snp_env5_peakbf5_win)[1]
snps_in_win[5,7] <- dim(snp_env6_peakbf5_win)[1]
snps_in_win[5,8] <- dim(snp_env7_peakbf5_win)[1]
snps_in_win[5,9] <- dim(snp_env8_peakbf5_win)[1]
snps_in_win[5,10] <- dim(snp_env9_peakbf5_win)[1]

snps_in_win[6,1] <- ">2"
snps_in_win[6,2] <- dim(snp_env1_peakbf2_win)[1]
snps_in_win[6,3] <- dim(snp_env2_peakbf2_win)[1]
snps_in_win[6,4] <- dim(snp_env3_peakbf2_win)[1]
snps_in_win[6,5] <- dim(snp_env4_peakbf2_win)[1]
snps_in_win[6,6] <- dim(snp_env5_peakbf2_win)[1]
snps_in_win[6,7] <- dim(snp_env6_peakbf2_win)[1]
snps_in_win[6,8] <- dim(snp_env7_peakbf2_win)[1]
snps_in_win[6,9] <- dim(snp_env8_peakbf2_win)[1]
snps_in_win[6,10] <- dim(snp_env9_peakbf2_win)[1]

snps_in_win[7,1] <- ">20 all win" 
snps_in_win[7,2] <- dim(env1_united_bf20)[1]
snps_in_win[7,3] <- dim(env2_united_bf20)[1]
snps_in_win[7,4] <- dim(env3_united_bf20)[1]
snps_in_win[7,5] <- dim(env4_united_bf20)[1]
snps_in_win[7,6] <- dim(env5_united_bf20)[1]
snps_in_win[7,7] <- dim(env6_united_bf20)[1]
snps_in_win[7,8] <- dim(env7_united_bf20)[1]
snps_in_win[7,9] <- dim(env8_united_bf20)[1]
snps_in_win[7,10] <- dim(env9_united_bf20)[1]

#Unique SNPs
snps_bf30_all <- rbind(env1_united_bf30,
                       env2_united_bf30,
                       env3_united_bf30,
                       env4_united_bf30,
                       env5_united_bf30,
                       env6_united_bf30,
                       env7_united_bf30,
                       env8_united_bf30,
                       env9_united_bf30)
snps_bf30_all_u <- as.data.frame(unique(snps_bf30_all$chr_snp))

snps_bf30 <- rbind(snp_env1_peakbf30_win,
                   snp_env2_peakbf30_win,
                   snp_env3_peakbf30_win,
                   snp_env4_peakbf30_win,
                   snp_env5_peakbf30_win,
                   snp_env6_peakbf30_win,
                   snp_env7_peakbf30_win,
                   snp_env8_peakbf30_win,
                   snp_env9_peakbf30_win)
snps_bf30_u <- as.data.frame(unique(snps_bf30$chr_snp))

snps_bf20 <- rbind(snp_env1_peakbf20_win,
                   snp_env2_peakbf20_win,
                   snp_env3_peakbf20_win,
                   snp_env4_peakbf20_win,
                   snp_env5_peakbf20_win,
                   snp_env6_peakbf20_win,
                   snp_env7_peakbf20_win,
                   snp_env8_peakbf20_win,
                   snp_env9_peakbf20_win)
snps_bf20_u <- as.data.frame(unique(snps_bf20$chr_snp))

snps_bf10 <- rbind(snp_env1_peakbf10_win,
                   snp_env2_peakbf10_win,
                   snp_env3_peakbf10_win,
                   snp_env4_peakbf10_win,
                   snp_env5_peakbf10_win,
                   snp_env6_peakbf10_win,
                   snp_env7_peakbf10_win,
                   snp_env8_peakbf10_win,
                   snp_env9_peakbf10_win)
snps_bf10_u <- as.data.frame(unique(snps_bf10$chr_snp))

snps_bf5 <- rbind(snp_env1_peakbf5_win,
                   snp_env2_peakbf5_win,
                   snp_env3_peakbf5_win,
                   snp_env4_peakbf5_win,
                   snp_env5_peakbf5_win,
                   snp_env6_peakbf5_win,
                   snp_env7_peakbf5_win,
                   snp_env8_peakbf5_win,
                   snp_env9_peakbf5_win)
snps_bf5_u <- as.data.frame(unique(snps_bf5$chr_snp))

snps_bf2 <- rbind(snp_env1_peakbf2_win,
                   snp_env2_peakbf2_win,
                   snp_env3_peakbf2_win,
                   snp_env4_peakbf2_win,
                   snp_env5_peakbf2_win,
                   snp_env6_peakbf2_win,
                   snp_env7_peakbf2_win,
                   snp_env8_peakbf2_win,
                   snp_env9_peakbf2_win)
snps_bf2_u <- as.data.frame(unique(snps_bf2$chr_snp))

snps_bf20_all <- rbind(env1_united_bf20,
                       env2_united_bf20,
                       env3_united_bf20,
                       env4_united_bf20,
                       env5_united_bf20,
                       env6_united_bf20,
                       env7_united_bf20,
                       env8_united_bf20,
                       env9_united_bf20)
snps_bf20_all_u <- as.data.frame(unique(snps_bf20_all$chr_snp))

#Input total
snps_in_win[1,11] <- dim(snps_bf30_all)[1]
snps_in_win[2,11] <- dim(snps_bf30)[1]
snps_in_win[3,11] <- dim(snps_bf20)[1]
snps_in_win[4,11] <- dim(snps_bf10)[1]
snps_in_win[5,11] <- dim(snps_bf5)[1]
snps_in_win[6,11] <- dim(snps_bf2)[1]
snps_in_win[7,11] <- dim(snps_bf30_all)[1]

#Input unique
snps_in_win[1,12] <- dim(snps_bf30_all_u)[1]
snps_in_win[2,12] <- dim(snps_bf30_u)[1]
snps_in_win[3,12] <- dim(snps_bf20_u)[1]
snps_in_win[4,12] <- dim(snps_bf10_u)[1]
snps_in_win[5,12] <- dim(snps_bf5_u)[1]
snps_in_win[6,12] <- dim(snps_bf2_u)[1]
snps_in_win[7,12] <- dim(snps_bf20_all_u)[1]

colnames(snps_in_win) <- c("log10BF","env1","env2","env3","env4","env5","env6","env7","env8","env9","total","unique")

snps_in_win


####################################################################################################################
#Summarize how many windows are retained
windows_in <- data.frame()

windows_in[1,1] <- "any"
windows_in[1,2] <- 19
windows_in[1,3] <- 14
windows_in[1,4] <- 30
windows_in[1,5] <- 24
windows_in[1,6] <- 12

windows_in[1,7] <- 11
windows_in[1,8] <- 23
windows_in[1,9] <- 18
windows_in[1,10] <- 19


windows_in[2,1] <- ">30"
windows_in[2,2] <- length(unique(snp_env1_peakbf30_win$win))
windows_in[2,3] <- length(unique(snp_env2_peakbf30_win$win))
windows_in[2,4] <- length(unique(snp_env3_peakbf30_win$win))
windows_in[2,5] <- length(unique(snp_env4_peakbf30_win$win))
windows_in[2,6] <- length(unique(snp_env5_peakbf30_win$win))
windows_in[2,7] <- length(unique(snp_env6_peakbf30_win$win))
windows_in[2,8] <- length(unique(snp_env7_peakbf30_win$win))
windows_in[2,9] <- length(unique(snp_env8_peakbf30_win$win))
windows_in[2,10] <- length(unique(snp_env9_peakbf30_win$win))

windows_in[3,1] <- ">20"
windows_in[3,2] <- length(unique(snp_env1_peakbf20_win$win))
windows_in[3,3] <- length(unique(snp_env2_peakbf20_win$win))
windows_in[3,4] <- length(unique(snp_env3_peakbf20_win$win))
windows_in[3,5] <- length(unique(snp_env4_peakbf20_win$win))
windows_in[3,6] <- length(unique(snp_env5_peakbf20_win$win))
windows_in[3,7] <- length(unique(snp_env6_peakbf20_win$win))
windows_in[3,8] <- length(unique(snp_env7_peakbf20_win$win))
windows_in[3,9] <- length(unique(snp_env8_peakbf20_win$win))
windows_in[3,10] <- length(unique(snp_env9_peakbf20_win$win))

windows_in[4,1] <- ">10"
windows_in[4,2] <- length(unique(snp_env1_peakbf10_win$win))
windows_in[4,3] <- length(unique(snp_env2_peakbf10_win$win))
windows_in[4,4] <- length(unique(snp_env3_peakbf10_win$win))
windows_in[4,5] <- length(unique(snp_env4_peakbf10_win$win))
windows_in[4,6] <- length(unique(snp_env5_peakbf10_win$win))
windows_in[4,7] <- length(unique(snp_env6_peakbf10_win$win))
windows_in[4,8] <- length(unique(snp_env7_peakbf10_win$win))
windows_in[4,9] <- length(unique(snp_env8_peakbf10_win$win))
windows_in[4,10] <- length(unique(snp_env9_peakbf10_win$win))

windows_in[5,1] <- ">5"
windows_in[5,2] <- length(unique(snp_env1_peakbf5_win$win))
windows_in[5,3] <- length(unique(snp_env2_peakbf5_win$win))
windows_in[5,4] <- length(unique(snp_env3_peakbf5_win$win))
windows_in[5,5] <- length(unique(snp_env4_peakbf5_win$win))
windows_in[5,6] <- length(unique(snp_env5_peakbf5_win$win))
windows_in[5,7] <- length(unique(snp_env6_peakbf5_win$win))
windows_in[5,8] <- length(unique(snp_env7_peakbf5_win$win))
windows_in[5,9] <- length(unique(snp_env8_peakbf5_win$win))
windows_in[5,10] <- length(unique(snp_env9_peakbf5_win$win))

windows_in[6,1] <- ">2"
windows_in[6,2] <- length(unique(snp_env1_peakbf2_win$win))
windows_in[6,3] <- length(unique(snp_env2_peakbf2_win$win))
windows_in[6,4] <- length(unique(snp_env3_peakbf2_win$win))
windows_in[6,5] <- length(unique(snp_env4_peakbf2_win$win))
windows_in[6,6] <- length(unique(snp_env5_peakbf2_win$win))
windows_in[6,7] <- length(unique(snp_env6_peakbf2_win$win))
windows_in[6,8] <- length(unique(snp_env7_peakbf2_win$win))
windows_in[6,9] <- length(unique(snp_env8_peakbf2_win$win))
windows_in[6,10] <- length(unique(snp_env9_peakbf2_win$win))
colnames(windows_in) <- c("log10BF",
                          "env1 win",
                          "env2 win",
                          "env3 win",
                          "env4 win",
                          "env5 win",
                          "env6 win",
                          "env7 win",
                          "env8 win",
                          "env9 win")
  
windows_in

write_csv(snps_in_win,"data/genomic_data/snps_in_win_bf.csv")
write_csv(windows_in,"data/genomic_data/windows_in_bf.csv")


#############################################################################################################
#Merge BF>30 + WZA BF>10
#Setup bf10 snps for merger
bf10_env1 <- snp_env1_peakbf10_win %>% select(chr_snp,Env,BF) %>% filter(BF<30)
bf10_env2 <- snp_env2_peakbf10_win %>% select(chr_snp,Env,BF) %>% filter(BF<30)
bf10_env3 <- snp_env3_peakbf10_win %>% select(chr_snp,Env,BF) %>% filter(BF<30)
bf10_env4 <- snp_env4_peakbf10_win %>% select(chr_snp,Env,BF) %>% filter(BF<30)
bf10_env5 <- snp_env5_peakbf10_win %>% select(chr_snp,Env,BF) %>% filter(BF<30)
bf10_env6 <- snp_env6_peakbf10_win %>% select(chr_snp,Env,BF) %>% filter(BF<30)
bf10_env7 <- snp_env7_peakbf10_win %>% select(chr_snp,Env,BF) %>% filter(BF<30)
bf10_env8 <- snp_env8_peakbf10_win %>% select(chr_snp,Env,BF) %>% filter(BF<30)
bf10_env9 <- snp_env9_peakbf10_win %>% select(chr_snp,Env,BF) %>% filter(BF<30)


#Filter baseline by BF=>30
snp_set_env1 <- rbind(env1_united_bf30,bf10_env1)
snp_set_env2 <- rbind(env2_united_bf30,bf10_env2)
snp_set_env3 <- rbind(env3_united_bf30,bf10_env3)
snp_set_env4 <- rbind(env4_united_bf30,bf10_env4)
snp_set_env5 <- rbind(env5_united_bf30,bf10_env5)
snp_set_env6 <- rbind(env6_united_bf30,bf10_env6)
snp_set_env7 <- rbind(env7_united_bf30,bf10_env7)
snp_set_env8 <- rbind(env8_united_bf30,bf10_env8)
snp_set_env9 <- rbind(env9_united_bf30,bf10_env9)



#Make unique set
# BF>30 + WZA BF>10
snp_set_all <- rbind(snp_set_env1,
                     snp_set_env2,
                     snp_set_env3,
                     snp_set_env4,
                     snp_set_env5,
                     snp_set_env6,
                     snp_set_env7,
                     snp_set_env8,
                     snp_set_env9)
snp_set_env <- as.data.frame(unique(snp_set_all$chr_snp))
colnames(snp_set_env) <- "chr_snp"

#BF>30
snp_set_all_bf30 <- rbind(env1_united_bf30,
                          env2_united_bf30,
                          env3_united_bf30,
                          env4_united_bf30,
                          env5_united_bf30,
                          env6_united_bf30,
                          env7_united_bf30,
                          env8_united_bf30,
                          env9_united_bf30
                          )
snp_set_env_bf30 <- as.data.frame(unique(snp_set_all_bf30$chr_snp))
colnames(snp_set_env_bf30) <- "chr_snp"

#BF>20
snp_set_all_bf20 <- rbind(env1_united_bf20,
                          env2_united_bf20,
                          env3_united_bf20,
                          env4_united_bf20,
                          env5_united_bf20,
                          env6_united_bf20,
                          env7_united_bf20,
                          env8_united_bf20,
                          env9_united_bf20
)
snp_set_env_bf20 <- as.data.frame(unique(snp_set_all_bf20$chr_snp))
colnames(snp_set_env_bf20) <- "chr_snp"

#WZA BF>10
snp_set_all_wza10 <- rbind(snp_env1_peakbf10_win,
                           snp_env2_peakbf10_win,
                           snp_env3_peakbf10_win,
                           snp_env4_peakbf10_win,
                           snp_env5_peakbf10_win,
                           snp_env6_peakbf10_win,
                           snp_env7_peakbf10_win,
                           snp_env8_peakbf10_win,
                           snp_env9_peakbf10_win
)
snp_set_env_wza10 <- as.data.frame(unique(snp_set_all_wza10$chr_snp))
colnames(snp_set_env_wza10) <- "chr_snp"



#Write out WZA BayPass SNP set

#write_csv(snp_set_env1,"data/genomic_data/snp_set_env1.csv")
#write_csv(snp_set_env2,"data/genomic_data/snp_set_env2.csv")
#write_csv(snp_set_env3,"data/genomic_data/snp_set_env3.csv")
#write_csv(snp_set_env4,"data/genomic_data/snp_set_env4.csv")
#write_csv(snp_set_env5,"data/genomic_data/snp_set_env5.csv")
#write_csv(snp_set_env6,"data/genomic_data/snp_set_env6.csv")
#write_csv(snp_set_env7,"data/genomic_data/snp_set_env7.csv")
#write_csv(snp_set_env8,"data/genomic_data/snp_set_env8.csv")
#write_csv(snp_set_env9,"data/genomic_data/snp_set_env9.csv")
#write_csv(snp_set_env,"data/genomic_data/snp_set_env.csv")


#Write only SNP set for BF>30 only

write_csv(env1_united_bf30,"data/genomic_data/snp_set_bf30_env1.csv")
write_csv(env2_united_bf30,"data/genomic_data/snp_set_bf30_env2.csv")
write_csv(env3_united_bf30,"data/genomic_data/snp_set_bf30_env3.csv")
write_csv(env4_united_bf30,"data/genomic_data/snp_set_bf30_env4.csv")
write_csv(env5_united_bf30,"data/genomic_data/snp_set_bf30_env5.csv")
write_csv(env6_united_bf30,"data/genomic_data/snp_set_bf30_env6.csv")
write_csv(env7_united_bf30,"data/genomic_data/snp_set_bf30_env7.csv")
write_csv(env8_united_bf30,"data/genomic_data/snp_set_bf30_env8.csv")
write_csv(env9_united_bf30,"data/genomic_data/snp_set_bf30_env9.csv")
write_csv(snp_set_env_bf30,"data/genomic_data/snp_set_bf30.csv")

#Write only SNP set for BF>20 only

write_csv(env1_united_bf20,"data/genomic_data/snp_set_bf20_env1.csv")
write_csv(env2_united_bf20,"data/genomic_data/snp_set_bf20_env2.csv")
write_csv(env3_united_bf20,"data/genomic_data/snp_set_bf20_env3.csv")
write_csv(env4_united_bf20,"data/genomic_data/snp_set_bf20_env4.csv")
write_csv(env5_united_bf20,"data/genomic_data/snp_set_bf20_env5.csv")
write_csv(env6_united_bf20,"data/genomic_data/snp_set_bf20_env6.csv")
write_csv(env7_united_bf20,"data/genomic_data/snp_set_bf20_env7.csv")
write_csv(env8_united_bf20,"data/genomic_data/snp_set_bf20_env8.csv")
write_csv(env9_united_bf20,"data/genomic_data/snp_set_bf20_env9.csv")
write_csv(snp_set_env_bf20,"data/genomic_data/snp_set_bf20.csv")

#Write out for WZA BF>10 only

write_csv(snp_env1_peakbf10_win,"data/genomic_data/snp_set_wza10_env1.csv")
write_csv(snp_env2_peakbf10_win,"data/genomic_data/snp_set_wza10_env2.csv")
write_csv(snp_env3_peakbf10_win,"data/genomic_data/snp_set_wza10_env3.csv")
write_csv(snp_env4_peakbf10_win,"data/genomic_data/snp_set_wza10_env4.csv")
write_csv(snp_env5_peakbf10_win,"data/genomic_data/snp_set_wza10_env5.csv")
write_csv(snp_env6_peakbf10_win,"data/genomic_data/snp_set_wza10_env6.csv")
write_csv(snp_env7_peakbf10_win,"data/genomic_data/snp_set_wza10_env7.csv")
write_csv(snp_env8_peakbf10_win,"data/genomic_data/snp_set_wza10_env8.csv")
write_csv(snp_env9_peakbf10_win,"data/genomic_data/snp_set_wza10_env9.csv")
write_csv(snp_set_env_wza10,"data/genomic_data/snp_set_wza10.csv")



