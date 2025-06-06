#############################################################################################################
## Calculate windows for WZA
## Author Daniel Anstett
## 
## Modified from Tom Booker WZA Vignette
## Last Modified April 16, 2024
#############################################################################################################
#Functions
#Calculate frequency
prop_A <- function(snp_table) {
  snp_prop_A<- snp_table %>% select (chr_snp)
  counter=2
  pop_num=1
  for (i in seq(1,dim(snp_table)[2]-1,2)) {
    
    snpA<-paste("V",i, sep="") #sets up string for snpA for pop num
    snpB<-paste("V",i+1, sep="") #sets up string for snpA for pop num
    P<-paste("P", pop_num, sep="") #sets up paper_ID (population ID)
    tmp<-snp_table %>% select(snpA,snpB) #graps all snps per paper_ID
    
    colnames(tmp)<-c("A", "B") #renames column headers
    
    snp_prop_A[,counter]<-tmp$A/(tmp$A + tmp$B) #calc proportion for A, assign to output
    colnames (snp_prop_A)[counter]<-P 
    
    counter<-counter+1
    pop_num<-pop_num+1
  }
  #snp_prop_A[is.na(snp_prop_A)] <- 0
  return(snp_prop_A)
}

#############################################################################################################
#Import libraries
library(tidyverse)
library(Kendall)

#Import full snp table for baseline
pop_order<-read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", header=F, sep="\t")
snp<-read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", header=F, sep=" ")
loci<-read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", header=F, sep="\t")
loci_win<-read_csv("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/loci_win.csv")

#Calculate frequency using prop_A
colnames(loci) <- c("chr","snp")
loci_united <- loci %>% unite(chr_snp,"chr","snp",sep="_")
snp_chr_snp <- cbind(loci_united,snp)
freq <- prop_A(snp_chr_snp)
all_data <-cbind(loci_win,freq) %>% mutate(chr_snp=paste(chr, snp, sep="_")) #add snp lables to rows
colnames(all_data)[4] <- "win"


#############################################################################################################
#############################################################################################################
#WZA data prep

# The first thing we'll do is calculate the average allele frequency across populations for each SNP and add it to the SNP info dataframe
#Calc row means
all_data<-all_data %>% mutate(p_bar = rowMeans(select(all_data, P1:P55), na.rm = TRUE))
all_data$q_bar <- 1 - all_data$p_bar

# take the minimum to get the MAF
all_data$MAF <- pmin(all_data$p_bar, all_data$q_bar)


#############################################################################################################

#Import BayPass Results
env1 <- read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/ENV_1_trim.tsv",header=F, sep="\t")
env2 <- read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/ENV_2_trim.tsv",header=F, sep="\t")
env3 <- read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/ENV_3_trim.tsv",header=F, sep="\t")
env4 <- read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/ENV_4_trim.tsv",header=F, sep="\t")
env5 <- read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/ENV_5_trim.tsv",header=F, sep="\t")

env6 <- read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/ENV_6_trim.tsv",header=F, sep="\t")
env7 <- read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/ENV_7_trim.tsv",header=F, sep="\t")
env8 <- read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/ENV_8_trim.tsv",header=F, sep="\t")
env9 <- read.table("/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/ENV_9_trim.tsv",header=F, sep="\t")

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

env1_united <- env1 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env2_united <- env2 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env3_united <- env3 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env4_united <- env4 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env5_united <- env5 %>% unite(chr_snp,"Chromosome","SNP",sep="_")

env6_united <- env6 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env7_united <- env7 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env8_united <- env8 %>% unite(chr_snp,"Chromosome","SNP",sep="_")
env9_united <- env9 %>% unite(chr_snp,"Chromosome","SNP",sep="_")

snps_env1_bf <- left_join(all_data,env1_united, by="chr_snp")
snps_env2_bf <- left_join(all_data,env2_united, by="chr_snp")
snps_env3_bf <- left_join(all_data,env3_united, by="chr_snp")
snps_env4_bf <- left_join(all_data,env4_united, by="chr_snp")
snps_env5_bf <- left_join(all_data,env5_united, by="chr_snp")

snps_env6_bf <- left_join(all_data,env6_united, by="chr_snp")
snps_env7_bf <- left_join(all_data,env7_united, by="chr_snp")
snps_env8_bf <- left_join(all_data,env8_united, by="chr_snp")
snps_env9_bf <- left_join(all_data,env9_united, by="chr_snp")

#Too large to store on github. Store locally
write_csv(snps_env1_bf, "/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/WZA_snps_env1_bf.csv")
write_csv(snps_env2_bf, "/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/WZA_snps_env2_bf.csv")
write_csv(snps_env3_bf, "/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/WZA_snps_env3_bf.csv")
write_csv(snps_env4_bf, "/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/WZA_snps_env4_bf.csv")
write_csv(snps_env5_bf, "/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/WZA_snps_env5_bf.csv")     

write_csv(snps_env6_bf, "/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/WZA_snps_env6_bf.csv")
write_csv(snps_env7_bf, "/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/WZA_snps_env7_bf.csv")
write_csv(snps_env8_bf, "/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/WZA_snps_env8_bf.csv")
write_csv(snps_env9_bf, "/Users/daniel_anstett/Dropbox/z_Documents/aLarge_files/M_gen/WZA_snps_env9_bf.csv")

